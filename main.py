from versions import get_versions
from ECCCW import error_correction_codewords
from alignment_patterns import alignment_pattern_locations
from itertools import product

# Get input and prepare data
qr_prompt = list(input())
mode_indicator = [0, 1, 0, 0]  # Byte mode
char_count_bits = [int(b) for b in format(len(qr_prompt), '08b')]
data_bits = []
for char in qr_prompt:
    data_bits.extend([int(b) for b in format(ord(char), '08b')])

# Combine: mode + count + data
bit_prompt = mode_indicator + char_count_bits + data_bits

# Fix version selection logic
version = 1
for ver, ecc in get_versions().items():
    if len(qr_prompt) <= ecc["L"]:  # Fixed: should be <= not <
        version = ver
        break

# Galois Field setup
GF_POLY = 0x11D
GF_SIZE = 256

# Precompute log and antilog tables
exp_table = [0] * (GF_SIZE * 2)
log_table = [0] * GF_SIZE

x = 1
for i in range(GF_SIZE - 1):
    exp_table[i] = x
    log_table[x] = i
    x <<= 1
    if x & 0x100:
        x ^= GF_POLY

# Extend exp_table to avoid overflow
for i in range(GF_SIZE - 1, GF_SIZE * 2):
    exp_table[i] = exp_table[i - (GF_SIZE - 1)]

# Galois field arithmetic
def gf_add(x, y):
    return x ^ y

def gf_sub(x, y):
    return x ^ y  # same as addition in GF(256)

def gf_mul(x, y):
    if x == 0 or y == 0:
        return 0
    return exp_table[log_table[x] + log_table[y]]

def gf_div(x, y):
    if y == 0:
        raise ZeroDivisionError()
    if x == 0:
        return 0
    return exp_table[(log_table[x] - log_table[y]) % 255]

# Generator polynomials
def generate_generator_poly(nsym):
    g = [1]
    for i in range(nsym):
        g = poly_mul(g, [1, exp_table[i]])
    return g

def poly_mul(p, q):
    res = [0] * (len(p) + len(q) - 1)
    for j in range(len(q)):
        for i in range(len(p)):
            res[i + j] ^= gf_mul(p[i], q[j])
    return res

# Reed-Solomon encoding
def rs_encode_msg(msg, nsym):
    gen = generate_generator_poly(nsym)
    msg_out = msg + [0] * nsym

    for i in range(len(msg)):
        coef = msg_out[i]
        if coef != 0:
            for j in range(len(gen)):
                msg_out[i + j] ^= gf_mul(gen[j], coef)

    return msg + msg_out[len(msg):]

# Convert bits to bytes and apply error correction
num_eccw = error_correction_codewords[version]["L"]
data_bytes = []
for i in range(0, len(bit_prompt), 8):
    byte_bits = bit_prompt[i:i+8]
    if len(byte_bits) == 8:
        data_bytes.append(int(''.join(map(str, byte_bits)), 2))

# Apply error correction
encoded_message = rs_encode_msg(data_bytes, num_eccw)

# Create QR code matrix
size = (version - 1) * 4 + 21
matrix = [[0 for _ in range(size)] for _ in range(size)]  # Initialize with 0s, not "#"

# Helper function for adding squares
def add_hollow_square(matrix, start_x, start_y, square_size, value):
    # Top and bottom borders
    for x in range(start_x, start_x + square_size):
        if 0 <= start_y < len(matrix) and 0 <= x < len(matrix[0]):
            matrix[start_y][x] = value
        if 0 <= start_y + square_size - 1 < len(matrix) and 0 <= x < len(matrix[0]):
            matrix[start_y + square_size - 1][x] = value
    
    # Left and right borders
    for y in range(start_y, start_y + square_size):
        if 0 <= y < len(matrix) and 0 <= start_x < len(matrix[0]):
            matrix[y][start_x] = value
        if 0 <= y < len(matrix) and 0 <= start_x + square_size - 1 < len(matrix[0]):
            matrix[y][start_x + square_size - 1] = value

# Add finder patterns
def add_finder_patterns(matrix):
    size = len(matrix)
    
    # Three finder patterns
    positions = [(0, 0), (size - 7, 0), (0, size - 7)]
    
    for x, y in positions:
        # 7x7 outer dark border
        add_hollow_square(matrix, x, y, 7, 1)
        # 5x5 inner light square
        add_hollow_square(matrix, x + 1, y + 1, 5, 0)
        # 3x3 inner dark square
        add_hollow_square(matrix, x + 2, y + 2, 3, 1)
        # Center dark module
        matrix[y + 3][x + 3] = 1

    # Add separators (white borders around finder patterns)
    for i in range(8):
        # Top-left
        matrix[7][i] = 0
        matrix[i][7] = 0
        # Top-right
        matrix[7][size - i - 1] = 0
        matrix[i][size - 8] = 0
        # Bottom-left
        matrix[size - 8][i] = 0
        matrix[size - i - 1][7] = 0

add_finder_patterns(matrix)

# Add timing patterns
for i in range(8, size - 8):
    matrix[6][i] = 1 if i % 2 == 0 else 0
    matrix[i][6] = 1 if i % 2 == 0 else 0

# Add alignment patterns
def add_alignment_patterns(matrix, version):
    if version == 1:
        return  # No alignment patterns for version 1
    
    locations = alignment_pattern_locations[version]
    alignment_positions = list(product(locations, repeat=2))
    
    for x, y in alignment_positions:
        # Skip if overlapping with finder patterns
        if ((x < 9 and y < 9) or 
            (x < 9 and y > size - 9) or 
            (x > size - 9 and y < 9)):
            continue
        
        # Add 5x5 alignment pattern
        add_hollow_square(matrix, x - 2, y - 2, 5, 1)
        add_hollow_square(matrix, x - 1, y - 1, 3, 0)
        matrix[y][x] = 1

add_alignment_patterns(matrix, version)

# Add dark module
matrix[(4 * version) + 9][8] = 1

# Data placement functions
def is_function_module(matrix, row, col, version):
    """Check if a module is a function module (not for data)"""
    size = len(matrix)
    
    # Check bounds
    if row < 0 or row >= size or col < 0 or col >= size:
        return True
    
    # Finder patterns (8x8 including separators)
    if ((row <= 8 and col <= 8) or  # Top-left
        (row <= 8 and col >= size - 8) or  # Top-right  
        (row >= size - 8 and col <= 8)):  # Bottom-left
        return True
    
    # Timing patterns
    if row == 6 or col == 6:
        return True
    
    # Dark module
    if row == ((4 * version) + 9) and col == 8:
        return True
    
    # Format information areas
    if ((row == 8 and (col <= 8 or col >= size - 8)) or
        (col == 8 and (row <= 8 or row >= size - 7))):
        return True

    # Version information (for versions 7+)
    if version >= 7:
        if ((size - 11 <= row <= size - 9) and (0 <= col <= 5)) or \
           ((0 <= row <= 5) and (size - 11 <= col <= size - 9)):
            return True
    
    # Check alignment patterns
    if version > 1:
        locations = alignment_pattern_locations[version]
        for ax, ay in product(locations, repeat=2):
            # Skip alignment patterns that overlap with finder patterns
            if ((ax < 9 and ay < 9) or 
                (ax < 9 and ay > size - 9) or 
                (ax > size - 9 and ay < 9)):
                continue
            
            # Check if current position is within alignment pattern
            if abs(row - ay) <= 2 and abs(col - ax) <= 2:
                return True
    
    return False

def add_padding_to_data(data_bits, total_capacity):
    """Add padding to data to fill the QR code capacity"""
    padded_bits = data_bits[:]
    
    # Step 1: Add up to 4 terminator bits (0000) if there's room
    terminator_bits = min(4, total_capacity - len(padded_bits))
    padded_bits.extend([0] * terminator_bits)
    
    # Step 2: Pad to make length multiple of 8 (byte boundary)
    while len(padded_bits) % 8 != 0 and len(padded_bits) < total_capacity:
        padded_bits.append(0)
    
    # Step 3: Add alternating pad bytes (236, 17) until capacity is reached
    pad_bytes = [236, 17]  # 11101100, 00010001 in binary
    pad_index = 0
    
    while len(padded_bits) < total_capacity:
        pad_byte = pad_bytes[pad_index % 2]
        pad_bits = [int(bit) for bit in format(pad_byte, '08b')]
        
        remaining_capacity = total_capacity - len(padded_bits)
        bits_to_add = min(8, remaining_capacity)
        padded_bits.extend(pad_bits[:bits_to_add])
        
        pad_index += 1
    
    return padded_bits

def place_data_in_matrix(matrix, data_bits, version):
    """Place data bits in the QR code matrix"""
    size = len(matrix)
    
    # Calculate total data capacity
    total_capacity = 0
    for row in range(size):
        for col in range(size):
            if not is_function_module(matrix, row, col, version):
                total_capacity += 1
    
    # Pad data to fill capacity
    if len(data_bits) < total_capacity:
        data_bits = add_padding_to_data(data_bits, total_capacity)
    
    bit_index = 0
    up = True  # Direction flag: True for up, False for down
    col = size - 1
    
    while col > 0 and bit_index < len(data_bits):
        # Skip timing pattern column
        if col == 6:
            col -= 1
            
        # Process current column pair
        for i in range(size):
            row = (size - 1 - i) if up else i
            
            # Process right column of the pair
            if bit_index < len(data_bits) and not is_function_module(matrix, row, col, version):
                matrix[row][col] = data_bits[bit_index]
                bit_index += 1
                
            # Process left column of the pair
            if bit_index < len(data_bits) and col > 0 and not is_function_module(matrix, row, col - 1, version):
                matrix[row][col - 1] = data_bits[bit_index]
                bit_index += 1
        
        col -= 2
        up = not up

# Convert encoded message to bits
def byte_list_to_bit_list(byte_list):
    bits = []
    for byte in byte_list:
        bits.extend([int(bit) for bit in format(byte, '08b')])
    return bits

final_bits = byte_list_to_bit_list(encoded_message)

# Place data in matrix
place_data_in_matrix(matrix, final_bits, version)

# Masking functions
def get_mask_pattern(mask_number, row, col):
    """Get mask pattern value for given mask number and position"""
    if mask_number == 0:
        return (row + col) % 2 == 0
    elif mask_number == 1:
        return row % 2 == 0
    elif mask_number == 2:
        return col % 3 == 0
    elif mask_number == 3:
        return (row + col) % 3 == 0
    elif mask_number == 4:
        return (row // 2 + col // 3) % 2 == 0
    elif mask_number == 5:
        return ((row * col) % 2) + ((row * col) % 3) == 0
    elif mask_number == 6:
        return (((row * col) % 2) + ((row * col) % 3)) % 2 == 0
    elif mask_number == 7:
        return (((row + col) % 2) + ((row * col) % 3)) % 2 == 0
    else:
        raise ValueError("Mask pattern must be 0-7")

def apply_mask(matrix, mask_number, version):
    """Apply mask pattern to data modules only"""
    size = len(matrix)
    masked_matrix = [row[:] for row in matrix]  # Deep copy
    
    for row in range(size):
        for col in range(size):
            # Only apply mask to data modules
            if not is_function_module(matrix, row, col, version):
                if get_mask_pattern(mask_number, row, col):
                    masked_matrix[row][col] = 1 - masked_matrix[row][col]
    
    return masked_matrix

def calculate_mask_penalty(matrix):
    """Calculate penalty score for mask pattern"""
    size = len(matrix)
    penalty = 0
    
    # Rule 1: Adjacent modules in row/column of same color
    for row in range(size):
        count = 1
        for col in range(1, size):
            if matrix[row][col] == matrix[row][col-1]:
                count += 1
            else:
                if count >= 5:
                    penalty += 3 + (count - 5)
                count = 1
        if count >= 5:
            penalty += 3 + (count - 5)
    
    for col in range(size):
        count = 1
        for row in range(1, size):
            if matrix[row][col] == matrix[row-1][col]:
                count += 1
            else:
                if count >= 5:
                    penalty += 3 + (count - 5)
                count = 1
        if count >= 5:
            penalty += 3 + (count - 5)
    
    # Rule 2: Block of modules of same color (2x2 or larger)
    for row in range(size - 1):
        for col in range(size - 1):
            if (matrix[row][col] == matrix[row][col + 1] == 
                matrix[row + 1][col] == matrix[row + 1][col + 1]):
                penalty += 3
    
    # Rule 3: 1:1:3:1:1 ratio pattern in row/column
    patterns = [
        [1, 0, 1, 1, 1, 0, 1, 0, 0, 0, 0],  # 10111010000
        [0, 0, 0, 0, 1, 0, 1, 1, 1, 0, 1]   # 00001011101
    ]
    
    for row in range(size):
        for col in range(size - 10):
            for pattern in patterns:
                if all(matrix[row][col + i] == pattern[i] for i in range(11)):
                    penalty += 40
    
    for col in range(size):
        for row in range(size - 10):
            for pattern in patterns:
                if all(matrix[row + i][col] == pattern[i] for i in range(11)):
                    penalty += 40
    
    # Rule 4: Proportion of dark modules
    total_modules = size * size
    dark_modules = sum(sum(row) for row in matrix)
    dark_ratio = dark_modules / total_modules
    
    ratio_percent = int(dark_ratio * 100)
    lower_multiple = (ratio_percent // 5) * 5
    upper_multiple = lower_multiple + 5
    
    penalty += min(abs(ratio_percent - lower_multiple), 
                  abs(ratio_percent - upper_multiple)) * 10
    
    return penalty

def find_best_mask(matrix, version):
    """Find the mask pattern with lowest penalty"""
    best_mask = 0
    best_penalty = float('inf')
    best_matrix = None
    
    for mask_num in range(8):
        masked_matrix = apply_mask(matrix, mask_num, version)
        penalty = calculate_mask_penalty(masked_matrix)
        
        if penalty < best_penalty:
            best_penalty = penalty
            best_mask = mask_num
            best_matrix = masked_matrix
    
    return best_mask, best_matrix

# Format information functions
def calculate_bch_code(data, generator_poly):
    """Calculate BCH error correction code"""
    data_poly = data << (generator_poly.bit_length() - 1)
    
    while data_poly.bit_length() >= generator_poly.bit_length():
        data_poly ^= generator_poly << (data_poly.bit_length() - generator_poly.bit_length())
    
    return data_poly

def generate_format_info(error_correction_level, mask_pattern):
    """Generate 15-bit format information"""
    ec_bits = {'L': 0b01, 'M': 0b00, 'Q': 0b11, 'H': 0b10}
    
    if isinstance(error_correction_level, str):
        ec_level = ec_bits[error_correction_level.upper()]
    else:
        ec_level = error_correction_level
    
    format_data = (ec_level << 3) | mask_pattern
    generator_poly = 0b10100110111  # BCH(15,5) generator
    bch_code = calculate_bch_code(format_data, generator_poly)
    
    format_info = (format_data << 10) | bch_code
    format_mask = 0b101010000010010
    format_info ^= format_mask
    
    return format_info

def place_format_info(matrix, format_info, version):
    """Place format information in the QR code"""
    size = len(matrix)
    format_bits = [(format_info >> i) & 1 for i in range(14, -1, -1)]
    
    # First copy around top-left finder pattern
    bit_index = 0
    
    # Horizontal strip in row 8
    positions_row8 = [0, 1, 2, 3, 4, 5, 7, 8]
    for col in positions_row8:
        if bit_index < 15:
            matrix[8][col] = format_bits[bit_index]
            bit_index += 1
    
    # Vertical strip in column 8
    positions_col8 = [7, 5, 4, 3, 2, 1, 0]
    for row in positions_col8:
        if bit_index < 15:
            matrix[row][8] = format_bits[bit_index]
            bit_index += 1
    
    # Second copy split between top-right and bottom-left
    bit_index = 0
    
    # Top-right horizontal
    for i in range(8):
        if bit_index < 15:
            matrix[8][size - 1 - i] = format_bits[bit_index]
            bit_index += 1
    
    # Bottom-left vertical
    for i in range(7):
        if bit_index < 15:
            matrix[size - 7 + i][8] = format_bits[bit_index]
            bit_index += 1

# Apply masking and add format info
best_mask, masked_matrix = find_best_mask(matrix, version)

# Add format information (after masking, but format info itself is not masked)
format_info = generate_format_info('L', best_mask)
place_format_info(masked_matrix, format_info, version)

# Display function
def print_qr_pretty(matrix):
    """Print QR code with proper quiet zone"""
    quiet_zone = 4
    size = len(matrix)
    empty_block = "\u2588" * 2  # Dark block
    full_block = "  "  # Light block

    # Top quiet zone
    for _ in range(quiet_zone):
        print(empty_block * (size + 2 * quiet_zone))

    for row in matrix:
        line = empty_block * quiet_zone  # Left quiet zone
        for cell in row:
            if cell == 1:
                line += full_block  # Light for 1
            else:  # 0 or any other value
                line += empty_block  # Dark for 0
        line += empty_block * quiet_zone  # Right quiet zone
        print(line)

    # Bottom quiet zone
    for _ in range(quiet_zone):
        print(empty_block * (size + 2 * quiet_zone))

# Print the final QR code
print_qr_pretty(masked_matrix)