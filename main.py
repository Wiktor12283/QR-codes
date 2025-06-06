from versions import get_versions
from ECCCW import error_correction_codewords
from alignment_patterns import alignment_pattern_locations
from itertools import product

version = 1
qr_prompt = list(input())
dec_prompt = [ord(char) for char in qr_prompt]
bit_prompt = ["0" + bin(dec_value)[2:] for dec_value in dec_prompt]

for ver, ecc in get_versions().items():
    if len(bit_prompt) < ecc["L"]:
        version = ver
        break


def generate_galois_field(data):
    gf = []
    for bits in data:
        gf.append(int(bits, 2))
    return gf


galois_field = generate_galois_field(bit_prompt)

# Log and Antilog tables
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


# Implementing galois field arythmetic
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


# Polynomial division
def rs_encode_msg(msg, nsym):
    gen = generate_generator_poly(nsym)
    msg_out = msg + [0] * nsym

    for i in range(len(msg)):
        coef = msg_out[i]
        if coef != 0:
            for j in range(len(gen)):
                msg_out[i + j] ^= gf_mul(gen[j], coef)

    # msg_out now holds the error correction codewords at the end
    return msg + msg_out[len(msg) :]


# Data analysis
num_eccw = error_correction_codewords[version]["L"]
encoded_message = rs_encode_msg(galois_field, num_eccw)
data_part = encoded_message[:-num_eccw]
ec_part = encoded_message[-num_eccw:]

# create qr code
matrix = [
    ["#" for _ in range((((version - 1) * 4) + 21))]
    for _ in range((((version - 1) * 4) + 21))
]


# Finder patterns
def add_hollow_square(matrix, top_left_x, top_left_y, size, value):
    for i in range(size):
        for j in range(size):
            if i == 0 or i == size - 1 or j == 0 or j == size - 1:
                matrix[top_left_y + i][top_left_x + j] = value


add_hollow_square(matrix, 0, 0, 7, 1)
add_hollow_square(matrix, len(matrix) - 7, 0, 7, 1)
add_hollow_square(matrix, 0, len(matrix) - 7, 7, 1)
add_hollow_square(matrix, 1, 1, 5, 0)
add_hollow_square(matrix, len(matrix) - 6, 1, 5, 0)
add_hollow_square(matrix, 1, len(matrix) - 6, 5, 0)
add_hollow_square(matrix, 2, 2, 3, 1)
matrix[3][3] = 1
add_hollow_square(matrix, len(matrix) - 5, 2, 3, 1)
matrix[3][-4] = 1
add_hollow_square(matrix, 2, len(matrix) - 5, 3, 1)
matrix[-4][3] = 1

for i in range(8):
    matrix[7][i] = 0
    matrix[i][7] = 0
    matrix[7][-i - 1] = 0
    matrix[i][-8] = 0
    matrix[-8][i] = 0
    matrix[-i - 1][7] = 0

# Timing patterns
for i in range(8, len(matrix) - 8):
    matrix[6][i] = 1 if i % 2 == 0 else 0
    matrix[i][6] = 1 if i % 2 == 0 else 0

# Alignment patterns

alignment_patterns = list(product(alignment_pattern_locations[version], repeat=2))

print("Alignment Patterns:", alignment_patterns)


def matrix_insert_alignment_patterns(matrix, alignment_patterns):
    for x, y in alignment_patterns:
        # Skip alignment patterns overlapping with finder patterns
        if (
            (x < 9 and y < 9)
            or (x < 9 and y > len(matrix) - 9)
            or (x > len(matrix) - 9 and y < 9)
        ):
            continue
        # Add alignment pattern
        add_hollow_square(matrix, x - 2, y - 2, 5, 1)
        add_hollow_square(matrix, x - 1, y - 1, 3, 0)
        matrix[y][x] = 1


matrix_insert_alignment_patterns(matrix, alignment_patterns)

# Adding dark module
matrix[(4 * version) + 9][8] = 1


# Placing data and error correction codewords
def place_data(matrix, data_bits):
    size = len(matrix)
    bit_idx = 0
    x = size - 1
    direction = -1  # Initial direction is up (-1)

    while x > 0:
        if x == 6:  # Skip vertical timing pattern
            x -= 1

        y_range = range(size - 1, -1, -1) if direction == -1 else range(size)

        for y in y_range:
            for dx in [0, -1]:  # Right column first, then left
                col = x + dx
                if 0 <= col < size and 0 <= y < size:
                    if matrix[y][col] is None:
                        if bit_idx < len(data_bits):
                            matrix[y][col] = int(data_bits[bit_idx])
                            bit_idx += 1
                        else:
                            matrix[y][col] = 0  # pad if no bits left

        x -= 2
        direction *= -1  # Reverse direction

    return matrix


matrix = place_data(matrix, data_part)

for line in matrix:
    print("".join(str(cell) for cell in line))
