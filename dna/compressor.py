base_to_bits = {
    'A': '00',
    'C': '01',
    'G': '10',
    'T': '11'
}
bits_to_base = {v: k for k, v in base_to_bits.items()}


def compress(dna_str: str) -> bytearray:
    
    binary_string = ''.join([base_to_bits[base] for base in dna_str])

    while len(binary_string) % 8 != 0:
        binary_string += '00'

    compressed = bytearray()
    for i in range(0, len(binary_string), 8):
        byte = binary_string[i:i+8]
        compressed.append(int(byte, 2))

    return compressed


def decompress(compressed_data: bytearray, original_length: int) -> str:
    binary_string = ''
    for byte in compressed_data:
        binary_string += f'{byte:08b}'

    dna_str = ''
    for i in range(0, original_length * 2, 2):
        bits = binary_string[i:i+2]
        dna_str += bits_to_base.get(bits, '?')  

    return dna_str
