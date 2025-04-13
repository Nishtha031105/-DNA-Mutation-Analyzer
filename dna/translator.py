complement = {
    'A': 'T',
    'T': 'A',
    'C': 'G',
    'G': 'C'
}

def get_complementary_sequence(input_dna):
    complementary_sequence = ""

    for base in input_dna:
        complementary_base = complement.get(base, base)  
        complementary_sequence += complementary_base

    return complementary_sequence