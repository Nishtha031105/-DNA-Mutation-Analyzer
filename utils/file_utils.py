import os
def read_dna_file(filepath: str) -> str:
    
    if not os.path.exists(filepath):
        raise FileNotFoundError(f"File '{filepath}' not found.")

    with open(filepath, 'r') as file:
        dna = file.read()

    input_dna = ''.join(dna.split()).upper()

    for char in input_dna:
        if char not in ['A', 'T', 'C', 'G']:
            raise ValueError(f"Invalid DNA character found: '{char}'")

    return input_dna

def write_dna_file(filename: str, data: str):
    with open(filename, "w") as f:
        f.write(data)
