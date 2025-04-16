import matplotlib.pyplot as plt
import numpy as np

def calculate_similarity(seq1, seq2):
    """Calculate percentage similarity between two sequences"""
    # Ensure sequences are the same length for direct comparison
    min_len = min(len(seq1), len(seq2))
    seq1 = seq1[:min_len]
    seq2 = seq2[:min_len]
    
    # Count matching positions
    matches = sum(n1 == n2 for n1, n2 in zip(seq1, seq2))
    
    # Calculate percentage
    similarity = (matches / min_len) * 100
    return similarity

def sliding_window_similarity(human_seq, other_seq, window_size=50):
    """Calculate similarity across sliding windows"""
    similarities = []
    min_len = min(len(human_seq), len(other_seq))
    
    for i in range(0, min_len - window_size + 1):
        human_window = human_seq[i:i+window_size]
        other_window = other_seq[i:i+window_size]
        sim = calculate_similarity(human_window, other_window)
        similarities.append(sim)
    
    return similarities

def find_diagnostic_sites(human_seq, other_seq, species_name, max_sites=10):
    """Find species-specific diagnostic sites"""
    min_len = min(len(human_seq), len(other_seq))
    
    print(f"\nDiagnostic Sites where {species_name} differs from human:")
    print("=" * 50)
    
    differences = []
    for pos in range(min_len):
        if human_seq[pos] != other_seq[pos]:
            differences.append(f"Position {pos+1}: Human={human_seq[pos]}, {species_name}={other_seq[pos]}")
            if len(differences) >= max_sites:
                differences.append(f"... and more differences")
                break
    
    for diff in differences:
        print(f"  {diff}")

def clean_sequence(seq):
    """Clean a DNA sequence by removing whitespace and normalizing to uppercase"""
    return ''.join(s.upper() for s in seq if s.upper() in 'ATGC')

def main():
    # Human COI reference sequence
    human_coi = """
    GTCCTACTATCCATGCAGGTATCTTCTATCTTTGGGGCATGAGCGGGCATAGTAGGCACAGCCCTAAGCCTCCTCATTCG
    AGCCGAGCTGGGCCAGCCAGGCAACCTTCTAGGTAACGACCACATCTACAACGTTATCGTCACAGCCCATGCATTCGTAA
    TAATCTTCTTCATAGTAATACCAATAATAATCGGAGGCTTCGGAAACTGACTAGTCCCCCTTATAATCGGTGCCCCCGAC
    ATAGCATTCCCACGAATAAATAACATAAGCTTCTGACTCCTCCCTCCATCCTTTCTCCTTCTTCTCGCATCCTCCGGAGT
    AGAAGCTGGCGCAGGAACAGGCTGAACAGTCTACCCTCCCCTAGCAGGAAACTACTCCCACCCTGGAGCCTCCGTAGACC
    TGGCAATCTTCTCCCTCCACCTGGCAGGTATTTCCTCCATCCTAGGAGCAATTAACTTCATCACCACAGCTATCAACATA
    AAACCCCCTGCAATATCCCAGTATCAAACTCCCCTATTCGTCTGATCAGTCCTAATTACCGCCGTCCTACTCCTCCTGTC
    CCTGCCCGTCCTCGCTGCAGGAATCACAATGCTGCTCACAGACCGCAACCTTAACACCACCTTCTTCGACCCGGCAGGAG
    GAGGAGACCCAGTCCTGTACCAACACCTATTCT
    """
    human_coi = clean_sequence(human_coi)
    
    print("DNA Sequence Comparison Tool")
    print("==========================")
    print("\nThis tool compares a DNA sequence from your species to human DNA.")
    
    # Get species name
    species_name = input("\nEnter the name of your species: ")
    
    # Get sequence input
    print(f"\nEnter the DNA sequence for {species_name}:")
    print("(A, T, G, C only - other characters will be removed)")
    species_seq = input()
    species_seq = clean_sequence(species_seq)
    
    # Check if sequence is valid
    if not species_seq:
        print("Error: The sequence provided is empty or contains no valid DNA bases (A, T, G, C).")
        return
    
    # Calculate overall similarity
    similarity = calculate_similarity(human_coi, species_seq)
    print(f"\nOverall DNA Similarity: {similarity:.2f}%")
    
    # Find diagnostic sites
    find_diagnostic_sites(human_coi, species_seq, species_name)
    
    # Create sliding window similarity plot
    window_similarities = sliding_window_similarity(human_coi, species_seq)
    
    plt.figure(figsize=(10, 6))
    plt.plot(window_similarities)
    plt.title(f"DNA Similarity: Human vs {species_name}")
    plt.xlabel("Position (50 bp window)")
    plt.ylabel("Percent Similarity")
    plt.ylim(0, 100)
    plt.axhline(y=90, color='g', linestyle='--', alpha=0.3, label="High conservation (90%)")
    plt.axhline(y=50, color='r', linestyle='--', alpha=0.3, label="Low conservation (50%)")
    plt.legend()
    
    try:
        plt.savefig("dna_comparison.png")
        print("\nPlot saved as 'dna_comparison.png'")
    except:
        print("\nCouldn't save plot file")
    
    try:
        plt.show()
    except:
        print("Couldn't display plot - you may need to run this in an environment that supports matplotlib visualization")

if __name__ == "__main__":
    main()
