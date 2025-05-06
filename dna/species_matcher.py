import matplotlib.pyplot as plt
import numpy as np

def calculate_similarity(seq1, seq2):
    """Calculate percentage similarity between two sequences"""
    min_len = min(len(seq1), len(seq2))
    seq1 = seq1[:min_len]
    seq2 = seq2[:min_len]
    
    matches = sum(n1 == n2 for n1, n2 in zip(seq1, seq2))
    
    similarity = (matches / min_len) * 100
    return similarity

def dp_sequence_alignment(seq1, seq2):
    """
    Use dynamic programming to find optimal sequence alignment
    Returns similarity score and aligned sequences
    """
    match_score = 1
    mismatch_penalty = -1
    gap_penalty = -2
    
    m, n = len(seq1), len(seq2)
    dp = [[0 for _ in range(n+1)] for _ in range(m+1)]
    
    for i in range(m+1):
        dp[i][0] = i * gap_penalty
    for j in range(n+1):
        dp[0][j] = j * gap_penalty
    
    for i in range(1, m+1):
        for j in range(1, n+1):
            match = dp[i-1][j-1] + (match_score if seq1[i-1] == seq2[j-1] else mismatch_penalty)
            delete = dp[i-1][j] + gap_penalty
            insert = dp[i][j-1] + gap_penalty
            
            # Take the maximum score
            dp[i][j] = max(match, delete, insert)
    
    max_possible_score = max(m, n) * match_score
    similarity = (dp[m][n] / max_possible_score) * 100
    
    similarity = max(0, min(100, similarity))
    
    return similarity, dp

def sliding_window_similarity(human_seq, other_seq, window_size=50):
    """Calculate similarity across sliding windows using dynamic programming"""
    similarities = []
    min_len = min(len(human_seq), len(other_seq))
    
    # Create a memoization cache to avoid redundant calculations
    memo = {}
    
    for i in range(0, min_len - window_size + 1):
        human_window = human_seq[i:i+window_size]
        other_window = other_seq[i:i+window_size]
        
        # Check if we've already computed this window pair
        window_pair = (human_window, other_window)
        if window_pair in memo:
            sim = memo[window_pair]
        else:
            # Use DP for this window calculation
            sim, _ = dp_sequence_alignment(human_window, other_window)
            memo[window_pair] = sim
            
        similarities.append(sim)
    
    return similarities

def find_diagnostic_sites(human_seq, other_seq, species_name, max_sites=10):
    """Find species-specific diagnostic sites using DP alignment information"""
    min_len = min(len(human_seq), len(other_seq))
    
    # Get alignment information from DP
    _, dp_matrix = dp_sequence_alignment(human_seq[:min_len], other_seq[:min_len])
    
    print(f"\nDiagnostic Sites where {species_name} differs from human:")
    print("=" * 50)
    
    # Find positions with significant differences based on DP scores
    differences = []
    for pos in range(min_len):
        if human_seq[pos] != other_seq[pos]:
            # The difference at this position impacts alignment score
            differences.append(f"Position {pos+1}: Human={human_seq[pos]}, {species_name}={other_seq[pos]}")
            if len(differences) >= max_sites:
                differences.append(f"... and more differences")
                break
    
    for diff in differences:
        print(f"  {diff}")

def clean_sequence(seq):
    """Clean a DNA sequence by removing whitespace and normalizing to uppercase"""
    return ''.join(s.upper() for s in seq if s.upper() in 'ATGC')

def show():
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
    
    # Calculate overall similarity (now using DP)
    similarity, _ = dp_sequence_alignment(human_coi, species_seq)
    print(f"\nOverall DNA Similarity (DP-based): {similarity:.2f}%")
    
    # For comparison, also show the original method
    original_similarity = calculate_similarity(human_coi, species_seq)
    print(f"Original DNA Similarity: {original_similarity:.2f}%")
    
    # Find diagnostic sites
    find_diagnostic_sites(human_coi, species_seq, species_name)
    
    # Create sliding window similarity plot using DP approach
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