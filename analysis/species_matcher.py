from Bio import SeqIO

def load_human_cytochrome_b(file_path):
    """Load human Cytochrome B sequence from a text file (FASTA format)."""
    try:
        with open(file_path, "r") as f:
            for record in SeqIO.parse(f, "fasta"):
                if "Homo sapiens" in record.description or "human" in record.description.lower():
                    return str(record.seq)
    except FileNotFoundError:
        raise FileNotFoundError(f"File '{file_path}' not found. Ensure it's in the same directory as this script.")
    raise ValueError("Human Cytochrome B sequence not found in the file.")

def needleman_wunsch(seq1, seq2, match_score=1, mismatch_penalty=-1, gap_penalty=-1):
    """Needleman-Wunsch global alignment algorithm using dynamic programming."""
    m, n = len(seq1), len(seq2)
    
    # Initialize scoring matrix
    score = [[0] * (n + 1) for _ in range(m + 1)]

    # Initialize first row and column
    for i in range(m + 1):
        score[i][0] = i * gap_penalty
    for j in range(n + 1):
        score[0][j] = j * gap_penalty

    # Fill the scoring matrix
    for i in range(1, m + 1):
        for j in range(1, n + 1):
            match = score[i - 1][j - 1] + (match_score if seq1[i - 1] == seq2[j - 1] else mismatch_penalty)
            delete = score[i - 1][j] + gap_penalty
            insert = score[i][j - 1] + gap_penalty
            score[i][j] = max(match, delete, insert)

    # Traceback to compute alignment
    aligned_seq1 = ""
    aligned_seq2 = ""
    i, j = m, n

    while i > 0 and j > 0:
        current = score[i][j]
        diagonal = score[i - 1][j - 1]
        up = score[i - 1][j]
        left = score[i][j - 1]

        if current == diagonal + (match_score if seq1[i - 1] == seq2[j - 1] else mismatch_penalty):
            aligned_seq1 = seq1[i - 1] + aligned_seq1
            aligned_seq2 = seq2[j - 1] + aligned_seq2
            i -= 1
            j -= 1
        elif current == up + gap_penalty:
            aligned_seq1 = seq1[i - 1] + aligned_seq1
            aligned_seq2 = "-" + aligned_seq2
            i -= 1
        else:
            aligned_seq1 = "-" + aligned_seq1
            aligned_seq2 = seq2[j - 1] + aligned_seq2
            j -= 1

    # Fill the rest
    while i > 0:
        aligned_seq1 = seq1[i - 1] + aligned_seq1
        aligned_seq2 = "-" + aligned_seq2
        i -= 1
    while j > 0:
        aligned_seq1 = "-" + aligned_seq1
        aligned_seq2 = seq2[j - 1] + aligned_seq2
        j -= 1

    return aligned_seq1, aligned_seq2

def calculate_similarity(seq1, seq2):
    """Calculate similarity percentage using Needleman-Wunsch alignment."""
    aligned1, aligned2 = needleman_wunsch(seq1, seq2)

    # Count matches
    matches = sum(a == b for a, b in zip(aligned1, aligned2) if a != '-' and b != '-')
    similarity = (matches / max(len(seq1), len(seq2))) * 100
    return similarity

def main():
    # File path
    human_cytb_file = "data/cytochrome_human.txt"  # Change as per your setup
    input_dna = input("Enter the DNA sequence to compare: ").strip().upper()
    
    # Validate input
    if not all(base in "ATCG" for base in input_dna):
        print("Error: Input sequence must only contain A, T, C, or G.")
        return

    try:
        human_cytb = load_human_cytochrome_b(human_cytb_file)

        similarity = calculate_similarity(input_dna, human_cytb)
        print(f"\nInput DNA is {similarity:.2f}% similar to human Cytochrome B (using DP-based alignment).")

    except Exception as e:
        print(f"Error: {e}")

if __name__ == "__main__":
    main()
