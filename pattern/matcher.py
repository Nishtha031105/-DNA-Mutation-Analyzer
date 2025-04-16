def search_in_hbb_transcript(input_seq: str):
    try:
        # Load reference HBB sequence
        with open("data/mutation_sequences/Chromosome11.txt", "r") as f:
            ref_seq = f.read().replace("\n", "").strip().upper()
        # Find first ATG in the input sequence
        start_idx = input_seq.find("ATG")
        if start_idx == -1:
            raise ValueError("Start codon 'ATG' not found in input sequence.")
        pattern = input_seq[start_idx:].upper()
        start_idx2=ref_seq.find("ATG")
        reference_seq=ref_seq[start_idx2:]
        # KMP setup
        def compute_lps(pattern):
            lps = [0] * len(pattern)
            length = 0
            i = 1
            while i < len(pattern):
                if pattern[i] == pattern[length]:
                    length += 1
                    lps[i] = length
                    i += 1
                else:
                    if length != 0:
                        length = lps[length - 1]
                    else:
                        lps[i] = 0
                        i += 1
            return lps

        def kmp_search(text, pattern):
            lps = compute_lps(pattern)
            i = j = 0
            while i < len(text):
                if pattern[j] == text[i]:
                    i += 1
                    j += 1
                if j == len(pattern):
                    return True  # Exact match found
                elif i < len(text) and pattern[j] != text[i]:
                    if j != 0:
                        j = lps[j - 1]
                    else:
                        i += 1
            return False  # No match

        if kmp_search(reference_seq, pattern):
            return -1  # Exact match found

        # If not matched, find all mismatch positions
        mismatch_indices = []
        for i in range(min(len(pattern), len(reference_seq))):
            if pattern[i] != reference_seq[i]:
                mismatch_indices.append(i)

        return mismatch_indices

    except Exception as e:
        print(f"Error: {e}")
        return None
    
result = search_in_hbb_transcript("ACATTTGCTTCTGACACAACTGTGTTCACTAGCAACCTCAAACAGACACCATGGTGCATCTGACTCCTGTGGAGAAGACTGCCGTTACTGCCCTGTGG")
print(result)
