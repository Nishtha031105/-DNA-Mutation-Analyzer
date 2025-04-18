import sys
import os
sys.path.append(os.path.abspath(os.path.join(os.path.dirname(__file__), '..')))

def detect_tp53_mutations(mismatch_indices):
    """
    Takes a list of mismatch indices and prints known mutation messages
    for the TP53 gene based on known mutation positions.
    """
    # Define mutation rules: absolute index → {mutated codon: message}
    messages = []
    MUTATION_RULES = {
        36: {
            "CAT": "In Codon 13: Normal (CCT → Proline) has been mutated to CAT which is Histidine. Associated with early-onset cancers."
        },
        
        114: {
            "TAT": "In Codon 39: Normal (CAT → Histidine) has been mutated to TAT which is Tyrosine. Found in Li-Fraumeni syndrome families."
        },
        156: {
            "TAG": "In Codon 53: Normal (TGG → Tryptophan) has been mutated to TAG which is STOP. Truncated p53 protein associated with aggressive tumors."
        },
        
        272: {
            "CAT": "In Codon 91: Normal (CGT → Arginine) has been mutated to CAT which is Histidine. Common mutation in multiple cancer types (R273H)."
        },
      
        336: {
            "TGT": "In Codon 113: Normal (CGT → Arginine) has been mutated to TGT which is Cysteine. Disrupts p53 protein conformation."
        },
        
        457: {
            "TTT": "In Codon 153: Normal (TTC → Phenylalanine) has been mutated to TTT which is still Phenylalanine. Silent mutation, no amino acid change."
        }
    }

    # Load reference sequence (same as in matcher)
    with open("C:\\Users\\Nishtha\\Desktop\\DNA_MUTATION_ANALYZER\\-DNA-Mutation-Analyzer\\data\\mutation_sequences\\chromosome17.txt", "r") as f:
        full_ref_seq = f.read().replace("\n", "").strip().upper()

    start = full_ref_seq.find("ATG")
    if start == -1:
        return "Start codon not found in reference sequence."

    # Trim reference to start from ATG
    ref_seq = full_ref_seq[start:]

    # Iterate through each mismatch and process at the codon level
    for idx in mismatch_indices:
        codon_start = idx - (idx % 3)  # Align to the start of the codon
        codon = ref_seq[codon_start:codon_start + 3]  # Get the reference codon
        if idx in MUTATION_RULES:
            # If the index is a known mutation, we proceed
            for mutated_codon, msg in MUTATION_RULES[idx].items():
                messages.append(msg)
        else:
            messages.append(f"Index {idx}: Unknown mutation site, codon {codon}.")
    return messages

def search_in_tp53_transcript(input_seq: str):
    try:
        # Load reference TP53 sequence
        with open("C:\\Users\\Nishtha\\Desktop\\DNA_MUTATION_ANALYZER\\-DNA-Mutation-Analyzer\\data\\mutation_sequences\\chromosome17.txt", "r") as f:
            ref_seq = f.read().replace("\n", "").strip().upper()

        # Find first ATG in the input sequence
        start_idx = input_seq.find("ATG")
        if start_idx == -1:
            return "Start codon 'ATG' not found in input sequence."
        pattern = input_seq[start_idx:].upper()

        start_idx2 = ref_seq.find("ATG")
        reference_seq = ref_seq[start_idx2:]

        # KMP setup (for exact match)
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
            return ["✅ Congratulations! The input sequence is normal."]

        # If not matched, find all mismatch positions
        mismatch_indices = []
        print("⚠ Uh Oh! Mutations detected.")
        for i in range(0, min(len(pattern), len(reference_seq)), 3):  # Check in chunks of 3
            if pattern[i:i+3] != reference_seq[i:i+3]:  # Compare codon by codon
                mismatch_indices.append(i)

        return detect_tp53_mutations(mismatch_indices)

    except Exception as e:
        print(f"Error: {e}")
        return None
        
example_dna = "ATGGAGGAGCCGCAGTCAGATCCTAGCGTCGAGCCCCATCTGAGTCAGGAAACATTTTCAGACCTATGGAAGAAATCGGTAAGAGGTGCGTGTTTGTGCCT"  
def get_result2(input_dna):
    return search_in_tp53_transcript(input_dna)
