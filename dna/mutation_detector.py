import sys
import os
sys.path.append(os.path.abspath(os.path.join(os.path.dirname(__file__), '..')))

def detect_hbb_mutations(mismatch_indices):
    """
    Takes a list of mismatch indices and prints known mutation messages
    for the HBB gene based on known mutation positions.
    """
    # Define mutation rules: absolute index → {expected/actual codon: message}
    messages = []
    MUTATION_RULES = {
        18: {
            "AAG": "In Codon 6: Normal (GAG → Glutamic acid) has been mutated to AAG that is Lysine. Haemoglobin C mutation detected!"
        },
        19: {
            "GTG": "In Codon 6: Normal (GAG → Glutamic acid) has been mutated to GTG that is Valine. Sickle cell mutation detected!"
        },
        78: {
            "AAG": "In Codon 26: Normal (GAG → Glutamic acid) has been mutated to AAG which is Lysine. Hemoglobin E mutation detected!"
        },
        102: {
        "TGG": "Codon 34: Normal (AGG → Arginine) has been mutated to TGG which is Tryptophan. Possible β-thalassemia mutation (rare)."
        },
        117: {
            "TAG": "In Codon 39: Normal (CAG → Glutamine) has been mutated to TAG which is STOP. β-Thalassemia nonsense mutation detected!"
        },
        123: {
        "TAA": "Codon 41: Normal (CAA → Glutamine) has been mutated to TAA (STOP). Known β⁰-thalassemia mutation."
        },
        165: {
        "TAG": "Codon 55: Normal (CAG → Glutamine) has been mutated to TAG (STOP codon). β⁰-thalassemia (early truncation)."
        },
        195: {
        "TGG": "Codon 65: Normal (CGG → Arginine) has been mutated to TGG which is Tryptophan. Variant hemoglobin, likely unstable."
        },
        204: {
        "TAA": "Codon 68: Normal (CAA → Glutamine) has been mutated to TAA (STOP). Severe β⁰-thalassemia mutation."
        },
        234: {
        "AAC": "Codon 78: Normal (GAC → Aspartic acid) mutated to AAC which is Asparagine. Possible mild hemoglobin variant."
        },
        273: {
        "GAG": "Codon 91: Normal (AAG → Lysine) has been mutated to GAG which is Glutamic acid. Rare hemoglobin variant (Hb Kansas)."
        }
    }

    # Load reference sequence (same as in matcher)
    with open("C:\\Users\\Nishtha\\Desktop\\DNA_MUTATION_ANALYZER\\-DNA-Mutation-Analyzer\\data\\mutation_sequences\\Chromosome11.txt", "r") as f:
        full_ref_seq = f.read().replace("\n", "").strip().upper()

    start = full_ref_seq.find("ATG")
    if start == -1:
        return "Start codon not found in reference sequence."

    # Trim reference to start from ATG
    ref_seq = full_ref_seq[start:]

    for idx in mismatch_indices:
        if idx in MUTATION_RULES:
            # If the index is a known mutation, we proceed
            for expected_codon, msg in MUTATION_RULES[idx].items():
                messages.append(msg)
        else:
            messages.append(f"Index {idx}: Unknown mutation site.")
    return messages

#matcher.py ka code

def search_in_hbb_transcript(input_seq: str):
    try:
        # Load reference HBB sequence
        with open("C:\\Users\\Nishtha\\Desktop\\DNA_MUTATION_ANALYZER\\-DNA-Mutation-Analyzer\\data\\mutation_sequences\\Chromosome11.txt", "r") as f:
            ref_seq = f.read().replace("\n", "").strip().upper()
        # Find first ATG in the input sequence
        start_idx = input_seq.find("ATG")
        if start_idx == -1:
            return "Start codon 'ATG' not found in input sequence."
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
            return ["✅ Congratulations! The input sequence is normal."]

        # If not matched, find all mismatch positions
        mismatch_indices = []
        print("⚠️ Uh Oh! Mutations detected.")
        for i in range(min(len(pattern), len(reference_seq))):
            if pattern[i] != reference_seq[i]:
                mismatch_indices.append(i)
        return detect_hbb_mutations(mismatch_indices)


    except Exception as e:
        print(f"Error: {e}")
        return None
        
example_dna="ACATTTGCTTCTGACACAACTGTGTTCACTAGCAACCTCAAACAGACACCATGGTGCATCTGACTCCTGTGGAGAAGTCTGCCGTTACTGCCCTGTGG"  
def get_result(input_dna):
    return search_in_hbb_transcript(input_dna)