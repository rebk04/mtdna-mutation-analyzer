from mutation_utils import find_affected_codon, apply_mutation, mutation_classification, translate_codon
import math

def chi_square(observed, expected):
    chi_result = 0
    for key in observed:
        chi_result += ((observed[key] - expected[key]) **2) / expected[key]
    return chi_result

def p_value(chi_result):
    p_value = math.exp(-chi_result / 2)
    return p_value

def analyze_orf(orf):
    mutation_counts = {
        "Silent": 0,
        "Nonsense": 0,
        "Missense": 0
    }

    missense_positions = {}

    for i in range(len(orf["sequence"])):
        for j in ["A", "T", "C", "G"]:
            if j == orf["sequence"][i]:
                continue
            original_codon, codon_start, base_pos = find_affected_codon(selected_orf=orf, mutation_position=i + 1)
            mutated_codon = apply_mutation(original_codon=original_codon, mutation_base_position=base_pos,
                                           mutating_base=j)
            original_aa = translate_codon(original_codon)
            mutated_aa = translate_codon(mutated_codon)
            mutation_type = mutation_classification(original_aa, mutated_aa)
            mutation_counts[mutation_type] += 1

            if mutation_type == "Missense":
                if i+1 not in missense_positions:
                    missense_positions[i+1] = 0
                missense_positions[i+1] += 1

    return mutation_counts, missense_positions


