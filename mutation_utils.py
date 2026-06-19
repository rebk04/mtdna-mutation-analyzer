import re

from exceptions import (MutationFormatError,
    MutationPositionError)
from constants import CODON_TABLE

def translate_codon(codon):
    return CODON_TABLE.get(codon, "Unknown")

def parse_mutation_input(mutation_str):
    mutation_str = mutation_str.strip().upper()

    pattern = r"^([ACTG])(\d+)([ACTG])$"
    match = re.match(pattern, mutation_str)

    if not match:
        raise MutationFormatError("Invalid mutation format! Use format like A25G")

    original_base = match.group(1)
    mutation_position = int(match.group(2))
    mutating_base = match.group(3)

    return original_base, mutating_base, mutation_position

def validate_mutation(selected_orf, original_base, mutation_position):
    orf_seq = selected_orf["sequence"]
    if mutation_position < 1 or mutation_position > len(orf_seq):
        raise MutationPositionError("Mutation position is outside of selected ORF!")

    if orf_seq[mutation_position -1] != original_base:
        raise MutationPositionError("Original base does not match ORF sequence!")

def find_affected_codon(selected_orf, mutation_position): #use paper and pen for 9th base if you forget.
    mutation_index = mutation_position - 1 #0 base convertion (for indexing). #if mutation_pos = 9 -> 8
    codon_number_index = mutation_index // 3 #8//3 = 2nd indexed codon (3rd codon correct).
    codon_starting_index = (codon_number_index*3) # 2nd indexed codon*3 -> affected codon first index.
    mutation_base_position = mutation_index%3 #8%3 = 2 -> codon's 3rd base which is correct for 9th base.

    original_codon = selected_orf["sequence"][codon_starting_index: codon_starting_index + 3]

    return original_codon, codon_starting_index, mutation_base_position

def apply_mutation(original_codon, mutation_base_position, mutating_base):
    if mutation_base_position == 0:
        mutated_codon = mutating_base + original_codon[1:3]
    elif mutation_base_position == 1:
        mutated_codon = original_codon[0] + mutating_base + original_codon[2]
    else:
        mutated_codon = original_codon[0:2] + mutating_base

    return mutated_codon

def mutation_classification(original_aa, mutated_aa):
    if mutated_aa == "STOP":
        return "Nonsense"
    elif original_aa == mutated_aa:
        return "Silent"
    else:
        return "Missense"
