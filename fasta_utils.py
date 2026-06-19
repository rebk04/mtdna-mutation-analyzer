from exceptions import DNAValidationError

def load_fasta(filename):
    with open(filename, "r") as file:
        lines = file.readlines()
    sequence_lines = []

    for line in lines:
        line = line.strip() #Cleaning.

        if not line:
            continue

        if line.startswith(">"):
            continue #Pass the line if the line who starts with '>'.

        sequence_lines.append(line)

    sequence = "".join(sequence_lines)
    sequence = sequence.upper()
    return sequence

def validate_dna(sequence):
    valid_bases = {"A","T","C","G","N"}
    for index, base in enumerate(sequence):
        if base not in valid_bases:
            raise DNAValidationError(f"Invalid base '{base}' at position {index + 1} ")
    return True