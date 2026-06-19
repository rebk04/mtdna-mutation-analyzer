class MutationRecord:
    def __init__(self, orf, mutation, original_codon, mutated_codon, original_amino_acid, mutated_amino_acid, mutation_type):
        self.orf = orf
        self.mutation = mutation
        self.original_codon = original_codon
        self.mutated_codon = mutated_codon
        self.original_amino_acid = original_amino_acid
        self.mutated_amino_acid = mutated_amino_acid
        self.mutation_type = mutation_type

    def to_dict(self):
        return {
            "orf": self.orf,
            "mutation": self.mutation,
            "original codon": self.original_codon,
            "mutated codon": self.mutated_codon,
            "original amino acid": self.original_amino_acid,
            "mutated amino acid": self.mutated_amino_acid,
            "mutation type": self.mutation_type
        }
