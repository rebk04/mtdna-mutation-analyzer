import csv

def export_mutations(mutation_history):
    with open("mutations.csv", "w", newline="") as file:
        fieldnames = [
            "orf",
            "mutation",
            "original codon",
            "mutated codon",
            "original amino acid",
            "mutated amino acid",
            "mutation type"
        ]

        writer = csv.DictWriter(file, fieldnames=fieldnames)
        writer.writeheader()

        writer.writerows(record.to_dict() for record in mutation_history)
        print("Mutation history exported to mutations.csv")
