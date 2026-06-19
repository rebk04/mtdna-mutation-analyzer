from exceptions import (
    MutationFormatError,
    MutationPositionError,
    ORFSelectionError
)

from fasta_utils import (load_fasta, validate_dna)
from orf_utils import find_orfs
from mutation_utils import (
    translate_codon,
    parse_mutation_input,
    validate_mutation,
    find_affected_codon,
    apply_mutation,
    mutation_classification
)

from storage import export_mutations
from models import MutationRecord
from analysis import chi_square, p_value, analyze_orf

mutation_history = []

def analyze_and_record_mutation(orfs):
    selected_orf_name = input("Choose an ORF:").strip().upper()

    selected_orf = None

    for orf in orfs:
        if orf["name"] == selected_orf_name:
            selected_orf = orf
            break

    if selected_orf is None:
        raise ORFSelectionError("Invalid ORF selected!")

    mutation_str = input("Enter your mutation string: ")

    original_base, mutating_base, mutation_position = parse_mutation_input(mutation_str)
    validate_mutation(selected_orf, original_base, mutation_position)
    original_codon, codon_starting_index, mutation_base_position = find_affected_codon(selected_orf, mutation_position)
    mutated_codon = apply_mutation(original_codon, mutation_base_position, mutating_base)

    original_aa = translate_codon(original_codon)
    mutated_aa = translate_codon(mutated_codon)

    mutation_type = mutation_classification(original_aa, mutated_aa)

    mutation_history.append(
        MutationRecord(
            selected_orf["name"],
            mutation_str,
            original_codon,
            mutated_codon,
            original_aa,
            mutated_aa,
            mutation_type
        )
    )

    print("Original codon and its name:", original_codon, original_aa)
    print("Codon starting index: ", codon_starting_index)
    print("Position inside codon: ", mutation_base_position)
    print("Mutated codon and its name: ", mutated_codon, mutated_aa)
    print("Mutation type:", mutation_type)


dna = load_fasta("mtDNA.fasta")
validate_dna(dna)
orfs = find_orfs(dna)
fake_orf = {"sequence": dna, "start": 1, "stop": len(dna), "name": "FULL_DNA"}
expected, _ = analyze_orf(fake_orf)


while True:
    print("-" * 50)
    print("See DNA information options, press 0")
    print("Show detected ORFs, press 1")
    print("Add and analyse a mutation, press 2")
    print("Show mutation history, press 3")
    print("Export mutation history to CSV, press 4")
    print("Choose an ORF for batch analysis, press 5")
    print("Exit, press 6")

    choice = input("Choose option: ")

    if choice == "0":
        print("Show DNA length, press 0")
        print("Show first 100 bases, press 1")
        print("Show selected sequence range, press 2")
        print("Show GC content, press 3")
        print("Back to main menu, press 4")

        sub_choice = input("Input a number: ")
        if sub_choice == "0":
            print("DNA length:", len(dna), "bases")
        elif sub_choice == "1":
            print("First 100 bases:", dna[:100])
        elif sub_choice == "2":
            try:
                starting = int(input("What is your starting index of base? "))
                ending = int(input("What is your ending index of base? "))

                if starting < 1 or ending > len(dna) or starting > ending:
                    print("Invalid range!")
                else:
                    print("Sequence for selected range:", dna[starting - 1:ending])  # 0-based indexing conversion.
            except ValueError:
                print("Invalid input! Start and end positions must be numbers.")
        elif sub_choice == "3":
            base_counts = {base: dna.count(base) for base in {"A", "T", "C", "G", "N"}}
            valid_base_count = len(dna) - base_counts["N"]

            if valid_base_count == 0:
                print("GC content cannot be calculated because there are no valid A/T/C/G bases.")
            else:
                gc_content = ((base_counts["G"] + base_counts["C"]) / valid_base_count) * 100
                print(f"GC content: {gc_content:.2f}%")

    elif choice == "1":
        sorted_orfs = sorted(orfs, key=lambda orf: orf["start"])

        for orf in sorted_orfs:
            print(orf["name"], orf["start"], orf["stop"], "Length:", len(orf["sequence"]))

    elif choice == "2":
        try:
            analyze_and_record_mutation(orfs)
        except MutationFormatError as e:
            print(e)
        except MutationPositionError as e:
            print(e)
        except ORFSelectionError as e:
            print(e)

    elif choice == "3":
        if not mutation_history:
            print("No mutations recorded yet")
            add_choice = input("Would you like to add a mutation? yes/no: ").strip().upper()
            if add_choice == "YES":
                try:
                    analyze_and_record_mutation(orfs)
                except MutationFormatError as e:
                    print(e)
                except MutationPositionError as e:
                    print(e)
                except ORFSelectionError as e:
                    print(e)
        else:
            for record in mutation_history:
                print("ORF:", record.orf)
                print("Mutation:", record.mutation)
                print("Original codon:", record.original_codon)
                print("Mutated codon:", record.mutated_codon)
                print("Original AA:", record.original_amino_acid)
                print("Mutated AA:", record.mutated_amino_acid)
                print("Type:", record.mutation_type)
                print("-" * 30)

    elif choice == "4":
        if not mutation_history:
            print("No mutation history to export.")
            add_choice = input("Would you like to add a mutation? yes/no: ").strip().upper()
            if add_choice == "YES":
                try:
                    analyze_and_record_mutation(orfs)
                    print("Mutation added successfully!")
                    export_choice = input("Would you like to export your mutation? yes/no: ").strip().upper()
                    if export_choice == "YES":
                        export_mutations(mutation_history)
                except MutationFormatError as e:
                    print(e)
                except MutationPositionError as e:
                    print(e)
                except ORFSelectionError as e:
                    print(e)
        else:
            export_mutations(mutation_history)

    elif choice == "5":
        sorted_orfs = sorted(orfs, key=lambda orf: orf["start"])

        for orf in sorted_orfs:
            print(orf["name"], orf["start"], orf["stop"], "Length:", len(orf["sequence"]))

        selected_orf_name = input("Select an ORF to begin batch analysis: ").strip().upper()
        selected_orf = None

        for orf in orfs:
            if orf["name"] == selected_orf_name:
                selected_orf = orf
                break

        if selected_orf is None:
            raise ORFSelectionError("Invalid ORF selected!")

        mutation_counts, missense_positions = analyze_orf(selected_orf)
        for key, value in mutation_counts.items():
            print(key, value)

        #Normalization cause size of sequence is too big compared to selected orf:
        total_observed = sum(mutation_counts.values())
        total_expected = sum(expected.values())

        normalized_expected = {
            key: (expected[key] / total_expected) * total_observed
            for key in expected
        }

        chi_result = chi_square(mutation_counts, normalized_expected)
        print("Chi-square result:", chi_result)
        p = p_value(chi_result)
        print("P-value:", p)
        if p < 0.05:
            print("Statistically significant difference from full DNA distribution!")
        else:
            print("No statistically significant difference!")

        sub_choice = input("Would you like to also see the hotspot position? YES/NO").strip().upper()

        if sub_choice == "YES":
            sorted_positions = sorted(missense_positions.items(), key=lambda x: x[1], reverse=True)
            print(f"The hotspot position is: {sorted_positions[0][0]}, Missense count: {sorted_positions[0][1]}")

    elif choice == "6":
        print("Exiting program.")
        break
