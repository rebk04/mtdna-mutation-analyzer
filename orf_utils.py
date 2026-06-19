from decorators import measure_time

@measure_time
def find_orfs(sequence):
    orfs = []
    orf_count = 1
    i = 0

    while i < len(sequence) - 2:
        codon = sequence[i:i + 3]

        if codon == "ATG":
            start_pos = i
            found_stop = False
            for j in range(start_pos + 3, len(sequence) -2, 3):
                stop_codon = sequence[j:j + 3]
                if stop_codon == "TAA" or stop_codon == "TAG" or stop_codon == "TGA":
                    orf_seq = sequence[start_pos: j + 3]
                    if len(orf_seq) >= 90:
                        orfs.append({
                            "name": f"ORF_{orf_count}",
                            "start": start_pos + 1,
                            "stop": j + 3,
                            "sequence": sequence[start_pos:j + 3]
                        })
                        orf_count += 1
                    i = j + 3
                    found_stop = True
                    break
            if not found_stop:
                i += 1
        else:
            i +=1
    return orfs
