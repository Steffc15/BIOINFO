# Ex 2: detect codon frequencies in all reading frames, sort them,
# and plot the 10 most frequent codons.
from collections import Counter

import matplotlib.pyplot as plt

from ex0 import translate_codon
from ex1 import OUT_FILE, read_fasta


def reverse_complement(seq):
    return seq[::-1].translate(str.maketrans("ACGT", "TGCA"))


def codon_frequencies(seq):
    """Count codons in all 6 reading frames (3 forward + 3 reverse complement)."""
    counts = Counter()
    for strand in (seq, reverse_complement(seq)):
        for frame in range(3):
            for i in range(frame, len(strand) - 2, 3):
                codon = strand[i:i + 3]
                if set(codon) <= set("ACGT"):
                    counts[codon] += 1
    return counts


def top_codons(n=10):
    header, seq = read_fasta(OUT_FILE)
    counts = codon_frequencies(seq)
    return header, counts, counts.most_common(n)


if __name__ == "__main__":
    header, counts, top10 = top_codons()
    total = sum(counts.values())

    print(f"Genome: {header}")
    print(f"Total codons (6 frames): {total}\n")
    print("All codons, sorted from highest to lowest:")
    for codon, c in counts.most_common():
        print(f"  {codon} ({translate_codon(codon)}): {c}  ({c / total:.2%})")

    labels = [f"{c}\n({translate_codon(c)})" for c, _ in top10]
    values = [v for _, v in top10]
    plt.figure(figsize=(10, 5))
    plt.bar(labels, values, color="steelblue")
    plt.title("Top 10 most frequent codons - Hepatitis delta virus (NC_001653.2)")
    plt.xlabel("Codon (amino acid)")
    plt.ylabel("Frequency")
    plt.tight_layout()
    plt.savefig("top10_codons.png")
    plt.show()
