import math
import os
import random
import urllib.request

import matplotlib.pyplot as plt

ACCESSION = "NC_001653.2"
HERE = os.path.dirname(os.path.abspath(__file__))
FASTA = os.path.join(HERE, "hdv_genome.fasta")
URL = ("https://eutils.ncbi.nlm.nih.gov/entrez/eutils/efetch.fcgi"
       f"?db=nuccore&id={ACCESSION}&rettype=fasta&retmode=text")

N_SAMPLES = 10
SAMPLE_MIN, SAMPLE_MAX = 100, 3000
AVG_NT_WEIGHT = 650
LADDER = [100, 200, 300, 400, 500, 600, 700, 800, 900, 1000, 1500, 2000, 3000]


def read_fasta():
    if not os.path.exists(FASTA):
        with urllib.request.urlopen(URL) as resp:
            with open(FASTA, "w") as f:
                f.write(resp.read().decode())
    with open(FASTA) as f:
        lines = f.read().splitlines()
    return lines[0][1:], "".join(l.strip() for l in lines[1:]).upper()


def take_samples(seq, n, rng):
    max_len = min(SAMPLE_MAX, len(seq))
    samples = []
    for _ in range(n):
        length = rng.randint(SAMPLE_MIN, max_len)
        start = rng.randint(0, len(seq) - length)
        samples.append(seq[start:start + length])
    return samples


def migration_distance(length, gel_length=100.0):
    lo, hi = math.log10(min(LADDER)), math.log10(max(LADDER))
    rel = (hi - math.log10(length)) / (hi - lo)
    return 5 + rel * (gel_length - 10)


def draw_gel(samples, title):
    lanes = [("Ladder", LADDER)] + [(f"S{i + 1}", [len(s)]) for i, s in enumerate(samples)]

    fig, ax = plt.subplots(figsize=(12, 7))
    ax.set_facecolor("#111122")
    for x, (name, sizes) in enumerate(lanes):
        ax.add_patch(plt.Rectangle((x - 0.35, -2), 0.7, 2, color="#444466"))
        for size in sizes:
            y = migration_distance(size)
            alpha = 0.45 + 0.55 * size / SAMPLE_MAX
            color = "#ffaa33" if name == "Ladder" else "#33ffcc"
            ax.hlines(y, x - 0.32, x + 0.32, color=color, linewidth=5, alpha=alpha)
            if name == "Ladder":
                ax.text(-0.55, y, f"{size}", ha="right", va="center", color="white", fontsize=8)
            else:
                ax.text(x, y + 2.5, f"{size} bp", ha="center", color="white", fontsize=8)

    ax.set_xlim(-1.2, len(lanes) - 0.5)
    ax.set_ylim(102, -4)
    ax.set_xticks(range(len(lanes)))
    ax.set_xticklabels([name for name, _ in lanes])
    ax.set_yticks([])
    ax.set_ylabel("Migration direction  (-)  →  (+)")
    ax.set_title(title)
    plt.tight_layout()
    plt.show()


if __name__ == "__main__":
    header, seq = read_fasta()
    print(f"Sequence: {header}")
    print(f"Length  : {len(seq)} nt\n")

    samples = take_samples(seq, N_SAMPLES, random.Random())

    print("Sample | Length (bp) | Weight (Da) | Distance (% of gel)")
    for i, s in enumerate(samples):
        print(f"  S{i + 1:<4} | {len(s):11d} | {len(s) * AVG_NT_WEIGHT:11d} | "
              f"{migration_distance(len(s)):6.1f}")

    order = sorted(range(N_SAMPLES), key=lambda i: len(samples[i]))
    print("\nFastest to slowest: " + " > ".join(f"S{i + 1}" for i in order))

    draw_gel(samples, f"Simulated gel electrophoresis - {header}")
