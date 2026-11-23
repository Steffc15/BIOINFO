"""Lab 7 - Restriction enzyme digestion.

Input : - the recognition site (and cleavage position) of 5 restriction enzymes
        - a bacterial DNA sequence from NCBI (E. coli K-12 MG1655, NC_000913.3)
Output: - number of cleavages for each enzyme
        - cleavage positions and fragment lengths
        - simulated agarose gel electrophoresis (one lane per enzyme)

Usage: python ex1.py [accession] [start] [stop]
"""
import math
import os
import sys
import urllib.request

import matplotlib.pyplot as plt

ACCESSION = sys.argv[1] if len(sys.argv) > 1 else "NC_000913.3"
SEQ_START = int(sys.argv[2]) if len(sys.argv) > 2 else 1
SEQ_STOP = int(sys.argv[3]) if len(sys.argv) > 3 else 10000

HERE = os.path.dirname(os.path.abspath(__file__))
FASTA = os.path.join(HERE, f"{ACCESSION}_{SEQ_START}-{SEQ_STOP}.fasta")
URL = ("https://eutils.ncbi.nlm.nih.gov/entrez/eutils/efetch.fcgi"
       f"?db=nuccore&id={ACCESSION}&rettype=fasta&retmode=text"
       f"&seq_start={SEQ_START}&seq_stop={SEQ_STOP}")

# name: (recognition site 5'->3', cut position on the top strand,
#        cut position on the bottom strand), both counted from the 5' end
ENZYMES = {
    "EcoRI":   ("GAATTC", 1, 1),   # G^AATTC
    "BamHI":   ("GGATCC", 1, 1),   # G^GATCC
    "HindIII": ("AAGCTT", 1, 1),   # A^AGCTT
    "TaqI":    ("TCGA",   1, 1),   # T^CGA
    "HaeIII":  ("GGCC",   2, 2),   # GG^CC (blunt)
}

COMPLEMENT = str.maketrans("ACGT", "TGCA")
LADDER = [100, 200, 300, 500, 750, 1000, 1500, 2000, 3000, 5000, 10000]
GEL_MIN, GEL_MAX = 50, 12000


def reverse_complement(seq):
    return seq.translate(COMPLEMENT)[::-1]


def read_fasta():
    if not os.path.exists(FASTA):
        with urllib.request.urlopen(URL) as resp:
            data = resp.read().decode()
        with open(FASTA, "w") as f:
            f.write(data)
    with open(FASTA) as f:
        lines = f.read().splitlines()
    return lines[0][1:], "".join(l.strip() for l in lines[1:]).upper()


def find_sites(seq, site):
    positions, i = [], seq.find(site)
    while i != -1:
        positions.append(i)
        i = seq.find(site, i + 1)
    return positions


def cleavage_positions(seq, site, cut_top, cut_bottom):
    """Returns the top-strand cut positions: a cut at p means between seq[p-1] and seq[p]."""
    cuts = {i + cut_top for i in find_sites(seq, site)}
    rc = reverse_complement(site)
    if rc != site:  # non-palindromic site: also search it on the bottom strand
        cuts |= {i + len(site) - cut_bottom for i in find_sites(seq, rc)}
    return sorted(c for c in cuts if 0 < c < len(seq))


def fragments(seq_len, cuts):
    bounds = [0] + cuts + [seq_len]
    return [(bounds[i] + 1, bounds[i + 1], bounds[i + 1] - bounds[i])
            for i in range(len(bounds) - 1)]


def ends_type(site, cut_top):
    overhang = len(site) - 2 * cut_top
    if overhang == 0:
        return "blunt ends"
    return f"5' sticky ends ({overhang} nt overhang)" if overhang > 0 else \
        f"3' sticky ends ({-overhang} nt overhang)"


def migration(length, gel_length=100.0):
    lo, hi = math.log10(GEL_MIN), math.log10(GEL_MAX)
    length = min(max(length, GEL_MIN), GEL_MAX)
    return 4 + (hi - math.log10(length)) / (hi - lo) * (gel_length - 8)


def draw_gel(lanes, title):
    lanes = [("Ladder", LADDER)] + lanes
    fig, ax = plt.subplots(figsize=(12, 8))
    ax.set_facecolor("#111122")
    max_len = max(max(sizes) for _, sizes in lanes)

    for x, (name, sizes) in enumerate(lanes):
        ax.add_patch(plt.Rectangle((x - 0.35, -3), 0.7, 3, color="#444466"))
        for size in sizes:
            y = migration(size)
            # brightness ~ amount of DNA in the band (longer fragments bind more dye)
            alpha = 0.35 + 0.65 * math.log10(size) / math.log10(max_len)
            color = "#ffaa33" if name == "Ladder" else "#33ffcc"
            ax.hlines(y, x - 0.32, x + 0.32, color=color, linewidth=4, alpha=min(alpha, 1))
            if name == "Ladder":
                ax.text(-0.55, y, f"{size}", ha="right", va="center", color="white", fontsize=8)

    ax.set_xlim(-1.3, len(lanes) - 0.5)
    ax.set_ylim(102, -5)
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
    print(f"Length  : {len(seq)} bp\n")

    lanes, all_cuts = [], set()
    for name, (site, cut_top, cut_bottom) in ENZYMES.items():
        cuts = cleavage_positions(seq, site, cut_top, cut_bottom)
        frags = fragments(len(seq), cuts)
        all_cuts |= set(cuts)
        lanes.append((name, [length for _, _, length in frags]))

        bottom = reverse_complement(site)[::-1]
        print("=" * 64)
        print(f"{name}   5'-{site[:cut_top]}^{site[cut_top:]}-3'   "
              f"3'-{bottom[:len(site) - cut_bottom]}^{bottom[len(site) - cut_bottom:]}-5'   "
              f"({ends_type(site, cut_top)})")
        print(f"Number of cleavages: {len(cuts)}")
        print(f"Cleavage positions : {cuts if cuts else '-'}")
        print(f"Fragments ({len(frags)}):")
        print("   #  |  Start |   End  | Length (bp)")
        for k, (start, end, length) in enumerate(frags, start=1):
            print(f"  {k:3d} | {start:6d} | {end:6d} | {length:6d}")
        print()

    combined = fragments(len(seq), sorted(all_cuts))
    lanes.append(("All 5", [length for _, _, length in combined]))
    print("=" * 64)
    print(f"Digestion with all 5 enzymes: {len(all_cuts)} cleavages, "
          f"{len(combined)} fragments")

    draw_gel(lanes, f"Simulated restriction digest - {header[:60]}")
