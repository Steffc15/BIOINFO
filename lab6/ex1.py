"""Lab 6 - Transposable elements.

An artificial DNA sequence (200-400 bp) is built in which 3-6 transposons are
simulated. Each transposon is modelled as a DNA transposon:

    5'-[left TIR]--------- body ---------[right TIR]-3'
        right TIR = reverse complement of left TIR  (terminal inverted repeats)

Three cases are simulated, each in its own sequence:
    1. isolated    - transposons that do not touch each other
    2. overlapping - two transposons share a region (+ one isolated element)
    3. nested      - one transposon inserted inside another (+ one isolated)

The detector does not know where the elements are. It scans every window of
plausible transposon length and reports the windows whose ends form a pair of
terminal inverted repeats. Afterwards it classifies the detected elements as
isolated / overlapping / nested and compares them with the ground truth.
"""
import random
import sys

import matplotlib.pyplot as plt
from matplotlib.patches import Rectangle

SEED = int(sys.argv[1]) if len(sys.argv) > 1 else 7
ALPHABET = "ACGT"
COMPLEMENT = {"A": "T", "T": "A", "C": "G", "G": "C"}

SEQ_MIN, SEQ_MAX = 200, 400
TIR_LEN = 10
TE_MIN, TE_MAX = 30, 45
DETECT_MIN, DETECT_MAX = 2 * TIR_LEN + 5, 100
OVERLAP_LEN = 15


def reverse_complement(seq):
    return "".join(COMPLEMENT[b] for b in reversed(seq))


def random_dna(n, rng):
    return "".join(rng.choice(ALPHABET) for _ in range(n))


def random_tir(rng):
    while True:
        tir = random_dna(TIR_LEN, rng)
        if tir != reverse_complement(tir):
            return tir


def make_transposon(rng, left_tir=None):
    length = rng.randint(TE_MIN, TE_MAX)
    left = left_tir or random_tir(rng)
    return left + random_dna(length - 2 * TIR_LEN, rng) + reverse_complement(left)


# ---------------------------------------------------------------- building blocks
# Every block returns (sequence, [(name, start, end), ...]) with 0-based,
# inclusive coordinates relative to the block.

def block_single(rng, name):
    te = make_transposon(rng)
    return te, [(name, 0, len(te) - 1)]


def block_overlapping(rng, name_a, name_b):
    """TE A ends with a region of OVERLAP_LEN bp that is also the start of TE B."""
    tir_a = random_tir(rng)
    len_a = rng.randint(TE_MIN, TE_MAX)
    shared = random_dna(OVERLAP_LEN - TIR_LEN, rng) + reverse_complement(tir_a)
    te_a = tir_a + random_dna(len_a - TIR_LEN - OVERLAP_LEN, rng) + shared

    tir_b = shared[:TIR_LEN]
    len_b = rng.randint(TE_MIN, TE_MAX)
    te_b = shared + random_dna(len_b - OVERLAP_LEN - TIR_LEN, rng) + reverse_complement(tir_b)

    start_b = len(te_a) - OVERLAP_LEN
    seq = te_a[:start_b] + te_b
    return seq, [(name_a, 0, len(te_a) - 1), (name_b, start_b, len(seq) - 1)]


def block_nested(rng, name_outer, name_inner):
    """TE inner is inserted in the body of TE outer (between its TIRs)."""
    outer = make_transposon(rng)
    inner = make_transposon(rng)
    pos = rng.randint(TIR_LEN + 2, len(outer) - TIR_LEN - 2)
    seq = outer[:pos] + inner + outer[pos:]
    return seq, [(name_outer, 0, len(seq) - 1), (name_inner, pos, pos + len(inner) - 1)]


def assemble(blocks, total_len, rng):
    """Places the blocks in a random host sequence, separated by random spacers."""
    used = sum(len(b[0]) for b in blocks)
    free = total_len - used
    cuts = sorted(rng.sample(range(1, free), len(blocks)))
    spacers = [b - a for a, b in zip([0] + cuts, cuts + [free])]

    seq, truth = "", []
    for spacer, (block_seq, elements) in zip(spacers, blocks):
        seq += random_dna(spacer, rng)
        offset = len(seq)
        seq += block_seq
        truth += [(n, s + offset, e + offset) for n, s, e in elements]
    seq += random_dna(spacers[-1], rng)
    return seq, truth


def tir_is_maximal(seq, start, end):
    """True if the inverted repeat cannot be extended outwards by one base."""
    if start == 0 or end == len(seq) - 1:
        return True
    return COMPLEMENT[seq[start - 1]] != seq[end + 1]


def build_case(kind, rng):
    while True:
        total_len = rng.randint(SEQ_MIN, SEQ_MAX)
        if kind == "isolated":
            blocks = [block_single(rng, f"TE{i + 1}") for i in range(rng.randint(3, 4))]
        elif kind == "overlapping":
            blocks = [block_overlapping(rng, "TE1", "TE2"), block_single(rng, "TE3")]
        else:
            blocks = [block_nested(rng, "TE1 (outer)", "TE2 (inner)"), block_single(rng, "TE3")]
        rng.shuffle(blocks)
        if sum(len(b[0]) for b in blocks) + len(blocks) + 1 > total_len:
            continue
        seq, truth = assemble(blocks, total_len, rng)
        # the simulated boundaries must be the real (maximal) TIR boundaries
        if all(tir_is_maximal(seq, s, e) for _, s, e in truth):
            return seq, truth


# ---------------------------------------------------------------- detection

def detect_transposons(seq):
    found = []
    for start in range(len(seq)):
        left = seq[start:start + TIR_LEN]
        if len(left) < TIR_LEN:
            break
        target = reverse_complement(left)
        for length in range(DETECT_MIN, DETECT_MAX + 1):
            end = start + length - 1
            if end >= len(seq):
                break
            if seq[end - TIR_LEN + 1:end + 1] == target and tir_is_maximal(seq, start, end):
                found.append((start, end))
    return found


def classify(elements):
    labels = []
    for i, (s1, e1) in enumerate(elements):
        rel = []
        for j, (s2, e2) in enumerate(elements):
            if i == j:
                continue
            if s1 <= s2 and e2 <= e1:
                rel.append(f"contains #{j + 1}")
            elif s2 <= s1 and e1 <= e2:
                rel.append(f"nested in #{j + 1}")
            elif s1 <= e2 and s2 <= e1:
                rel.append(f"overlaps #{j + 1} ({min(e1, e2) - max(s1, s2) + 1} bp)")
        labels.append(", ".join(rel) if rel else "isolated")
    return labels


# ---------------------------------------------------------------- output

def print_case(title, seq, truth, detected):
    print("=" * 78)
    print(f"Case: {title}   |   length = {len(seq)} bp   |   simulated TEs = {len(truth)}")
    print("=" * 78)
    for i in range(0, len(seq), 60):
        print(f"{i + 1:4d}  {seq[i:i + 60]}")

    print("\nSimulated (ground truth):")
    for name, s, e in truth:
        print(f"  {name:<12} {s + 1:4d} - {e + 1:4d}  ({e - s + 1} bp)")

    print("\nDetected transposable elements:")
    print("   #  | Start |  End  | Len | Left TIR   | Right TIR  | Relation")
    for k, ((s, e), rel) in enumerate(zip(detected, classify(detected)), start=1):
        print(f"  {k:2d}  | {s + 1:5d} | {e + 1:5d} | {e - s + 1:3d} | "
              f"{seq[s:s + TIR_LEN]} | {seq[e - TIR_LEN + 1:e + 1]} | {rel}")

    truth_set = {(s, e) for _, s, e in truth}
    hits = len(truth_set & set(detected))
    print(f"\nCorrectly detected: {hits}/{len(truth_set)}   "
          f"False positives: {len(set(detected) - truth_set)}\n")


def draw(cases):
    fig, axes = plt.subplots(len(cases), 1, figsize=(12, 2.6 * len(cases)))
    for ax, (title, seq, truth, detected) in zip(axes, cases):
        ax.add_patch(Rectangle((0, -0.1), len(seq), 0.2, color="grey"))
        levels = []
        for k, (s, e) in enumerate(sorted(detected, key=lambda x: (x[0], -x[1]))):
            level = next((i for i, last in enumerate(levels) if last < s), len(levels))
            if level == len(levels):
                levels.append(e)
            levels[level] = e
            y = 0.35 + 0.5 * level
            ax.add_patch(Rectangle((s, y), e - s + 1, 0.35, color="#4c9be8", alpha=0.85))
            ax.add_patch(Rectangle((s, y), TIR_LEN, 0.35, color="#e8574c"))
            ax.add_patch(Rectangle((e - TIR_LEN + 1, y), TIR_LEN, 0.35, color="#e8574c"))
            ax.text((s + e) / 2, y + 0.17, f"{s + 1}-{e + 1}", ha="center", va="center",
                    fontsize=8, color="white")
        ax.set_xlim(0, len(seq))
        ax.set_ylim(-0.3, 0.45 + 0.5 * max(len(levels), 1) + 0.2)
        ax.set_yticks([])
        ax.set_title(f"{title} - {len(detected)} TEs detected (red = inverted repeats)")
    axes[-1].set_xlabel("Position (bp)")
    plt.tight_layout()
    plt.show()


if __name__ == "__main__":
    rng = random.Random(SEED)
    cases = []
    for kind in ("isolated", "overlapping", "nested"):
        seq, truth = build_case(kind, rng)
        detected = detect_transposons(seq)
        print_case(kind.upper(), seq, truth, detected)
        cases.append((kind.capitalize(), seq, truth, detected))
    draw(cases)
