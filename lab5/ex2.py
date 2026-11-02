from collections import Counter
from itertools import product

import matplotlib.pyplot as plt

from ex1 import read_fasta

ALPHABET = "ACGT"
UNIT_MIN, UNIT_MAX = 3, 6
REPEAT_MIN, REPEAT_MAX = 6, 10


def brute_force_patterns():
    for size in range(UNIT_MIN, UNIT_MAX + 1):
        for combo in product(ALPHABET, repeat=size):
            yield "".join(combo)


def find_tandem_repeats(seq, unit):
    found = []
    i = seq.find(unit * 2)
    while i != -1:
        copies = 2
        while seq.startswith(unit, i + copies * len(unit)):
            copies += 1
        found.append((i, copies))
        i = seq.find(unit * 2, i + copies * len(unit))
    return found


def is_primitive(unit):
    return (unit + unit).find(unit, 1) == len(unit)


def detect_repetitions(seq):
    results = []
    for unit in brute_force_patterns():
        if not is_primitive(unit):
            continue
        for pos, copies in find_tandem_repeats(seq, unit):
            total = copies * len(unit)
            if REPEAT_MIN <= total <= REPEAT_MAX:
                results.append((pos, unit, copies, seq[pos:pos + total]))
    return sorted(results)


if __name__ == "__main__":
    header, seq = read_fasta()
    print(f"Sequence: {header}")
    print(f"Length  : {len(seq)} nt")
    print(f"Patterns tested (brute force, {UNIT_MIN}-{UNIT_MAX} letters): "
          f"{sum(1 for _ in brute_force_patterns())}\n")

    reps = detect_repetitions(seq)
    print(f"Repetitions of {REPEAT_MIN}-{REPEAT_MAX} letters found: {len(reps)}\n")
    print("Position | Unit   | Copies | Repetition")
    for pos, unit, copies, region in reps:
        print(f"{pos:8d} | {unit:<6} | {copies:6d} | {region}")

    if reps:
        per_unit = Counter(unit for _, unit, _, _ in reps).most_common(15)
        plt.figure(figsize=(11, 5))
        plt.bar([u for u, _ in per_unit], [c for _, c in per_unit], color="steelblue")
        plt.title(f"Most frequent tandem repeats ({REPEAT_MIN}-{REPEAT_MAX} nt) - {header}")
        plt.xlabel("Repeated unit")
        plt.ylabel("Occurrences")
        plt.xticks(rotation=45)
        plt.tight_layout()
        plt.show()
