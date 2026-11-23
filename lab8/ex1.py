"""Lab 8 - Ex 1: global alignment of two DNA sequences (Needleman-Wunsch).

Builds the score matrix, reconstructs the optimal alignment by traceback and
shows the matrix with the traceback path.

Usage: python ex1.py [S1] [S2] [match] [mismatch] [gap]
"""
import sys

import matplotlib.pyplot as plt
import numpy as np

S1 = sys.argv[1].upper() if len(sys.argv) > 1 else "ACCGTGAAGCCAATAC"
S2 = sys.argv[2].upper() if len(sys.argv) > 2 else "AGCGTGCAGCCAATAC"
MATCH = int(sys.argv[3]) if len(sys.argv) > 3 else 1
MISMATCH = int(sys.argv[4]) if len(sys.argv) > 4 else -1
GAP = int(sys.argv[5]) if len(sys.argv) > 5 else -2


def score(a, b):
    return MATCH if a == b else MISMATCH


def needleman_wunsch(s1, s2):
    rows, cols = len(s1) + 1, len(s2) + 1
    m = [[0] * cols for _ in range(rows)]
    for i in range(rows):
        m[i][0] = i * GAP
    for j in range(cols):
        m[0][j] = j * GAP

    for i in range(1, rows):
        for j in range(1, cols):
            m[i][j] = max(m[i - 1][j - 1] + score(s1[i - 1], s2[j - 1]),
                          m[i - 1][j] + GAP,
                          m[i][j - 1] + GAP)
    return m


def traceback(m, s1, s2):
    i, j = len(s1), len(s2)
    a1, a2, path = [], [], [(i, j)]
    while i > 0 or j > 0:
        if i > 0 and j > 0 and m[i][j] == m[i - 1][j - 1] + score(s1[i - 1], s2[j - 1]):
            a1.append(s1[i - 1]); a2.append(s2[j - 1]); i -= 1; j -= 1
        elif i > 0 and m[i][j] == m[i - 1][j] + GAP:
            a1.append(s1[i - 1]); a2.append("-"); i -= 1
        else:
            a1.append("-"); a2.append(s2[j - 1]); j -= 1
        path.append((i, j))
    return "".join(reversed(a1)), "".join(reversed(a2)), path[::-1]


def print_matrix(m, s1, s2):
    print("        " + "".join(f"{c:>4}" for c in "-" + s2))
    for i, row in enumerate(m):
        print(f"   {('-' + s1)[i]:>4} " + "".join(f"{v:4d}" for v in row))


def draw(m, s1, s2, path, a1, a2, bars):
    fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(14, 7),
                                   gridspec_kw={"width_ratios": [1.15, 1]})
    data = np.array(m)
    im = ax1.imshow(data, cmap="viridis")
    fig.colorbar(im, ax=ax1, fraction=0.046, pad=0.04, label="Score")
    ax1.set_xticks(range(len(s2) + 1))
    ax1.set_xticklabels("-" + s2)
    ax1.set_yticks(range(len(s1) + 1))
    ax1.set_yticklabels("-" + s1)
    ax1.xaxis.tick_top()
    ax1.set_xlabel("S2"); ax1.set_ylabel("S1")
    ax1.set_title("Score matrix", pad=25)

    path_grid = np.zeros_like(data, dtype=float)
    for i, j in path:
        path_grid[i, j] = 1
    ax2.imshow(path_grid, cmap="Reds", vmin=0, vmax=1.3)
    ax2.set_xticks(np.arange(-0.5, len(s2) + 1, 1), minor=True)
    ax2.set_yticks(np.arange(-0.5, len(s1) + 1, 1), minor=True)
    ax2.grid(which="minor", color="lightgrey", linewidth=0.5)
    ax2.plot([j for _, j in path], [i for i, _ in path], color="black", linewidth=1.5)
    ax2.set_xticks(range(len(s2) + 1))
    ax2.set_xticklabels("-" + s2)
    ax2.set_yticks(range(len(s1) + 1))
    ax2.set_yticklabels("-" + s1)
    ax2.xaxis.tick_top()
    ax2.set_title("Traceback path", pad=25)

    fig.text(0.5, 0.06, f"{a1}\n{bars}\n{a2}", ha="center", family="monospace", fontsize=13)
    fig.suptitle(f"Needleman-Wunsch  (match={MATCH}, mismatch={MISMATCH}, gap={GAP})  "
                 f"score = {m[-1][-1]}")
    plt.tight_layout(rect=(0, 0.13, 1, 0.95))
    plt.show()


if __name__ == "__main__":
    m = needleman_wunsch(S1, S2)
    a1, a2, path = traceback(m, S1, S2)
    bars = "".join("|" if x == y else (" " if "-" in (x, y) else ".") for x, y in zip(a1, a2))
    matches = bars.count("|")

    print(f"S1: {S1}  ({len(S1)} nt)")
    print(f"S2: {S2}  ({len(S2)} nt)")
    print(f"Scoring: match = {MATCH}, mismatch = {MISMATCH}, gap = {GAP}\n")
    print("Score matrix:")
    print_matrix(m, S1, S2)

    print("\nOptimal global alignment:")
    print(f"  {a1}\n  {bars}\n  {a2}")
    print(f"\nMatches   : {matches}")
    print(f"Mismatches: {bars.count('.')}")
    print(f"Gaps      : {a1.count('-') + a2.count('-')}")
    print(f"Length    : {len(a1)}")
    print(f"Similarity: {matches / len(a1):.0%}")
    print(f"Score     : {m[-1][-1]}")

    draw(m, S1, S2, path, a1, a2, bars)
