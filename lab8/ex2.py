"""Lab 8 - Ex 2: promoter pattern detection with a sliding window.

For every window of 30 bp the application computes:
    - C+G %                 = (count(C) + count(G)) / window length * 100
    - Kappa Index of Coincidence (IC): the sequence is compared with copies of
      itself shifted by u = 1..N-1 positions; for each shift the percentage of
      identical positions is computed, and IC is the mean of these percentages.

The two signals form the "pattern" of the sequence. The center of weight of each
signal is computed and drawn on the chart.

Usage: python ex2.py [sequence] [window]
"""
import sys

import matplotlib.pyplot as plt

S = (sys.argv[1].upper() if len(sys.argv) > 1 else
     "CGGACTGATCTATCT" + "AC" * 27 + "CGTAGCATCTATCGATCTATCTAGCGATCTA")
WINDOW = int(sys.argv[2]) if len(sys.argv) > 2 else 30


def cg_content(seq):
    return (seq.count("C") + seq.count("G")) / len(seq) * 100


def kappa_ic(seq):
    n = len(seq)
    total = 0.0
    for u in range(1, n):
        matches = sum(1 for i in range(n - u) if seq[i] == seq[i + u])
        total += matches / (n - u) * 100
    return total / (n - 1)


def center_of_weight(values):
    total = sum(values)
    return sum(i * v for i, v in enumerate(values)) / total if total else 0.0


if __name__ == "__main__":
    if len(S) < WINDOW:
        sys.exit(f"Sequence is shorter than the window ({len(S)} < {WINDOW})")

    windows = [S[i:i + WINDOW] for i in range(len(S) - WINDOW + 1)]
    cg = [cg_content(w) for w in windows]
    ic = [kappa_ic(w) for w in windows]

    print(f"Sequence: {S}")
    print(f"Length  : {len(S)} nt   Window: {WINDOW} nt   Windows: {len(windows)}")
    print(f"Whole sequence: C+G = {cg_content(S):.2f}%   IC = {kappa_ic(S):.2f}\n")

    print("  Pos | Window                          |  C+G %  |   IC")
    for i, (w, c, k) in enumerate(zip(windows, cg, ic), start=1):
        print(f"  {i:3d} | {w} | {c:6.2f}  | {k:6.2f}")

    cg_center, ic_center = center_of_weight(cg), center_of_weight(ic)
    print(f"\nCenter of weight  C+G: {cg_center + 1:.2f}   IC: {ic_center + 1:.2f}")

    x = range(1, len(windows) + 1)
    fig, (ax1, ax2) = plt.subplots(2, 1, figsize=(12, 8))

    ax1.plot(x, cg, color="tab:blue", label="C+G %")
    ax1.plot(x, ic, color="tab:red", label="Kappa IC")
    ax1.set_title(f"Promoter pattern (sliding window = {WINDOW} nt)")
    ax1.set_xlabel("Window start position")
    ax1.set_ylabel("Value")
    ax1.legend(); ax1.grid(alpha=0.3)

    ax2.scatter(cg, ic, c=list(x), cmap="viridis", s=18, label="windows")
    ax2.scatter([sum(cg) / len(cg)], [sum(ic) / len(ic)], color="red", s=120, marker="X",
                label=f"center of weight (C+G={sum(cg) / len(cg):.2f}, IC={sum(ic) / len(ic):.2f})")
    ax2.set_title("Pattern in the (C+G %, IC) plane")
    ax2.set_xlabel("C+G %")
    ax2.set_ylabel("Kappa IC")
    ax2.legend(); ax2.grid(alpha=0.3)

    for ax, center, color, name in ((ax1, cg_center, "tab:blue", "C+G"),
                                    (ax1, ic_center, "tab:red", "IC")):
        ax.axvline(center + 1, color=color, linestyle="--", alpha=0.7,
                   label=f"{name} center = {center + 1:.1f}")
    ax1.legend()

    plt.tight_layout()
    plt.show()
