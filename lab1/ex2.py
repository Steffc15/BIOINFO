import matplotlib.pyplot as plt

from ex1 import S, find_alphabet


def relative_frequencies(sequence):
    alphabet = find_alphabet(sequence)
    return {symbol: sequence.count(symbol) / len(sequence) for symbol in alphabet}


if __name__ == "__main__":
    freqs = relative_frequencies(S)

    print(f"Sequence: {S}  ({len(S)} symbols)\n")
    print("Symbol | Count | Relative frequency")
    for symbol, f in freqs.items():
        print(f"  {symbol}    |  {S.count(symbol):3d}  | {f:.4f}  ({f:.2%})")
    print(f"\nSum of frequencies: {sum(freqs.values()):.2f}")

    plt.figure(figsize=(6, 4))
    plt.bar(freqs.keys(), [f * 100 for f in freqs.values()],
            color=['#ff9999', '#66b3ff', '#99ff99', '#ffcc99'][:len(freqs)],
            edgecolor='black')
    plt.title(f"Relative frequencies in {S}")
    plt.xlabel("Symbol")
    plt.ylabel("Percentage (%)")
    plt.tight_layout()
    plt.show()
