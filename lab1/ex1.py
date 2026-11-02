S = "ATTGTCCCAATCTGTTG"


def find_alphabet(sequence):
    alphabet = []
    for symbol in sequence:
        if symbol not in alphabet:
            alphabet.append(symbol)
    return sorted(alphabet)


if __name__ == "__main__":
    alphabet = find_alphabet(S)
    print(f"Sequence: {S}")
    print(f"Length  : {len(S)}")
    print(f"Alphabet: {{{', '.join(alphabet)}}}  ({len(alphabet)} symbols)")
