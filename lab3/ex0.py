def translate_codon(codon):
    
    genetic_code = {
        'ATA':'I', 'ATC':'I', 'ATT':'I', 'ATG':'M',
        'ACA':'T', 'ACC':'T', 'ACG':'T', 'ACT':'T',
        'AAC':'N', 'AAT':'N', 'AAA':'K', 'AAG':'K',
        'AGC':'S', 'AGT':'S', 'AGA':'R', 'AGG':'R',
        'CTA':'L', 'CTC':'L', 'CTG':'L', 'CTT':'L',
        'CCA':'P', 'CCC':'P', 'CCG':'P', 'CCT':'P',
        'CAC':'H', 'CAT':'H', 'CAA':'Q', 'CAG':'Q',
        'CGA':'R', 'CGC':'R', 'CGG':'R', 'CGT':'R',
        'GTA':'V', 'GTC':'V', 'GTG':'V', 'GTT':'V',
        'GCA':'A', 'GCC':'A', 'GCG':'A', 'GCT':'A',
        'GAC':'D', 'GAT':'D', 'GAA':'E', 'GAG':'E',
        'GGA':'G', 'GGC':'G', 'GGG':'G', 'GGT':'G',
        'TCA':'S', 'TCC':'S', 'TCG':'S', 'TCT':'S',
        'TTC':'F', 'TTT':'F', 'TTA':'L', 'TTG':'L',
        'TAC':'Y', 'TAT':'Y', 'TAA':'_', 'TAG':'_',
        'TGC':'C', 'TGT':'C', 'TGA':'_', 'TGG':'W',
    }
    
    # Normalize input: make uppercase and convert RNA 'U' to DNA 'T'
    codon = codon.upper().replace('U', 'T')
    
    # Return the amino acid, or '?' if the codon is invalid
    return genetic_code.get(codon, '?')

def translate_sequence(seq, frame=0):
    """Translate a DNA/RNA sequence codon by codon, starting at the given frame (0, 1, 2)."""
    seq = seq.upper().replace('U', 'T')
    return ''.join(translate_codon(seq[i:i + 3]) for i in range(frame, len(seq) - 2, 3))


if __name__ == "__main__":
    # Example usage:
    print(f"AUG translates to: {translate_codon('AUG')}") # Output: M (Methionine)
    # Example from the board: AAU CCC -> N (Asparagine), P (Proline)
    print(f"AAUCCC translates to: {translate_sequence('AAUCCC')}")