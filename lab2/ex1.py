import math
import matplotlib.pyplot as plt
import numpy as np

def calculate_simple_tm(sequence):
    seq = sequence.upper()
    gc_count = seq.count('G') + seq.count('C')
    at_count = seq.count('A') + seq.count('T')
    return 4 * gc_count + 2 * at_count

def calculate_salt_tm(sequence, na_conc, constant=81.5):
 
    seq = sequence.upper()
    length = len(seq)
    if length == 0:
        return 0
    
    gc_count = seq.count('G') + seq.count('C')
    gc_percent = (gc_count / length) * 100
    
    tm = constant + 16.6 * math.log10(na_conc) + 0.41 * gc_percent - (600 / length)
    return tm

def main():
    S = "TAACGCTA"
    
    na_concentration = 0.05 
    
    simple_tm = calculate_simple_tm(S)
    salt_tm_81_5 = calculate_salt_tm(S, na_concentration, constant=81.5)
    salt_tm_8_15 = calculate_salt_tm(S, na_concentration, constant=8.15) 
    
    print("="*40)
    print(f"Sequence: {S}")
    print(f"Length: {len(S)} bp")
    print(f"Simple Tm: {simple_tm} °C")
    print(f"Salt-adjusted Tm (Na+ = {na_concentration}M, const=81.5): {salt_tm_81_5:.2f} °C")
    print(f"Salt-adjusted Tm (Na+ = {na_concentration}M, const=8.15): {salt_tm_8_15:.2f} °C")
    print("="*40)

   
    plt.style.use('bmh') 
    fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(12, 5))

    nucleotides = ['A', 'T', 'C', 'G']
    counts = [S.upper().count(n) for n in nucleotides]
    colors = ['#ff9999','#66b3ff','#99ff99','#ffcc99']
    
    ax1.bar(nucleotides, counts, color=colors, edgecolor='black')
    ax1.set_title(f'Nucleotide Composition (Seq: {S})')
    ax1.set_xlabel('Nucleotide')
    ax1.set_ylabel('Count')
    ax1.set_yticks(range(0, max(counts) + 2))

    na_concs = np.linspace(0.01, 1.0, 100)
    
    tm_values = [calculate_salt_tm(S, conc, constant=81.5) for conc in na_concs]
    
    ax2.plot(na_concs, tm_values, color='#8E44AD', linewidth=2.5)
    ax2.set_title('Salt Adjusted Tm vs Na+ Concentration')
    ax2.set_xlabel('[Na+] Concentration (M)')
    ax2.set_ylabel('Melting Temperature Tm (°C)')
    
    plt.tight_layout()
    plt.show()

if __name__ == "__main__":
    main()