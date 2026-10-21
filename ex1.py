import urllib.request

ACCESSION = "NC_001653.2"
OUT_FILE = "hdv_genome.fasta"
URL = ("https://eutils.ncbi.nlm.nih.gov/entrez/eutils/efetch.fcgi"
       f"?db=nuccore&id={ACCESSION}&rettype=fasta&retmode=text")


def read_fasta(path):
    with open(path) as f:
        lines = f.read().splitlines()
    header = lines[0][1:]
    seq = ''.join(line.strip() for line in lines[1:])
    return header, seq.upper()


if __name__ == "__main__":
    with urllib.request.urlopen(URL) as resp:
        data = resp.read().decode()
    with open(OUT_FILE, "w") as f:
        f.write(data)
    header, seq = read_fasta(OUT_FILE)
    print(f"Downloaded: {header}")
    print(f"Length: {len(seq)} nt")
