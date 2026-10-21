import os
import random
import sys
import time
import urllib.request

ACCESSION = sys.argv[1] if len(sys.argv) > 1 else "NC_045512.2"
N_READS = int(sys.argv[2]) if len(sys.argv) > 2 else 2000
READ_MIN, READ_MAX = 100, 150
MIN_OVERLAP = 30
SEED = 42

HERE = os.path.dirname(os.path.abspath(__file__))
FASTA = os.path.join(HERE, f"{ACCESSION}.fasta")
URL = ("https://eutils.ncbi.nlm.nih.gov/entrez/eutils/efetch.fcgi"
       f"?db=nuccore&id={ACCESSION}&rettype=fasta&retmode=text")


def download_genome():
    if not os.path.exists(FASTA):
        with urllib.request.urlopen(URL) as resp:
            data = resp.read().decode()
        with open(FASTA, "w") as f:
            f.write(data)
    with open(FASTA) as f:
        lines = f.read().splitlines()
    return lines[0][1:], "".join(l.strip() for l in lines[1:]).upper()


def sample_reads(genome, n, rng):
    reads = []
    for _ in range(n):
        length = rng.randint(READ_MIN, READ_MAX)
        start = rng.randint(0, len(genome) - length)
        reads.append(genome[start:start + length])
    return reads


def remove_contained(reads):
    uniq = list(dict.fromkeys(reads))
    joined = "$".join(uniq)
    starts, pos = [], 0
    for r in uniq:
        starts.append(pos)
        pos += len(r) + 1
    kept = []
    for r, own in zip(uniq, starts):
        p = joined.find(r)
        while p != -1 and p == own:
            p = joined.find(r, p + 1)
        if p == -1:
            kept.append(r)
    return kept


def all_overlaps(reads):
    k = MIN_OVERLAP
    prefix_index = {}
    for j, r in enumerate(reads):
        prefix_index.setdefault(r[:k], []).append(j)
    overlaps = []
    for i, a in enumerate(reads):
        for start in range(1, len(a) - k + 1):
            for j in prefix_index.get(a[start:start + k], ()):
                if j != i:
                    olen = len(a) - start
                    if reads[j][:olen] == a[start:]:
                        overlaps.append((olen, i, j))
    return overlaps


def greedy_assemble(reads):
    reads = remove_contained(reads)
    overlaps = all_overlaps(reads)
    overlaps.sort(reverse=True)

    n = len(reads)
    succ, pred = [None] * n, [None] * n
    group = list(range(n))

    def find(x):
        while group[x] != x:
            group[x] = group[group[x]]
            x = group[x]
        return x

    for olen, i, j in overlaps:
        if succ[i] is None and pred[j] is None and find(i) != find(j):
            succ[i], pred[j] = (j, olen), i
            group[find(i)] = find(j)

    contigs = []
    for start in range(n):
        if pred[start] is None:
            seq, cur = [reads[start]], start
            while succ[cur] is not None:
                nxt, olen = succ[cur]
                seq.append(reads[nxt][olen:])
                cur = nxt
            contigs.append("".join(seq))
    contigs.sort(key=len, reverse=True)
    return contigs, n, len(overlaps)


def uncovered_regions(genome, reads):
    cov = [0] * (len(genome) + 1)
    for r in reads:
        p = genome.find(r)
        cov[p] += 1
        cov[p + len(r)] -= 1
    gaps, depth, gap_start = [], 0, None
    for i in range(len(genome)):
        depth += cov[i]
        if depth == 0 and gap_start is None:
            gap_start = i
        elif depth > 0 and gap_start is not None:
            gaps.append((gap_start, i))
            gap_start = None
    if gap_start is not None:
        gaps.append((gap_start, len(genome)))
    return gaps


def longest_repeat(genome, max_k=5000):
    def has_repeat(k):
        seen = set()
        for i in range(len(genome) - k + 1):
            w = genome[i:i + k]
            if w in seen:
                return True
            seen.add(w)
        return False
    lo, hi = 0, min(max_k, len(genome))
    while lo < hi:
        mid = (lo + hi + 1) // 2
        if has_repeat(mid):
            lo = mid
        else:
            hi = mid - 1
    return lo


def main():
    rng = random.Random(SEED)
    header, genome = download_genome()
    print(f"Genome : {header}")
    print(f"Length : {len(genome)} nt")

    reads = sample_reads(genome, N_READS, rng)
    avg = sum(map(len, reads)) / len(reads)
    coverage = sum(map(len, reads)) / len(genome)
    print(f"Reads  : {len(reads)} samples, {READ_MIN}-{READ_MAX} nt "
          f"(avg {avg:.1f}), coverage {coverage:.1f}x")
    print(f"First sample: {reads[0][:60]}...")

    t0 = time.time()
    contigs, n_kept, n_ov = greedy_assemble(reads)
    dt = time.time() - t0

    print("\n--- Reconstruction (greedy overlap) ---")
    print(f"Reads after removing contained/duplicate ones: {n_kept}")
    print(f"Overlaps >= {MIN_OVERLAP} nt found: {n_ov}")
    print(f"Contigs produced: {len(contigs)}   (1 = full genome rebuilt)")
    print(f"Time: {dt:.2f} s")

    correct = [c for c in contigs if c in genome]
    wrong = [c for c in contigs if c not in genome]
    print(f"Longest contig: {len(contigs[0])} nt "
          f"({100 * len(contigs[0]) / len(genome):.1f}% of genome)")
    print(f"Contigs that match the genome exactly: {len(correct)}")
    print(f"Mis-assembled contigs (wrong order, chimeras): {len(wrong)}")
    print(f"Total assembled length: {sum(map(len, contigs))} nt "
          f"vs genome {len(genome)} nt")
    if len(contigs) == 1 and contigs[0] == genome:
        print("=> The genome was rebuilt PERFECTLY.")

    gaps = uncovered_regions(genome, reads)
    print("\n--- Why it is (not) perfect ---")
    print(f"Genome regions never sampled by any read: {len(gaps)} "
          f"({sum(b - a for a, b in gaps)} nt missing)")
    for a, b in gaps[:5]:
        print(f"   gap {a}-{b} ({b - a} nt)")
    rep = longest_repeat(genome)
    print(f"Longest repeated region in genome: {rep} nt "
          f"(reads are {READ_MIN}-{READ_MAX} nt)")
    if rep >= READ_MIN:
        print("   -> repeats longer than a read: the assembler cannot know "
              "how the copies connect, so it collapses or mis-joins them.")


if __name__ == "__main__":
    main()
