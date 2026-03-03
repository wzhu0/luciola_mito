import os, subprocess

ASSEMBLY_DIR = "/home/wzhu/luciola/mito/01_assembly/output"
GFF_DIR      = "/home/wzhu/luciola/mito/02_annotation/output/mitos2"
OUT          = "/home/wzhu/luciola/mito/03_align/01_extracted/rrnL.fasta"

CORRECTED = {
    "AlVe-LuNov-34-m-2022-2":    ("327+(circular)",                       10976, 12226, False),
    "AlVe-LuNov-36-m-2022-4":    ("6282+(circular)",                       7862,  9112, False),
    "AlVe-LuNov-38-m-2022-6":    ("373+(circular)",                        7862,  9112, False),
    "MoDo-LuNov-31-m-2022-1":    ("712-(circular)",                       10723, 11973, False),
    "MoDo-LuNov-32-m-2022-2":    ("207-(circular)",                        4310,  5560, False),
    "ItPie-LuLus-20-m-2022-5":   ("189-(circular)",                       16215,   995, True),  # end=995 (21bp overlap)
    "PoAlf-LuLus-65-m-2022-11":  ("543_539_545_127_15539_1615-(circular)", 1851,  3101, False),
    "PoAlf-LuLus-66-m-2022-12":  ("6247+(circular)",                       1970,  3220, False),
    "PoAlf-LuLus-68-m-2022-14":  ("293+(circular)",                        6951,  8201, False),
}

GENOME_LEN = {
    "AlVe-LuNov-34-m-2022-2":   16470,
    "AlVe-LuNov-36-m-2022-4":   16471,
    "AlVe-LuNov-38-m-2022-6":   16471,
    "MoDo-LuNov-31-m-2022-1":   16469,
    "MoDo-LuNov-32-m-2022-2":   16469,
    "ItPie-LuLus-20-m-2022-5":  16470,
    "PoAlf-LuLus-65-m-2022-11": 16525,
    "PoAlf-LuLus-66-m-2022-12": 16525,
    "PoAlf-LuLus-68-m-2022-14": 16525,
}

def faidx(fa, contig, start, end):
    region = "{}:{}-{}".format(contig, start, end)
    r = subprocess.run(["samtools", "faidx", fa, region],
                       stdout=subprocess.PIPE, stderr=subprocess.PIPE)
    lines = [l for l in r.stdout.decode().strip().split("\n") if not l.startswith(">")]
    return "".join(lines)

def get_gff_rrnl(sample):
    gff = "{}/{}/result.gff".format(GFF_DIR, sample)
    with open(gff) as f:
        for line in f:
            if "\trRNA\t" in line and "rrnL" in line:
                parts = line.strip().split("\t")
                return parts[0], int(parts[3]), int(parts[4])
    raise ValueError("rrnL not found in GFF for {}".format(sample))

samples = sorted(os.listdir(GFF_DIR))

with open(OUT, "w") as out:
    for sample in samples:
        fa = "{}/{}/animal_mt.K105.complete.graph1.1.path_sequence.fasta".format(ASSEMBLY_DIR, sample)
        if not os.path.exists(fa):
            print("  SKIP (no assembly): {}".format(sample))
            continue

        if sample in CORRECTED:
            contig, start, end, wraps = CORRECTED[sample]
            genome_len = GENOME_LEN[sample]
            if wraps:
                seq = faidx(fa, contig, start, genome_len) + faidx(fa, contig, 1, end)
                print("  CORRECTED (wrap): {} {} bp".format(sample, len(seq)))
            else:
                seq = faidx(fa, contig, start, end)
                print("  CORRECTED: {} {} bp".format(sample, len(seq)))
        else:
            contig, start, end = get_gff_rrnl(sample)
            seq = faidx(fa, contig, start, end)
            print("  GFF: {} {} bp".format(sample, len(seq)))

        out.write(">{}\n{}\n".format(sample, seq))

print("\nDone. Written to {}".format(OUT))
