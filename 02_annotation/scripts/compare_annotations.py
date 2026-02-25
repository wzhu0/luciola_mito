#!/usr/bin/env python3
import os

PCG_MITOS2 = {'nad1','nad2','nad3','nad4','nad4l','nad5','nad6',
              'cox1','cox2','cox3','atp6','atp8','cob'}  # cob not cytb!

PCG_MITOFINDER = {'ND1','ND2','ND3','ND4','ND4L','ND5','ND6',
                  'COX1','COX2','COX3','ATP6','ATP8','CYTB'}

NAME_MAP = {
    'nad1':'ND1','nad2':'ND2','nad3':'ND3','nad4':'ND4','nad4l':'ND4L',
    'nad5':'ND5','nad6':'ND6','cox1':'COX1','cox2':'COX2','cox3':'COX3',
    'atp6':'ATP6','atp8':'ATP8','cob':'CYTB'
}

WORK_DIR = os.path.expanduser("~/luciola/mito")
SAMPLE_LIST = os.path.join(WORK_DIR, "utils/sample_list.txt")
MITOS2_DIR = os.path.join(WORK_DIR, "02_annotation/output/mitos2")
MF_DIR = os.path.join(WORK_DIR, "02_annotation/output/mitofinder")

def parse_mitos2_gff(gff_file):
    genes = {}
    with open(gff_file) as f:
        for line in f:
            if line.startswith('#'): continue
            parts = line.strip().split('\t')
            if len(parts) < 9: continue
            if parts[2] != 'gene': continue
            attrs = parts[8]
            name = None
            for attr in attrs.split(';'):
                if attr.startswith('Name='):
                    name = attr.replace('Name=','').strip()
                    break
            if name and name.lower() in PCG_MITOS2:
                genes[name.lower()] = (int(parts[3]), int(parts[4]), parts[6])
    return genes

def parse_mitofinder_gff(gff_file):
    genes = {}
    with open(gff_file) as f:
        for line in f:
            if line.startswith('#'): continue
            parts = line.strip().split('\t')
            if len(parts) < 9: continue
            if parts[2] != 'gene': continue
            name_field = parts[8]
            if not name_field.startswith('Name='):
                continue
            # name is first word after Name=, before ' gene'
            name = name_field.replace('Name=','').split()[0].strip()
            if name in PCG_MITOFINDER:
                # keep only first occurrence (avoid partials overwriting full)
                if name not in genes:
                    genes[name] = (int(parts[3]), int(parts[4]), parts[6])
    return genes

with open(SAMPLE_LIST) as f:
    samples = [os.path.basename(l.strip()) for l in f if l.strip()]

print(f"{'SAMPLE':<45} {'GENE':<8} {'M2_START':>9} {'M2_END':>9} {'M2_STR':>6} {'MF_START':>9} {'MF_END':>9} {'MF_STR':>6} {'ΔSTART':>8} {'ΔEND':>8} {'MATCH':>6}")
print("-" * 130)

for sample in samples:
    m2_gff = os.path.join(MITOS2_DIR, sample, "result.gff")
    mf_gff = os.path.join(MF_DIR, sample,
                f"{sample}_MitoFinder_arwen_Final_Results",
                f"{sample}_mtDNA_contig.gff")

    if not os.path.exists(m2_gff) or not os.path.exists(mf_gff):
        print(f"{sample:<45} MISSING FILES")
        continue

    m2_genes = parse_mitos2_gff(m2_gff)
    mf_genes = parse_mitofinder_gff(mf_gff)

    for m2_name in sorted(PCG_MITOS2):
        mf_name = NAME_MAP[m2_name]
        m2 = m2_genes.get(m2_name)
        mf = mf_genes.get(mf_name)

        if m2 and mf:
            d_start = mf[0] - m2[0]
            d_end = mf[1] - m2[1]
            match = "OK" if abs(d_start) <= 3 and abs(d_end) <= 3 else "DIFF"
            print(f"{sample:<45} {mf_name:<8} {m2[0]:>9} {m2[1]:>9} {m2[2]:>6} {mf[0]:>9} {mf[1]:>9} {mf[2]:>6} {d_start:>8} {d_end:>8} {match:>6}")
        elif m2 and not mf:
            print(f"{sample:<45} {mf_name:<8} {m2[0]:>9} {m2[1]:>9} {m2[2]:>6} {'MISSING':>9} {'':>9} {'':>6} {'':>8} {'':>8} {'MISS_MF':>6}")
        elif mf and not m2:
            print(f"{sample:<45} {mf_name:<8} {'MISSING':>9} {'':>9} {'':>6} {mf[0]:>9} {mf[1]:>9} {mf[2]:>6} {'':>8} {'':>8} {'MISS_M2':>6}")
