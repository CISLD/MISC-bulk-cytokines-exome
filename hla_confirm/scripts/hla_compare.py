#!/usr/bin/env python3
"""HLA class I: OptiType (exome) vs T1K (exome, bulk RNA flare and remission), two-field resolution,
with manuscript patient codes.

usage: hla_compare.py <project dir> <out.tsv>
"""
import csv
import os
import sys

P, out = sys.argv[1:3]
key = {r["wes_code"]: r for r in csv.DictReader(open(os.path.join(P, "wes", "patient_key.tsv")), delimiter="\t")}
opt = {r["code"]: r for r in csv.DictReader(open(os.path.join(P, "wes", "results", "hla_class_I.tsv")), delimiter="\t")}

def two(a):
    a = a.replace("HLA-", "")
    g, f = a.split("*")
    return g + "*" + ":".join(f.split(":")[:2])

def t1k(path):
    res = {}
    if not os.path.exists(path):
        return None
    for r in csv.reader(open(path), delimiter="\t"):
        gene = r[0].replace("HLA-", "")
        if gene in ("A", "B", "C"):
            res[gene] = sorted(two(r[i]) for i in (2, 5) if len(r) > i and r[i] not in ("", ".")) or []
            if len(res[gene]) == 1:
                res[gene] = res[gene] * 2
    return res

rows = []
for w in sorted(opt):
    o = opt[w]
    calls = {"OptiType_WES": {g: sorted([o[g + "1"], o[g + "2"]]) for g in "ABC"},
             "T1K_WES": t1k(os.path.join(P, "hla_confirm", "t1k", f"WES_{w}", f"WES_{w}_genotype.tsv")),
             "T1K_RNA_flare": t1k(os.path.join(P, "hla_confirm", "t1k", f"RNA_{w}_flare", f"RNA_{w}_flare_genotype.tsv")),
             "T1K_RNA_remission": t1k(os.path.join(P, "hla_confirm", "t1k", f"RNA_{w}_remission", f"RNA_{w}_remission_genotype.tsv"))}
    mc = key.get(w, {}).get("manuscript_code", "?")
    for g in "ABC":
        ref = calls["OptiType_WES"][g]
        r = [mc if mc != "-" else "outside main cohort", w, g]
        agree = []
        for k in ("OptiType_WES", "T1K_WES", "T1K_RNA_flare", "T1K_RNA_remission"):
            c = calls[k]
            if c is None or g not in c:
                r.append("NA")
                continue
            r.append(" / ".join(c[g]))
            if k != "OptiType_WES":
                agree.append(c[g] == ref)
        r.append("all agree" if agree and all(agree) else ("differs" if agree else "NA"))
        rows.append(r)
with open(out, "w", newline="") as fh:
    w_ = csv.writer(fh, delimiter="\t", lineterminator="\n")
    w_.writerow(["patient", "wes_code", "locus", "OptiType_WES", "T1K_WES", "T1K_RNA_flare", "T1K_RNA_remission", "concordance"])
    w_.writerows(rows)
for r in rows:
    if r[-1] != "all agree":
        print("\t".join(r))
print("loci:", len(rows), "| all agree:", sum(r[-1] == "all agree" for r in rows))
