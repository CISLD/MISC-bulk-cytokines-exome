#!/usr/bin/env python3
"""Per-patient coding/splice variants in the Table S9 gene panel, from bcbio GATK VCFs
(SnpEff 4.3t ANN, gnomAD exomes AF) after `bcftools norm -m -any` and PASS/alt filtering.

Input per sample: TSV from
  bcftools query -f '%CHROM\t%POS\t%REF\t%ALT\t[%GT]\t[%AD]\t[%DP]\t[%GQ]\t%INFO/ANN\t%INFO/gnomAD_AF\t%INFO/af_exac_all\t%INFO/clinvar_sig\n'

usage: wes_panel_variants.py <panel_genes.tsv> <out.tsv> <code>=<query.tsv.gz> ...
"""
import csv, gzip, sys

LOF = {'stop_gained', 'frameshift_variant', 'splice_acceptor_variant', 'splice_donor_variant',
       'start_lost', 'stop_lost', 'exon_loss_variant', 'transcript_ablation'}
MISSENSE = {'missense_variant', 'inframe_insertion', 'inframe_deletion', 'disruptive_inframe_insertion',
            'disruptive_inframe_deletion', 'protein_protein_contact', 'initiator_codon_variant'}
SPLICE_REGION = {'splice_region_variant'}
RANK = {**{t: 0 for t in LOF}, **{t: 1 for t in MISSENSE}, **{t: 2 for t in SPLICE_REGION}}
CLASS = {0: 'LoF', 1: 'missense/inframe', 2: 'splice_region'}


def fnum(x):
    try:
        return float(x)
    except (TypeError, ValueError):
        return None


def main():
    panel_path, out_path, *inputs = sys.argv[1:]
    panel = {}
    for r in csv.reader(open(panel_path), delimiter='\t'):
        panel.setdefault(r[0], set()).add(r[1])
    seen_genes = set()
    out = csv.writer(open(out_path, 'w'), delimiter='\t', lineterminator='\n')
    out.writerow(['code', 'chrom', 'pos', 'ref', 'alt', 'gt', 'zygosity', 'dp', 'gq', 'alt_fraction',
                  'qual_pass', 'gene', 'panel_category', 'effect', 'variant_class', 'transcript',
                  'hgvs_c', 'hgvs_p', 'gnomad_exome_af', 'exac_af', 'clinvar_bcbio'])
    for spec in inputs:
        code, path = spec.split('=', 1)
        for line in gzip.open(path, 'rt'):
            f = line.rstrip('\n').split('\t')
            chrom, pos, ref, alt, gt, ad, dp, gq, ann, gaf, eaf, cv = f[:12]
            best = None
            for e in ann.split(','):
                a = e.split('|')
                if len(a) < 11:
                    continue
                gene = a[3]
                if gene in panel:
                    seen_genes.add(gene)
                else:
                    continue
                effs = a[1].split('&')
                r = min((RANK[t] for t in effs if t in RANK), default=None)
                if r is None:
                    continue
                # prefer the allele-matched entry, then the most severe class
                key = (0 if a[0] == alt else 1, r)
                if best is None or key < best[0]:
                    best = (key, gene, a[1], r, a[6], a[9], a[10])
            if best is None:
                continue
            _, gene, eff, r, tx, hc, hp = best
            alleles = gt.replace('|', '/').split('/')
            zyg = 'hom' if alleles.count('1') == 2 else ('het' if '1' in alleles else 'other')
            ads = [fnum(x) for x in ad.split(',')]
            af = None
            if len(ads) >= 2 and None not in ads[:2] and sum(ads[:2]) > 0:
                af = ads[1] / sum(ads[:2])
            dpv, gqv = fnum(dp), fnum(gq)
            ok = (dpv is not None and dpv >= 10 and gqv is not None and gqv >= 20 and
                  (zyg != 'het' or (af is not None and 0.2 <= af <= 0.8)))
            out.writerow([code, chrom, pos, ref, alt, gt, zyg, dp, gq,
                          '' if af is None else f'{af:.2f}', ok, gene, ';'.join(sorted(panel[gene])),
                          eff, CLASS[r], tx, hc, hp, gaf if gaf != '.' else '', eaf if eaf != '.' else '',
                          cv if cv != '.' else ''])
    missing = sorted(set(panel) - seen_genes)
    with open(out_path.replace('.tsv', '_genes_without_annotation.txt'), 'w') as fh:
        fh.write('\n'.join(missing) + '\n')
    print(f'panel genes never seen in any ANN entry: {len(missing)}')


if __name__ == '__main__':
    main()
