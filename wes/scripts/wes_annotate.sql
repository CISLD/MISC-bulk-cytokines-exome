-- Join panel variants to the in-house GRCh37 frequency database (8,814 cases) and ClinVar 20260615
-- (lifted to GRCh37). Paths are substituted by wes_scan.sb.
CREATE TABLE pv AS
SELECT * FROM read_csv('@IN@', delim = '\t', header = true, all_varchar = true);

ATTACH '@SNV@' AS snv (READ_ONLY);
ATTACH '@CV@'  AS cv  (READ_ONLY);

CREATE TABLE c38 AS
SELECT variation_id, any_value(ref) AS ref38, any_value(alt) AS alt38
FROM cv.clinvar_grch38 GROUP BY variation_id;

COPY (
  SELECT pv.*,
         s.obs AS inhouse_cases, s.hom AS inhouse_hom,
         c.variation_id AS clinvar_id, c.cls AS clinvar_class, c.clnsig AS clinvar_sig,
         c.stars AS clinvar_stars, c.clnrevstat AS clinvar_review,
         -- a RECOVER_SWAPPED_REF_ALT record: GRCh37 ALT is the GRCh38 REF (same or opposite strand)
         coalesce((k.ref38 = c.alt AND k.alt38 = c.ref)
               OR (translate(k.ref38, 'ACGT', 'TGCA') = c.alt AND translate(k.alt38, 'ACGT', 'TGCA') = c.ref),
                  false) AS clinvar_swapped
  FROM pv
  LEFT JOIN snv.variants s
         ON s.chrom = pv.chrom AND s.pos = CAST(pv.pos AS BIGINT) AND s.ref = pv.ref AND s.alt = pv.alt
  LEFT JOIN cv.clinvar_grch37 c
         ON c.chrom = pv.chrom AND c.pos = CAST(pv.pos AS BIGINT) AND c.ref = pv.ref AND c.alt = pv.alt
  LEFT JOIN c38 k ON k.variation_id = c.variation_id
  ORDER BY pv.code, pv.chrom, CAST(pv.pos AS BIGINT)
) TO '@OUT@' (HEADER, DELIMITER '\t');
