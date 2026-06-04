import argparse
import sys
from collections import defaultdict, Counter
import pysam
import pandas as pd

__version__ = "0.1.1"

PLUS_STRAND_AMBIGUOUS  = frozenset({("G", "A")})
MINUS_STRAND_AMBIGUOUS = frozenset({("C", "T")})
ANY_STRAND_AMBIGUOUS   = PLUS_STRAND_AMBIGUOUS | MINUS_STRAND_AMBIGUOUS

def _ambiguous_set_for(strand):
    if strand == "+":
        return PLUS_STRAND_AMBIGUOUS
    if strand == "-":
        return MINUS_STRAND_AMBIGUOUS
    return ANY_STRAND_AMBIGUOUS

def is_conversion_ambiguous(vtype, ref, alt, strand=None):
    r, a = ref.upper(), alt.upper()
    rules = _ambiguous_set_for(strand)
    if vtype == "snv":
        return (r, a) in rules
    if r and a and r[0] != a[0]:
        return (r[0], a[0]) in rules
    return False

def parse_variants(path):
    sep = "," if path.lower().endswith(".csv") else "\t"
    df = pd.read_csv(path, sep=sep, dtype={"chrom": str}, comment="#")
    required = {"name", "chrom", "pos", "ref", "alt", "type"}
    missing = required - set(df.columns)
    if missing:
        kind = "CSV" if sep == "," else "TSV"
        sys.exit(f"Variants {kind} missing columns: {sorted(missing)}")
    df["pos"] = df["pos"].astype(int)
    df["type"] = df["type"].str.lower()
    bad = df[~df["type"].isin({"snv", "ins", "del"})]
    if len(bad):
        sys.exit(f"Unsupported variant types: {bad['type'].unique().tolist()}")
    if "strand" not in df.columns:
        df["strand"] = ""
    df["strand"] = df["strand"].fillna("").astype(str).str.strip()
    bad_s = df[~df["strand"].isin({"", "+", "-"})]
    if len(bad_s):
        sys.exit(f"Unsupported strand values (use '+', '-', or leave blank): "
                 f"{bad_s['strand'].unique().tolist()}")
    return df

def call_snv(bam, chrom, pos1, ref, alt, min_baseq, min_mapq):
    counts = defaultdict(Counter)
    reads_seen = defaultdict(set)
    pos0 = pos1 - 1
    for col in bam.pileup(
        chrom, pos0, pos0 + 1,
        truncate=True,
        min_base_quality=min_baseq,
        min_mapping_quality=min_mapq,
        ignore_overlaps=False,
        stepper="nofilter",
    ):
        if col.reference_pos != pos0:
            continue
        for pr in col.pileups:
            read = pr.alignment
            if read.is_secondary or read.is_supplementary or read.is_unmapped:
                continue
            try:
                sm = read.get_tag("SM")
            except KeyError:
                sm = "NA"
            if read.query_name in reads_seen[sm]:
                continue
            reads_seen[sm].add(read.query_name)
            if pr.is_del:
                counts[sm]["del"] += 1
            elif pr.query_position is None:
                counts[sm]["skip"] += 1
            else:
                base = read.query_sequence[pr.query_position].upper()
                counts[sm][base] += 1
    rows = []
    for sm, c in counts.items():
        depth = sum(c.values())
        ref_n = c.get(ref.upper(), 0)
        alt_n = c.get(alt.upper(), 0)
        other = depth - ref_n - alt_n - c.get("del", 0) - c.get("skip", 0)
        vaf = (alt_n / depth) if depth else 0.0
        rows.append({
            "sample": sm,
            "depth": depth,
            "ref_count": ref_n,
            "alt_count": alt_n,
            "del_count": c.get("del", 0),
            "other_count": other,
            "vaf": vaf,
        })
    return rows

def call_insertion(bam, chrom, pos1, ref, alt, min_mapq):
    expected_ins = alt[len(ref):].upper()
    pos0 = pos1 - 1
    per_sample = defaultdict(lambda: {
        "spanning": 0,
        "ins_any": 0,
        "ins_exact": 0,
        "ins_lengths": Counter(),
        "ins_seqs": Counter(),
        "reads": set(),
    })
    for read in bam.fetch(chrom, pos0, pos0 + 1):
        if read.is_secondary or read.is_supplementary or read.is_unmapped:
            continue
        if read.mapping_quality < min_mapq:
            continue
        try:
            sm = read.get_tag("SM")
        except KeyError:
            sm = "NA"
        if read.query_name in per_sample[sm]["reads"]:
            continue
        per_sample[sm]["reads"].add(read.query_name)

        aligned = read.get_aligned_pairs(matches_only=False)
        anchor_idx = None
        for i, (qpos, rpos) in enumerate(aligned):
            if rpos == pos0 and qpos is not None:
                anchor_idx = i
                break
        if anchor_idx is None:
            continue
        per_sample[sm]["spanning"] += 1
        ins_bases = []
        for qpos, rpos in aligned[anchor_idx + 1:]:
            if rpos is None and qpos is not None:
                ins_bases.append(read.query_sequence[qpos].upper())
            else:
                break
        if ins_bases:
            ins_seq = "".join(ins_bases)
            per_sample[sm]["ins_any"] += 1
            per_sample[sm]["ins_lengths"][len(ins_seq)] += 1
            per_sample[sm]["ins_seqs"][ins_seq] += 1
            if ins_seq == expected_ins:
                per_sample[sm]["ins_exact"] += 1
    rows = []
    for sm, d in per_sample.items():
        depth = d["spanning"]
        ins_any = d["ins_any"]
        ins_exact = d["ins_exact"]
        top_seqs = ", ".join(
            f"{seq}:{n}" for seq, n in d["ins_seqs"].most_common(3)
        )
        rows.append({
            "sample": sm,
            "depth": depth,
            "ref_count": depth - ins_any,
            "alt_count": ins_exact,
            "ins_any_count": ins_any,
            "ins_other_count": ins_any - ins_exact,
            "top_inserted_seqs": top_seqs,
            "vaf": (ins_exact / depth) if depth else 0.0,
        })
    return rows

def call_deletion(bam, chrom, pos1, ref, alt, min_mapq):
    expected_del_len = len(ref) - len(alt)
    pos0 = pos1 - 1
    per_sample = defaultdict(lambda: {
        "spanning": 0,
        "del_any": 0,
        "del_exact": 0,
        "del_lengths": Counter(),
        "reads": set(),
    })
    for read in bam.fetch(chrom, pos0, pos0 + 1):
        if read.is_secondary or read.is_supplementary or read.is_unmapped:
            continue
        if read.mapping_quality < min_mapq:
            continue
        try:
            sm = read.get_tag("SM")
        except KeyError:
            sm = "NA"
        if read.query_name in per_sample[sm]["reads"]:
            continue
        per_sample[sm]["reads"].add(read.query_name)

        aligned = read.get_aligned_pairs(matches_only=False)
        anchor_idx = None
        for i, (qpos, rpos) in enumerate(aligned):
            if rpos == pos0 and qpos is not None:
                anchor_idx = i
                break
        if anchor_idx is None:
            continue
        per_sample[sm]["spanning"] += 1
        del_len = 0
        for qpos, rpos in aligned[anchor_idx + 1:]:
            if qpos is None and rpos is not None:
                del_len += 1
            else:
                break
        if del_len > 0:
            per_sample[sm]["del_any"] += 1
            per_sample[sm]["del_lengths"][del_len] += 1
            if del_len == expected_del_len:
                per_sample[sm]["del_exact"] += 1
    rows = []
    for sm, d in per_sample.items():
        depth = d["spanning"]
        del_any = d["del_any"]
        del_exact = d["del_exact"]
        top_lens = ", ".join(
            f"{ln}:{n}" for ln, n in d["del_lengths"].most_common(3)
        )
        rows.append({
            "sample": sm,
            "depth": depth,
            "ref_count": depth - del_any,
            "alt_count": del_exact,
            "del_any_count": del_any,
            "del_other_count": del_any - del_exact,
            "top_deleted_lengths": top_lens,
            "vaf": (del_exact / depth) if depth else 0.0,
        })
    return rows


def main():
    ap = argparse.ArgumentParser(
        prog="genevar",
        description="BaseCode GeneVar — inspect BAM read support for known variant sites, per sample (split by SM tag).",
    )
    ap.add_argument("-i", "--bam", required=True, help="Input BAM.")
    ap.add_argument("-v", "--variants", required=True,
                    help="Variants TSV/CSV with columns: "
                         "name chrom pos ref alt type")
    ap.add_argument("-o", "--out", required=True,
                    help="Output TSV: per-(variant, sample) call counts.")
    ap.add_argument("--min-baseq", type=int, default=0,
                    help="Minimum base quality for SNV pileup (default: 0).")
    ap.add_argument("--min-mapq", type=int, default=0,
                    help="Minimum mapping quality (default: 0).")
    ap.add_argument("--no-check-conversion", dest="check_conversion",
                    action="store_false",
                    help=("Skip the BaseCode conversion-ambiguity check. By "
                          "default, GeneVar flags variants whose REF/ALT "
                          "collides with BaseCode chemistry."))
    ap.set_defaults(check_conversion=True)
    ap.add_argument("--version", action="version",
                    version=f"BaseCode GeneVar {__version__}")
    args = ap.parse_args()

    variants = parse_variants(args.variants)
    bam = pysam.AlignmentFile(args.bam, "rb")

    out_rows = []
    for _, v in variants.iterrows():
        strand = v.get("strand", "") or ""
        ambiguous = is_conversion_ambiguous(v["type"], v["ref"], v["alt"], strand)
        if args.check_conversion and ambiguous:
            if strand == "+":
                rule = "+ strand gene, G->A is the conversion direction"
            elif strand == "-":
                rule = "- strand gene, C->T is the conversion direction"
            else:
                rule = "strand unknown, flagging both G<->A and C<->T"
            print(
                f"Warning: variant {v['name']} ({v['chrom']}:{v['pos']} "
                f"{v['ref']}>{v['alt']}) is indistinguishable from BaseCode "
                f"chemistry noise ({rule}). Alt-supporting reads may reflect "
                f"a base modification rather than a true variant.",
                file=sys.stderr,
            )
        if v["type"] == "snv":
            rows = call_snv(
                bam, v["chrom"], int(v["pos"]), v["ref"], v["alt"],
                args.min_baseq, args.min_mapq,
            )
        elif v["type"] == "ins":
            rows = call_insertion(
                bam, v["chrom"], int(v["pos"]), v["ref"], v["alt"],
                args.min_mapq,
            )
        else:
            rows = call_deletion(
                bam, v["chrom"], int(v["pos"]), v["ref"], v["alt"],
                args.min_mapq,
            )
        for r in rows:
            r2 = {"variant": v["name"],
                  "chrom": v["chrom"], "pos": int(v["pos"]),
                  "ref": v["ref"], "alt": v["alt"], "type": v["type"],
                  "strand": strand}
            r2.update(r)
            if args.check_conversion:
                r2["conversion_ambiguous"] = ambiguous
            out_rows.append(r2)

    if not out_rows:
        sys.exit("No reads found at any variant position")

    df = pd.DataFrame(out_rows)
    leading = ["variant", "chrom", "pos", "ref", "alt", "type", "strand",
               "conversion_ambiguous", "sample",
               "depth", "ref_count", "alt_count", "vaf"]
    leading = [c for c in leading if c in df.columns]
    rest = [c for c in df.columns if c not in leading]
    df = df[leading + rest].sort_values(["variant", "sample"])
    df.to_csv(args.out, sep="\t", index=False)
    print(f"Wrote {len(df)} rows to {args.out}", file=sys.stderr)


if __name__ == "__main__":
    main()
