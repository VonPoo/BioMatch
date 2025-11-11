import os
import shutil
import subprocess
import pandas as pd
from .utils import ANALYSIS_SCRIPTS_DIR, SPECIES_AUTOSOMES, _normalize, run_biomatch_proc, safe_copy


def run_ntsm_eval(args) -> bool:
    if not args.count_db or not args.eval_result:
        print("Error: --count-db and --eval-result are required for ntsmEval")
        return False
    if os.path.abspath(args.count_db) == os.path.abspath(args.eval_result):
        print("Error: --count-db must be different from --eval-result")
        return False
    os.makedirs(args.eval_result, exist_ok=True)
    base = args.panel_name if getattr(args, "panel_name", None) else "eval"
    base = base[:-3] if base.endswith(".fa") else base
    out_csv = os.path.join(args.eval_result, f"{base}_BioMatch.csv")
    try:
        cmd = f"ntsmEval -a '{args.count_db}'/* > '{out_csv}' 2>/dev/null"
        print(f"Running ntsmEval: {cmd}")
        run_biomatch_proc(cmd)
        print(f"ntsmEval completed. Log saved to {out_csv}")
    except subprocess.CalledProcessError as e:
        print(f"Error running ntsmEval: {e}")
        return False
    try:
        df = pd.read_csv(out_csv, sep='\t')
        target_col = None
        if "same" in df.columns:
            target_col = "same"
        else:
            target_col = next((c for c in df.columns if c.lower() == "same"), None)
        if target_col is not None:
            filtered = df[df[target_col].astype(str) == "1"]
            filtered_match_path = os.path.join(args.eval_result, f"{base}_filtered_match.tsv")
            filtered.to_csv(filtered_match_path, sep='\t', index=False)
            after_path = out_csv.replace(".csv", ".after_filtered.csv")
            filtered.to_csv(after_path, index=False)
        else:
            print("Warning: 'same' column not found; skipped filtered output")
    except Exception as e:
        print(f"Warning: Failed to create filtered_match file: {e}")
    return True


def run_deepkin_eval(args) -> bool:
    if not args.match or not args.eval_result:
        print("Error: --match and --eval-result are required for deepkin evaluation")
        return False
    os.makedirs(args.eval_result, exist_ok=True)
    match_path = args.match
    is_dir = os.path.isdir(match_path)
    is_vcf = False
    is_plink2 = False
    if not is_dir:
        is_vcf = match_path.endswith(".vcf") or match_path.endswith(".vcf.gz")
    else:
        def any_with_ext(root, exts):
            for dp, dn, fnames in os.walk(root):
                for fn in fnames:
                    if any(fn.endswith(ext) for ext in exts):
                        return True
            return False
        has_vcf = any_with_ext(match_path, [".vcf", ".vcf.gz"])
        has_pgen = any_with_ext(match_path, [".pgen"])
        has_bed = any_with_ext(match_path, [".bed"])
        if has_vcf:
            is_vcf = True
        elif has_pgen:
            is_plink2 = True
        elif has_bed:
            is_vcf = False
        else:
            print("Error: Input directory does not contain supported VCF/PLINK files")
            return False
    file_type = "vcf" if is_vcf else "plink"

    if is_vcf:
        genomic_vcf = os.path.join(args.eval_result, "genomic.vcf.gz")
        non_genomic_vcf = os.path.join(args.eval_result, "non_genomic.vcf.gz")
        if is_dir:
            src_g = None
            src_ng = None
            for candidate in ("genomic.vcf.gz", "genomic.vcf"):
                p = os.path.join(match_path, candidate)
                if os.path.exists(p):
                    src_g = p
                    break
            for candidate in ("non_genomic.vcf.gz", "non_genomic.vcf"):
                p = os.path.join(match_path, candidate)
                if os.path.exists(p):
                    src_ng = p
                    break
            if src_g is None and src_ng is None:
                vcf_candidates = [os.path.join(match_path, fn) for fn in os.listdir(match_path) if fn.endswith(".vcf") or fn.endswith(".vcf.gz")]
                if not vcf_candidates:
                    print("Error: No VCF files found in the input directory")
                    return False
                if len(vcf_candidates) >= 2:
                    def classify(path: str) -> str:
                        nm = os.path.basename(path).lower()
                        if any(k in nm for k in ("non", "rna", "non_genomic")):
                            return "non"
                        return "gen"
                    gen_list = [p for p in vcf_candidates if classify(p) == "gen"]
                    non_list = [p for p in vcf_candidates if classify(p) == "non"]
                    if gen_list and non_list:
                        src_g = gen_list[0]
                        src_ng = non_list[0]
                    else:
                        src_g = vcf_candidates[0]
                        src_ng = vcf_candidates[1]
                else:
                    src_g = vcf_candidates[0]
            def ensure_gz(src, dst):
                if src.endswith(".vcf.gz"):
                    safe_copy(src, dst)
                elif src.endswith(".vcf"):
                    tmp_gz = os.path.join(args.eval_result, "tmp.vcf.gz")
                    run_biomatch_proc(f"bgzip -c '{src}' > '{tmp_gz}'")
                    shutil.move(tmp_gz, dst)
                else:
                    print(f"Error: Unsupported VCF file: {src}")
                    raise RuntimeError("Unsupported VCF")
            ensure_gz(src_g, genomic_vcf)
            if src_ng is not None:
                ensure_gz(src_ng, non_genomic_vcf)
        else:
            src = args.match
            if src.endswith(".vcf.gz"):
                safe_copy(src, genomic_vcf)
            elif src.endswith(".vcf"):
                tmp_gz = os.path.join(args.eval_result, "tmp.vcf.gz")
                run_biomatch_proc(f"bgzip -c '{src}' > '{tmp_gz}'")
                shutil.move(tmp_gz, genomic_vcf)
            else:
                print("Error: Unsupported VCF input. Provide .vcf or .vcf.gz")
                return False
    else:
        match_dir = match_path if is_dir else os.path.dirname(match_path)
        os.makedirs(args.eval_result, exist_ok=True)

        def find_plink2_sets(root):
            bases = []
            for dp, dn, fnames in os.walk(root):
                for fn in fnames:
                    if fn.endswith('.pgen'):
                        base = os.path.splitext(fn)[0]
                        full = os.path.join(dp, base)
                        if os.path.exists(full + '.pvar') and os.path.exists(full + '.psam'):
                            bases.append(full)
            return bases

        def find_plink1_bases(root):
            bases = []
            for dp, dn, fnames in os.walk(root):
                for fn in fnames:
                    if fn.endswith('.bed'):
                        base = os.path.splitext(fn)[0]
                        full = os.path.join(dp, base)
                        if os.path.exists(full + '.bim') and os.path.exists(full + '.fam'):
                            bases.append(full)
            return bases

        def convert_plink2(prefix_base: str, target_prefix: str):
            cmd = f"plink2 --pfile '{prefix_base}' --max-alleles 2 --make-bed --out '{os.path.join(args.eval_result, target_prefix)}'"
            run_biomatch_proc(cmd)

        converted_any = False
        if is_plink2:
            p2_bases = find_plink2_sets(match_dir)
            g_p2 = next((b for b in p2_bases if os.path.basename(b).lower() == 'genomic'), None)
            ng_p2 = next((b for b in p2_bases if os.path.basename(b).lower() == 'non_genomic'), None)
            if g_p2:
                convert_plink2(g_p2, 'genomic')
                converted_any = True
            if ng_p2:
                convert_plink2(ng_p2, 'non_genomic')
                converted_any = True
            if not converted_any and p2_bases:
                mapped = {"genomic": None, "non_genomic": None}
                for b in p2_bases:
                    name = os.path.basename(b).lower()
                    parent = os.path.basename(os.path.dirname(b)).lower()
                    if ("non" in name or "rna" in name or "non" in parent or "rna" in parent) and mapped["non_genomic"] is None:
                        mapped["non_genomic"] = b
                    elif mapped["genomic"] is None:
                        mapped["genomic"] = b
                if mapped["genomic"]:
                    convert_plink2(mapped["genomic"], 'genomic')
                    converted_any = True
                if mapped["non_genomic"]:
                    convert_plink2(mapped["non_genomic"], 'non_genomic')
                    converted_any = True
            if not converted_any:
                print("Error: No PLINK2 .pgen sets found")
                return False

        g_triplet_ok = all(os.path.exists(os.path.join(args.eval_result, f"genomic{ext}")) for ext in (".bed", ".bim", ".fam"))
        ng_triplet_ok = all(os.path.exists(os.path.join(args.eval_result, f"non_genomic{ext}")) for ext in (".bed", ".bim", ".fam"))
        if not g_triplet_ok and not ng_triplet_ok:
            bases = find_plink1_bases(match_dir)
            if not bases:
                print("Error: Missing PLINK files (.bed/.bim/.fam) in input")
                return False
            mapped = {"genomic": None, "non_genomic": None}
            for b in bases:
                name = os.path.basename(b).lower()
                parent = os.path.basename(os.path.dirname(b)).lower()
                if ("non" in name or "rna" in name or "non" in parent or "rna" in parent) and mapped["non_genomic"] is None:
                    mapped["non_genomic"] = b
                elif mapped["genomic"] is None:
                    mapped["genomic"] = b
            if mapped["genomic"]:
                for ext in (".bed", ".bim", ".fam"):
                    safe_copy(mapped["genomic"] + ext, os.path.join(args.eval_result, f"genomic{ext}"))
            if mapped["non_genomic"]:
                for ext in (".bed", ".bim", ".fam"):
                    safe_copy(mapped["non_genomic"] + ext, os.path.join(args.eval_result, f"non_genomic{ext}"))

    # Determine chromosome autosome count: prefer explicit override, otherwise map by species, fallback to 22
    chr_set = None
    if getattr(args, 'chromosomes', None):
        try:
            chr_set = int(args.chromosomes)
        except Exception:
            chr_set = None
    if chr_set is None and getattr(args, 'species', None):
        chr_set = SPECIES_AUTOSOMES.get(_normalize(args.species))
    if chr_set is None:
        chr_set = 22

    keep_base_val = (getattr(args, 'keep_base', None) or "")
    filter_script_path = os.path.join(ANALYSIS_SCRIPTS_DIR, "special_sequence_method_snp_filter.py")
    params_expr = (
        f"res <- list(chr_set={chr_set}, file_type='{file_type}', dir='{args.eval_result}', "
        f"keep_base='{keep_base_val}', script_filter_py='{filter_script_path}'); "
        f"save(res, file=file.path('{args.eval_result}','params.RData'))"
    )

    try:
        run_biomatch_proc(f"Rscript -e \"{params_expr}\"")
    except subprocess.CalledProcessError as e:
        print(f"Error creating params.RData: {e}")
        return False

    r_script = os.path.join(ANALYSIS_SCRIPTS_DIR, "VCF_Geno_nonGeno.R" if is_vcf else "PLINK_Geno_nonGeno.R")
    print(f"Running BioMatch evaluation: Rscript {r_script} {args.eval_result}")
    try:
        env = os.environ.copy()
        run_biomatch_proc(f"Rscript '{r_script}' '{args.eval_result}'", env=env)
        print(f"BioMatch evaluation completed. Results saved to {args.eval_result}")
        return True
    except subprocess.CalledProcessError as e:
        print(f"Error running BioMatch evaluation: {e}")
        return False