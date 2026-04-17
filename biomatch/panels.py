import os
import sys
import glob
import shutil
from .utils import KMER_REF_PANELS_DIR, ANALYSIS_SCRIPTS_DIR, SPECIES_AUTOSOMES, _normalize, run_biomatch_proc


def ensure_species_dir(species_name: str) -> str:
    existing = [d for d in os.listdir(KMER_REF_PANELS_DIR) if os.path.isdir(os.path.join(KMER_REF_PANELS_DIR, d))]
    match = None
    for d in existing:
        if _normalize(d) == _normalize(species_name):
            match = d
            break
    if match is None:
        match = species_name.strip()
    sp_dir = os.path.join(KMER_REF_PANELS_DIR, match)
    os.makedirs(sp_dir, exist_ok=True)
    return sp_dir


def find_panel_by_name(panel_name: str) -> str | None:
    target = _normalize(panel_name.replace(".fa", ""))
    for species_dir in sorted(os.listdir(KMER_REF_PANELS_DIR)):
        sp_path = os.path.join(KMER_REF_PANELS_DIR, species_dir)
        if not os.path.isdir(sp_path):
            continue
        for fn in os.listdir(sp_path):
            if fn.endswith(".fa"):
                base = fn[:-3]
                if _normalize(base) == target:
                    return os.path.join(sp_path, fn)
    return None


def generate_panel(args) -> bool:
    if not args.ref or not args.pop_vcf or not args.panel_name:
        print("Error: --ref, --pop-vcf, and --panel-name are required for panel generation")
        return False

    base_name = args.panel_name if args.panel_name.endswith('.fa') else args.panel_name + '.fa'
    base = base_name[:-3]
    species = base.split('_ref_')[0] if '_ref_' in base else base
    species_dir = ensure_species_dir(species)

    chr_set = args.chr_set if getattr(args, 'chr_set', None) else SPECIES_AUTOSOMES.get(_normalize(species))

    work_dir = args.out if getattr(args, 'out', None) else os.path.join(species_dir, f"tmp_{base}")
    os.makedirs(work_dir, exist_ok=True)
    input_dir = os.path.join(work_dir, "input")
    os.makedirs(input_dir, exist_ok=True)

    # 1) Prepare reference
    try:
        ref_src = args.ref
        if ref_src.endswith('.gz'):
            import gzip
            unzipped_ref = os.path.join(work_dir, os.path.basename(ref_src).replace('.gz', ''))
            print(f"Unzipping reference via Python: {ref_src} -> {unzipped_ref}")
            with gzip.open(ref_src, 'rb') as fin, open(unzipped_ref, 'wb') as fout:
                shutil.copyfileobj(fin, fout)
        else:
            unzipped_ref = os.path.join(work_dir, os.path.basename(ref_src))
            shutil.copy2(ref_src, unzipped_ref)
    except Exception as e:
        print(f"Error preparing reference: {e}")
        return False

    # 2) Filter population VCF via plink2
    try:
        filter_prefix = os.path.join(input_dir, "maf45_55_bi_allelic")
        p2_cmd = (
            f"plink2 --vcf '{args.pop_vcf}' --vcf-half-call missing --snps-only --max-alleles 2 "
            + (f"--chr-set {chr_set} " if chr_set else "")
            + "--allow-extra-chr --maf 0.45 --max-maf 0.55 --make-pgen --out '" + filter_prefix + "'"
        )
        print(f"Filtering VCF via PLINK2: {p2_cmd}")
        run_biomatch_proc(p2_cmd)
    except Exception as e:
        print(f"Error filtering VCF with plink2: {e}")
        return False

    # 3) Export filtered dataset back to VCF
    try:
        export_cmd = (
            f"plink2 --pfile '{filter_prefix}' "
            + (f"--chr-set {chr_set} " if chr_set else "")
            + "--allow-extra-chr --export vcf --out '" + filter_prefix + "'"
        )
        print(f"Exporting filtered dataset to VCF: {export_cmd}")
        run_biomatch_proc(export_cmd)
        filtered_vcf = filter_prefix + ".vcf"
    except Exception as e:
        print(f"Error exporting VCF with plink2: {e}")
        return False

    # 4) Rename reference headers using ref_map_new.py with VCF-aware prefix ('auto')
    try:
        refmap_script = os.path.join(ANALYSIS_SCRIPTS_DIR, "ref_map_new.py")
        renamed_ref = os.path.join(work_dir, os.path.splitext(os.path.basename(unzipped_ref))[0] + "_renamed.fna")
        refmap_cmd = (
            f"'{sys.executable}' '{refmap_script}' -ref '{unzipped_ref}' -out '{renamed_ref}' "
            f"-mode auto -vcf '{filtered_vcf}' -prefix auto -quiet"
        )
        print(f"Renaming reference headers with VCF-aware prefix: {refmap_cmd}")
        run_biomatch_proc(refmap_cmd)
    except Exception as e:
        print(f"Error in reference renaming: {e}")
        return False

    # 5) Generate k-mer reference panel using ntsmSiteGen
    try:
        ntsm_cmd = (
            f"ntsmSiteGen generate-sites name={species} ref='{renamed_ref}' vcf='{filtered_vcf}'"
        )
        print(f"Generating k-mer panel via ntsmSiteGen in {work_dir}: {ntsm_cmd}")
        run_biomatch_proc(ntsm_cmd, cwd=work_dir)
    except Exception as e:
        print(f"Error generating k-mer panel via ntsmSiteGen: {e}")
        return False

    # 6) Install final panel into species dir
    try:
        candidates = sorted(glob.glob(os.path.join(work_dir, f"{species}_n*.fa")))
        if not candidates:
            print(f"Error: No stratified panel files found under {work_dir}")
            return False
        target_sequences = 150000

        def count_lines(path: str) -> int:
            try:
                with open(path, 'r') as fh:
                    return sum(1 for _ in fh)
            except Exception:
                return 0

        scored = []
        for p in candidates:
            lines = count_lines(p)
            seqs = lines // 2
            scored.append((p, seqs, abs(seqs - target_sequences)))
        scored.sort(key=lambda x: (x[2], -x[1]))
        chosen_path, chosen_seqs, _ = scored[0]

        final_panel_path = os.path.join(species_dir, base_name)
        shutil.copy2(chosen_path, final_panel_path)
        print(f"Panel generated and installed: {final_panel_path}")
        print(f"Selected panel: {os.path.basename(chosen_path)} with ~{chosen_seqs} sequences (wc -l/2)")
    except Exception as e:
        print(f"Error installing final panel: {e}")
        return False

    try:
        shutil.rmtree(work_dir, ignore_errors=True)
    except Exception:
        pass
    return True