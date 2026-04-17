import os
from concurrent.futures import ThreadPoolExecutor, as_completed
from tqdm import tqdm
from .panels import find_panel_by_name
from .utils import run_biomatch_proc


def run_counting(args) -> bool:
    if not args.count or not args.count_db or not args.panel_name:
        print("Error: --count, --count-db, and --panel-name are required for counting")
        return False

    os.makedirs(args.count_db, exist_ok=True)

    panel_path = find_panel_by_name(args.panel_name)
    if panel_path is None or not os.path.exists(panel_path):
        print(f"Error: Panel not found for name: {args.panel_name}")
        return False

    exts = (".fa", ".fasta", ".fq", ".fastq", ".fa.gz", ".fasta.gz", ".fq.gz", ".fastq.gz")
    files = []
    if os.path.isdir(args.count):
        for fn in os.listdir(args.count):
            if fn.lower().endswith(exts):
                files.append(os.path.join(args.count, fn))
    else:
        files = [args.count] if str(args.count).lower().endswith(exts) else []

    if not files:
        print(f"Error: No FASTA/FASTQ files found in {args.count}")
        return False

    def sample_key(path: str) -> str:
        base = os.path.basename(path)
        for ext in (".fa.gz", ".fasta.gz", ".fq.gz", ".fastq.gz", ".fa", ".fasta", ".fq", ".fastq"):
            if base.endswith(ext):
                base = base[: -len(ext)]
                break
        b = base
        toks = ["_R1", "_R2", ".R1", ".R2", "-R1", "-R2", "_1", "_2", "-1", "-2"]
        for t in toks:
            if b.endswith(t):
                return b[: -len(t)]
        return b

    groups: dict[str, dict[str, str]] = {}
    for path in files:
        key = sample_key(path)
        entry = groups.setdefault(key, {})
        nm = os.path.basename(path)
        low = nm.lower()
        if any(tok in low for tok in ("_r1", ".r1", "-r1")) or (low.endswith("_1") or low.endswith("-1")):
            entry["R1"] = path
        elif any(tok in low for tok in ("_r2", ".r2", "-r2")) or (low.endswith("_2") or low.endswith("-2")):
            entry["R2"] = path
        else:
            entry.setdefault("SE", path)

    jobs: list[str] = []
    for key, entry in groups.items():
        out_path = os.path.join(args.count_db, f"{key}.counts.txt")
        if "R1" in entry and "R2" in entry:
            cmd = (
                f"ntsmCount -t 1 -s '{panel_path}' '{entry['R1']}' '{entry['R2']}' > '{out_path}' 2>/dev/null"
            )
        else:
            inp = entry.get("SE") or entry.get("R1") or entry.get("R2")
            cmd = (
                f"ntsmCount -t 1 -s '{panel_path}' '{inp}' > '{out_path}' 2>/dev/null"
            )
        jobs.append(cmd)

    print(f"Running counting on {len(jobs)} sample(s) with {args.threads} parallel jobs")
    errors: list[str] = []
    try:
        with ThreadPoolExecutor(max_workers=args.threads) as ex:
            futures = [ex.submit(run_biomatch_proc, cmd) for cmd in jobs]
            pbar = tqdm(total=len(futures), desc="Counting", unit="sample")
            for fut in as_completed(futures):
                try:
                    fut.result()
                except Exception as e:
                    errors.append(str(e))
                pbar.update(1)
            pbar.close()
        if errors:
            print(f"Counting finished with {len(errors)} error(s)")
            for msg in errors:
                print(f"  - {msg}")
            return False
        print(f"Counting completed. Results saved to {args.count_db}")
        return True
    except Exception as e:
        print(f"Error running counting: {e}")
        return False