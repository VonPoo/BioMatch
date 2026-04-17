#!/usr/bin/env python3
"""
BioMatch: A data-driven framework for comprehensive sample identification
Main script for handling all processing modes
"""

import os
import sys
import argparse
import subprocess
import shutil
import tempfile
from pathlib import Path
import pandas as pd
from concurrent.futures import ThreadPoolExecutor, as_completed
from tqdm import tqdm

# Package version —— 保持与 meta.yaml 和 setup.py 一致
__version__ = "1.0.0"

# Get the package directory
PACKAGE_DIR = os.path.dirname(os.path.abspath(__file__))
ANALYSIS_SCRIPTS_DIR = os.path.join(PACKAGE_DIR, "analysis_scripts")


# ---------------------------------------------------------------------------
# Resource discovery (conda-aware)
# ---------------------------------------------------------------------------

def _get_install_prefix() -> str:
    """Return the active install prefix ($CONDA_PREFIX, then sys.prefix)."""
    return os.environ.get("CONDA_PREFIX") or sys.prefix


def get_panels_root() -> str:
    """Prefer panels installed under $PREFIX/share/biomatch; fallback to package data."""
    prefix = _get_install_prefix()
    share_panels = os.path.join(prefix, "share", "biomatch", "kmer_ref_panels")
    if os.path.isdir(share_panels):
        return share_panels
    return os.path.join(PACKAGE_DIR, "kmer_ref_panels")


def get_resources_dir() -> str:
    """Locate biomatch resources (makefile + override .py scripts)."""
    prefix = _get_install_prefix()
    res = os.path.join(prefix, "share", "biomatch", "resources")
    if os.path.isdir(res):
        return res
    # Dev fallback: scripts shipped alongside the package source tree
    dev = os.path.join(PACKAGE_DIR, "..", "resources")
    return os.path.abspath(dev)


def get_ntsm_scripts_dir() -> str:
    """Locate ntsm's scripts directory (installed by the `ntsm` conda package)."""
    prefix = _get_install_prefix()
    ns = os.path.join(prefix, "bin", "ntsm-scripts")
    if not os.path.isdir(ns):
        raise FileNotFoundError(
            f"ntsm-scripts not found at {ns}. "
            "Please ensure ntsm is installed: conda install -c bioconda ntsm"
        )
    return ns


KMER_REF_PANELS_DIR = get_panels_root()


# ---------------------------------------------------------------------------
# Runtime pipeline staging (replaces the old post-link.sh logic)
# ---------------------------------------------------------------------------

def build_ntsm_stage(stage_dir: str) -> str:
    """
    Assemble a staging directory that contains:
      - All scripts from the ntsm package (symlinked or copied);
      - BioMatch's override scripts on top of them.

    Returns the absolute path of the staging directory.

    This replaces the previous post-link.sh that wrote into ntsm-scripts/.
    Nothing in the ntsm package is modified — everything happens under
    stage_dir, which is fully owned by biomatch at runtime.
    """
    stage = Path(stage_dir).resolve()
    stage.mkdir(parents=True, exist_ok=True)

    ntsm_dir = Path(get_ntsm_scripts_dir())
    res_dir = Path(get_resources_dir())

    # 1) Baseline: bring all ntsm scripts into the stage
    for src in ntsm_dir.iterdir():
        if not src.is_file():
            continue
        dst = stage / src.name
        if dst.exists() or dst.is_symlink():
            dst.unlink()
        try:
            os.symlink(src, dst)
        except OSError:
            shutil.copy2(src, dst)

    # 2) Override with biomatch's optimized versions
    override_names = ("00_extractSNPsfromVCF.py", "filterRepetiveSNP.py")
    for name in override_names:
        src = res_dir / name
        if not src.is_file():
            continue
        dst = stage / name
        if dst.exists() or dst.is_symlink():
            dst.unlink()
        shutil.copy2(src, dst)
        try:
            os.chmod(dst, 0o755)
        except OSError:
            pass

    # 3) Place biomatch makefile alongside scripts as well, for convenience
    mk = res_dir / "makefile.biomatch"
    if mk.is_file():
        shutil.copy2(mk, stage / "makefile")

    return str(stage)


def run_ntsm_make(target: str, work_dir: str, make_vars: dict | None = None) -> None:
    """
    Run `make` with biomatch's makefile against a freshly staged script dir.

    Args:
        target:    make target (e.g. "generate-sites", "generate-pca-rot-mat")
        work_dir:  directory where intermediate files are produced
        make_vars: dict of make variables (e.g. {"name": "...", "ref": "..."})
    """
    work = Path(work_dir).resolve()
    work.mkdir(parents=True, exist_ok=True)

    stage = build_ntsm_stage(str(work / ".biomatch_stage"))
    makefile = os.path.join(get_resources_dir(), "makefile.biomatch")

    if not os.path.isfile(makefile):
        raise FileNotFoundError(
            f"biomatch makefile not found at {makefile}. "
            "Please reinstall biomatch."
        )

    cmd = [
        "make", "-f", makefile,
        "-C", str(work),
        f"SCRIPT_DIR={stage}",
    ]
    for k, v in (make_vars or {}).items():
        cmd.append(f"{k}={v}")
    cmd.append(target)

    env = os.environ.copy()
    env["BIOMATCH_SCRIPT_DIR"] = stage
    subprocess.run(cmd, env=env, check=True)


# ---------------------------------------------------------------------------
# ANSI color support
# ---------------------------------------------------------------------------

COLOR_MODE = os.environ.get("BIOMATCH_COLOR", "auto").lower()

def set_color_mode_from_argv():
    global COLOR_MODE
    for i, tok in enumerate(sys.argv):
        if tok == "--color" and i + 1 < len(sys.argv):
            COLOR_MODE = sys.argv[i + 1].lower()
        elif tok.startswith("--color="):
            COLOR_MODE = tok.split("=", 1)[1].lower()

def supports_color() -> bool:
    force = COLOR_MODE
    if force == "always":
        return True
    if force == "never":
        return False
    try:
        return sys.stdout.isatty()
    except Exception:
        return False

COLOR = {
    "cyan": "\033[1;36m",
    "yellow": "\033[1;33m",
    "green": "\033[1;32m",
    "magenta": "\033[1;35m",
    "red": "\033[1;31m",
    "blue": "\033[1;34m",
    "white": "\033[1;37m",
    "bold": "\033[1m",
    "dim": "\033[2m",
    "reset": "\033[0m",
    "#8bd3dd": "\033[38;2;139;211;221m",
    "#f3d2c1": "\033[38;2;243;210;193m",
    "#f582ae": "\033[38;2;245;130;174m",
}

class ColorHelpFormatter(argparse.RawDescriptionHelpFormatter):
    def start_section(self, heading):
        if supports_color():
            colored_heading = f"{COLOR['yellow']}{heading}{COLOR['reset']}"
            super().start_section(colored_heading)
        else:
            super().start_section(heading)

    def _format_action_invocation(self, action):
        if supports_color():
            if not action.option_strings:
                metavar = self._format_args(action, action.dest.upper())
                return f"{COLOR['blue']}{metavar}{COLOR['reset']}"
            else:
                parts = [f"{COLOR['green']}{o}{COLOR['reset']}" for o in action.option_strings]
                return ', '.join(parts)
        return super()._format_action_invocation(action)

    def _format_action(self, action):
        result = super()._format_action(action)
        if supports_color():
            import re
            result = re.sub(r'\b([A-Z_]+)\b', f'{COLOR["blue"]}\\1{COLOR["reset"]}', result)
            result = re.sub(r'(default: )([^)]+)', f'\\1{COLOR["dim"]}\\2{COLOR["reset"]}', result)
            result = re.sub(r'(\.[a-z]+)', f'{COLOR["magenta"]}\\1{COLOR["reset"]}', result)
        return result

    def format_usage(self):
        usage = super().format_usage()
        if supports_color():
            usage = usage.replace("usage:", f"{COLOR['cyan']}usage:{COLOR['reset']}")
            usage = usage.replace("biomatch", f"{COLOR['bold']}biomatch{COLOR['reset']}")
            import re
            usage = re.sub(r'\[([^\]]+)\]', f'[{COLOR["dim"]}\\1{COLOR["reset"]}]', usage)
            usage = re.sub(r'\{([^}]+)\}', f'{{{COLOR["yellow"]}\\1{COLOR["reset"]}}}', usage)
        return usage


# ---------------------------------------------------------------------------
# Banner
# ---------------------------------------------------------------------------

BANNER_TEMPLATE = r"""
=============================================================================
 BioMatch — A data-driven framework for comprehensive sample identification
=============================================================================

    ██████╗  ██╗ ██████╗  ███╗   ███╗  █████╗ ████████╗ ███████╗ ██╗  ██╗
    ██╔══██╗ ██║ ██╔══██╗ ████╗ ████║ ██╔══██╗╚═ ██╔══╝ ██╔════╝ ██║  ██║
    ██████╔╝ ██║ ██║  ██║ ██╔████╔██║ ███████║   ██║    ██║      ███████║
    ██╔══██╗ ██║ ██║  ██║ ██║╚██╔╝██║ ██╔══██║   ██║    ██║      ██╔══██║
    ██████╔╝ ██║ ██████╔╝ ██║ ╚═╝ ██║ ██║  ██║   ██║    ███████╗ ██║  ██║
    ╚═════╝  ╚═╝ ╚═════╝  ╚═╝     ╚═╝ ╚═╝  ╚═╝   ╚═╝    ╚══════╝ ╚═╝  ╚═╝

 BioMatch {version}
 Project: BioMatch | Written by VonPoo <fengbobo927@163.com>
 Copyright (C) 2025 Zhejiang University
"""

def get_colored_banner():
    if not supports_color():
        return BANNER_TEMPLATE.format(version=__version__)
    return f"""
{COLOR['cyan']}=============================================================================
 BioMatch — A data-driven framework for comprehensive sample identification
============================================================================={COLOR['reset']}

{COLOR['#8bd3dd']}
    ██████╗  ██╗ ██████╗  ███╗   ███╗  █████╗ ████████╗ ███████╗ ██╗  ██╗
    ██╔══██╗ ██║ ██╔══██╗ ████╗ ████║ ██╔══██╗╚═ ██╔══╝ ██╔════╝ ██║  ██║
    ██████╔╝ ██║ ██║  ██║ ██╔████╔██║ ███████║   ██║    ██║      ███████║
    ██╔══██╗ ██║ ██║  ██║ ██║╚██╔╝██║ ██╔══██║   ██║    ██║      ██╔══██║
    ██████╔╝ ██║ ██████╔╝ ██║ ╚═╝ ██║ ██║  ██║   ██║    ███████╗ ██║  ██║
    ╚═════╝  ╚═╝ ╚═════╝  ╚═╝     ╚═╝ ╚═╝  ╚═╝   ╚═╝    ╚══════╝ ╚═╝  ╚═╝
    {COLOR['reset']}

{COLOR['#f3d2c1']} BioMatch {COLOR['bold']}{__version__}{COLOR['reset']}
{COLOR['#f582ae']} Project: BioMatch | Written by VonPoo <fengbobo927@163.com>
 Copyright (C) 2025 Zhejiang University{COLOR['reset']}
"""

def print_banner():
    print(get_colored_banner())


# ---------------------------------------------------------------------------
# Species map
# ---------------------------------------------------------------------------

SPECIES_AUTOSOMES = {
    "arabian_camel": 36, "cat": 18, "cattle": 29, "chicken": 38,
    "chimpanzee": 23, "chukar_partridge": 38, "cobitidae": 23,
    "darwin_finches": 39, "dog": 38, "domestic_yak": 29, "donkey": 30,
    "dugong": 28, "fox": 16, "giant_panda": 20, "goat": 29, "grivet": 29,
    "honey_bee": 15, "horse": 31, "human": 22, "lion": 18, "macaque": 20,
    "mallard": 39, "mouse": 19, "norway_rat": 20, "pig": 18, "rabbit": 21,
    "sheep": 26, "swan_goose": 39, "turkey": 39, "water_buffalo": 24,
    "zebra_finch": 39,
}


# ---------------------------------------------------------------------------
# Panel helpers
# ---------------------------------------------------------------------------

def _normalize(s: str) -> str:
    return s.strip().lower()

def find_panel_by_name(panel_name: str) -> str | None:
    target = _normalize(panel_name.replace(".fa", ""))
    if not os.path.isdir(KMER_REF_PANELS_DIR):
        return None
    for species_dir in sorted(os.listdir(KMER_REF_PANELS_DIR)):
        sp_path = os.path.join(KMER_REF_PANELS_DIR, species_dir)
        if not os.path.isdir(sp_path):
            continue
        for fn in os.listdir(sp_path):
            if fn.endswith(".fa") and _normalize(fn[:-3]) == target:
                return os.path.join(sp_path, fn)
    return None

def ensure_species_dir(species_name: str) -> str:
    os.makedirs(KMER_REF_PANELS_DIR, exist_ok=True)
    existing = [d for d in os.listdir(KMER_REF_PANELS_DIR)
                if os.path.isdir(os.path.join(KMER_REF_PANELS_DIR, d))]
    match = next((d for d in existing if _normalize(d) == _normalize(species_name)), None)
    if match is None:
        match = species_name.strip()
    sp_dir = os.path.join(KMER_REF_PANELS_DIR, match)
    os.makedirs(sp_dir, exist_ok=True)
    return sp_dir


# ---------------------------------------------------------------------------
# Argument parser
# ---------------------------------------------------------------------------

def setup_parser():
    desc = "BioMatch: A data-driven framework for comprehensive sample identification"
    if supports_color():
        desc = f"{COLOR['cyan']}{desc}{COLOR['reset']}"

    parser = argparse.ArgumentParser(
        description=desc,
        formatter_class=ColorHelpFormatter,
        epilog="""
Processing Modes:
----------------
1. Generate Panel Only:
   biomatch --gen-panel --ref *.fa --pop-vcf *.vcf --chr-set * --panel-name **_ref_version

2. Generate Panel, Count & Eval:
   biomatch --gen-panel --ref *.fa --pop-vcf *.vcf --chr-set * --panel-name **_ref_version --count **/Dir --count-db count_result_path --eval-result eval_result_path

3. Count & Eval using Existing Panel:
   biomatch --panel-name **_ref_version --count **/Dir --count-db count_result_path -t 20 --eval-result eval_result_path

4. BioMatch Eval on VCF/PLINK:
   biomatch --match **.vcf/** --species * [or --chr N] --eval-result eval_result_path

5. BioMatch Eval with Base-Set Keep (recommended):
   biomatch --match **.vcf/** --species * [or --chr N] --keep-base <bases> --eval-result eval_result_path

6. Default Eval on Count Results:
   biomatch --count-db count_result_path --eval-result eval_result_path
        """
    )

    panel_group = parser.add_argument_group("Panel Generation Options")
    panel_group.add_argument("--gen-panel", action="store_true", help="Generate a k-mer panel")
    panel_group.add_argument("--ref", help="Reference genome FASTA file")
    panel_group.add_argument("--pop-vcf", help="Population VCF file")
    panel_group.add_argument("--panel-name", help="Name for the panel to generate or use")
    panel_group.add_argument("--chr-set", type=int, help="Autosome count for PLINK2 chr-set (e.g., 22)")
    panel_group.add_argument("--out", help="Output directory for intermediate results")

    count_group = parser.add_argument_group("Counting Options")
    count_group.add_argument("--count", help="FASTA/FASTQ file or directory")
    count_group.add_argument("--count-db", help="Path to save counting results")
    count_group.add_argument("-t", "--threads", type=int, default=10, help="Number of parallel jobs")

    eval_group = parser.add_argument_group("Evaluation Options")
    eval_group.add_argument("--match", help="VCF file or PLINK prefix for BioMatch evaluation")
    eval_group.add_argument("--species", help="Species name for autosome count")
    eval_group.add_argument("--chromosomes", "--chr", dest="chromosomes", type=int,
                            help="Override autosome chromosome count")
    eval_group.add_argument("--keep-base",
        help="Allowed bases set for retention (e.g., ATC or A,T,C).")
    eval_group.add_argument("--eval-result", help="Path to save evaluation results")

    parser.add_argument("--color", choices=["auto", "always", "never"],
                        default=os.environ.get("BIOMATCH_COLOR", "auto"),
                        help="Colorized help output mode")
    parser.add_argument("--list-panels", action="store_true",
                        help="List built-in panel names")
    parser.add_argument("--version", action="version", version=f"BioMatch {__version__}")

    return parser


# ---------------------------------------------------------------------------
# Lightweight process helper
# ---------------------------------------------------------------------------

def run_biomatch_proc(command: str, check: bool = True,
                      env: dict | None = None, cwd: str | None = None):
    """Run a shell command using bash -lc."""
    return subprocess.run(["bash", "-lc", command], check=check, env=env, cwd=cwd)


# ---------------------------------------------------------------------------
# Mode delegates
# ---------------------------------------------------------------------------

def generate_panel(args):
    from .panels import generate_panel as _generate_panel
    return _generate_panel(args)

def run_counting(args):
    from .counting import run_counting as _run_counting
    return _run_counting(args)

def run_ntsm_eval(args):
    from .evaluation import run_ntsm_eval as _run_ntsm_eval
    return _run_ntsm_eval(args)

def run_deepkin_eval(args):
    from .evaluation import run_deepkin_eval as _run_deepkin_eval
    return _run_deepkin_eval(args)

def run_default_eval(args):
    from .evaluation import run_ntsm_eval as _run_ntsm_eval
    return _run_ntsm_eval(args)


# ---------------------------------------------------------------------------
# Main
# ---------------------------------------------------------------------------

def main():
    # 启动 banner
    print_banner()

    set_color_mode_from_argv()
    parser = setup_parser()
    args = parser.parse_args()

    os.environ["BIOMATCH_COLOR"] = args.color
    globals()["COLOR_MODE"] = args.color

    if getattr(args, "list_panels", False):
        print("Available panel names:")
        try:
            names = []
            if os.path.isdir(KMER_REF_PANELS_DIR):
                for species_dir in sorted(os.listdir(KMER_REF_PANELS_DIR)):
                    sp_path = os.path.join(KMER_REF_PANELS_DIR, species_dir)
                    if not os.path.isdir(sp_path):
                        continue
                    for fn in sorted(os.listdir(sp_path)):
                        if fn.endswith('.fa'):
                            names.append(os.path.splitext(fn)[0])
            if names:
                for nm in sorted(names):
                    print(f" - {nm}")
            else:
                print("(no panels found)")
        except Exception as e:
            print(f"Error listing panels: {e}")
        return 0

    if args.gen_panel:
        if not generate_panel(args):
            return 1
        if args.count:
            if not run_counting(args):
                return 1
            if args.eval_result and not run_ntsm_eval(args):
                return 1

    elif args.match:
        if not run_deepkin_eval(args):
            return 1

    elif args.panel_name and args.count:
        if not run_counting(args):
            return 1
        if args.eval_result and not run_ntsm_eval(args):
            return 1

    elif args.count_db and args.eval_result:
        if not run_ntsm_eval(args):
            return 1

    else:
        parser.print_help()
        return 1

    return 0


if __name__ == "__main__":
    sys.exit(main())
