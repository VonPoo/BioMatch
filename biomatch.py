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
import glob
from pathlib import Path
import pandas as pd
from concurrent.futures import ThreadPoolExecutor, as_completed
from tqdm import tqdm

# Package version
__version__ = "0.4.4"

# Get the package directory
PACKAGE_DIR = os.path.dirname(os.path.abspath(__file__))
ANALYSIS_SCRIPTS_DIR = os.path.join(PACKAGE_DIR, "analysis_scripts")

# Prefer panels installed under conda share if available; fallback to package data
def get_panels_root():
    conda_prefix = os.environ.get("CONDA_PREFIX")
    if conda_prefix:
        share_panels = os.path.join(conda_prefix, "share", "biomatch", "kmer_ref_panels")
        if os.path.isdir(share_panels):
            return share_panels
    return os.path.join(PACKAGE_DIR, "kmer_ref_panels")

KMER_REF_PANELS_DIR = get_panels_root()

# ANSI color support
COLOR_MODE = os.environ.get("BIOMATCH_COLOR", "auto").lower()

def set_color_mode_from_argv():
    global COLOR_MODE
    # Pre-scan argv to respect --color before parser renders help
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
}

# Add truecolor hex mappings used by banner to avoid KeyError
COLOR.update({
    "#8bd3dd": "\033[38;2;139;211;221m",
    "#f3d2c1": "\033[38;2;243;210;193m",
    "#f582ae": "\033[38;2;245;130;174m",
})

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
                # Positional arguments
                metavar = self._format_args(action, action.dest.upper())
                return f"{COLOR['blue']}{metavar}{COLOR['reset']}"
            else:
                # Optional arguments
                parts = []
                for option_string in action.option_strings:
                    parts.append(f"{COLOR['green']}{option_string}{COLOR['reset']}")
                return ', '.join(parts)
        else:
            return super()._format_action_invocation(action)

    def _format_action(self, action):
        # Get the basic formatted action
        result = super()._format_action(action)
        
        if supports_color():
            # Color metavars (like COUNT, PANEL_NAME, etc.)
            import re
            result = re.sub(r'\b([A-Z_]+)\b', f'{COLOR["blue"]}\\1{COLOR["reset"]}', result)
            
            # Color default values
            result = re.sub(r'(default: )([^)]+)', f'\\1{COLOR["dim"]}\\2{COLOR["reset"]}', result)
            
            # Color file paths and extensions
            result = re.sub(r'(\.[a-z]+)', f'{COLOR["magenta"]}\\1{COLOR["reset"]}', result)
        
        return result

    def format_usage(self):
        usage = super().format_usage()
        if supports_color():
            # Color the usage line
            usage = usage.replace("usage:", f"{COLOR['cyan']}usage:{COLOR['reset']}")
            usage = usage.replace("biomatch", f"{COLOR['bold']}biomatch{COLOR['reset']}")
            
            # Color optional and required parts
            import re
            usage = re.sub(r'\[([^\]]+)\]', f'[{COLOR["dim"]}\\1{COLOR["reset"]}]', usage)
            usage = re.sub(r'\{([^}]+)\}', f'{{{COLOR["yellow"]}\\1{COLOR["reset"]}}}', usage)
        return usage

# Startup banner with color support
def get_colored_banner():
    """Generate colored banner based on color mode"""
    if not supports_color():
        return BANNER_TEMPLATE.format(version=__version__)
    
    # Colored version of the banner
    colored_banner = f"""
{COLOR['cyan']}=============================================================================
 BioMatch — A data-driven framework for comprehensive sample identification
============================================================================={COLOR['reset']}

{COLOR.get('#8bd3dd', COLOR['cyan'])}
    ██████╗  ██╗ ██████╗  ███╗   ███╗  █████╗ ████████╗ ███████╗ ██╗  ██╗
    ██╔══██╗ ██║ ██╔══██╗ ████╗ ████║ ██╔══██╗╚═ ██╔══╝ ██╔════╝ ██║  ██║
    ██████╔╝ ██║ ██║  ██║ ██╔████╔██║ ███████║   ██║    ██║      ███████║
    ██╔══██╗ ██║ ██║  ██║ ██║╚██╔╝██║ ██╔══██║   ██║    ██║      ██╔══██║
    ██████╔╝ ██║ ██████╔╝ ██║ ╚═╝ ██║ ██║  ██║   ██║    ███████╗ ██║  ██║
    ╚═════╝  ╚═╝ ╚═════╝  ╚═╝     ╚═╝ ╚═╝  ╚═╝   ╚═╝    ╚══════╝ ╚═╝  ╚═╝
    {COLOR['reset']}

{COLOR.get('#f3d2c1', COLOR['yellow'])} BioMatch {COLOR['bold']}{__version__}{COLOR['reset']}
{COLOR.get('#f582ae', COLOR['magenta'])} Project: BioMatch | Written by VonPoo <fengbobo927@163.com>
 Copyright (C) 2025 Zhejiang University{COLOR['reset']}
"""
    return colored_banner

# Plain banner template for non-color mode
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
 Project: Biomatch | Written by VonPoo <fengbobo927@163.com>
 Copyright (C) 2025 Zhejiang University
"""

# Dynamic banner property
BANNER = get_colored_banner()

# Species to autosome count mapping (extendable)
SPECIES_AUTOSOMES = {
    "arabian_camel": 36,
    "cat": 18,
    "cattle": 29,
    "chicken": 38,
    "chimpanzee": 23,
    "chukar_partridge": 38,
    "cobitidae": 23,
    "darwin_finches": 39,
    "dog": 38,
    "domestic_yak": 29,
    "donkey": 30,
    "dugong": 28,
    "fox": 16,
    "giant_panda": 20,
    "goat": 29,
    "grivet": 29,
    "honey_bee": 15,
    "horse": 31,
    "human": 22,
    "lion": 18,
    "macaque": 20,
    "mallard": 39,
    "mouse": 19,
    "norway_rat": 20,
    "pig": 18,
    "rabbit": 21,
    "sheep": 26,
    "swan_goose": 39,
    "turkey": 39,
    "water_buffalo": 24,
    "zebra_finch": 39,
    # add more as needed
}

def print_banner():
    """Print the banner with color support"""
    print(get_colored_banner())

def _normalize(s: str) -> str:
    return s.strip().lower()

def find_panel_by_name(panel_name: str) -> str:
    """Case-insensitive search for a panel across species dirs.
    Accepts base name without .fa. Returns absolute path to panel .fa or None.
    """
    target = _normalize(panel_name.replace(".fa", ""))
    # Search all species subdirs
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

def ensure_species_dir(species_name: str) -> str:
    """Ensure species directory exists under panels root and return its path."""
    # Use original casing if exists; else create with Title case first letter
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

def setup_parser():
    """Set up command line argument parser"""
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

Keep-Base (Base-Set) Mode:
--------------------------
- Purpose: Recommended for special sequencing methods where base changes may occur, like WGBS.
           Specify an allowed base set and retain variants whose alleles only
           contain allowed bases, approximating "keep unchanged bases".
- Examples: --keep-base ATC (keeps A/T/C, excludes any G-containing alleles)
            --keep-base A,T (comma-separated equivalent)

6. Default Eval on Count Results:
   biomatch --count-db count_result_path --eval-result eval_result_path
        """
    )
    
    # Panel generation options
    panel_title = "Panel Generation Options"
    if supports_color():
        panel_title = f"{COLOR['yellow']}{panel_title}{COLOR['reset']}"
    panel_group = parser.add_argument_group(panel_title)
    panel_group.add_argument("--gen-panel", action="store_true", help="Generate a k-mer panel")
    panel_group.add_argument("--ref", help="Reference genome FASTA file")
    panel_group.add_argument("--pop-vcf", help="Population VCF file")
    panel_group.add_argument("--panel-name", help="Name for the panel to generate or use")
    panel_group.add_argument("--chr-set", type=int, help="Autosome count for PLINK2 chr-set (e.g., 22)")
    panel_group.add_argument("--out", help="Output directory for intermediate results (defaults to species panel dir)")
    
    # Counting options
    count_title = "Counting Options"
    if supports_color():
        count_title = f"{COLOR['green']}{count_title}{COLOR['reset']}"
    count_group = parser.add_argument_group(count_title)
    count_group.add_argument("--count", help="FASTA/FASTQ file or directory containing files for identification")
    count_group.add_argument("--count-db", help="Path to save counting results")
    count_group.add_argument("-t", "--threads", type=int, default=10, help="Number of parallel jobs for counting")
    
    # Evaluation options
    eval_title = "Evaluation Options"
    if supports_color():
        eval_title = f"{COLOR['magenta']}{eval_title}{COLOR['reset']}"
    eval_group = parser.add_argument_group(eval_title)
    eval_group.add_argument("--match", help="VCF file or PLINK format file prefix for BioMatch evaluation")
    eval_group.add_argument("--species", help="Species name for autosome count (e.g., human, cattle)")
    eval_group.add_argument("--chromosomes", "--chr", dest="chromosomes", type=int, help="Override autosome chromosome count directly (e.g., 22)")
    # Keep-base: base-set semantics; retain only variants whose alleles are composed of allowed bases
    # Accepts compact form (e.g., ATC) or comma-separated (e.g., A,T,C). Case-insensitive.
    eval_group.add_argument(
        "--keep-base",
        help=(
            "Allowed bases set for retention (e.g., ATC or A,T,C). "
            "Removes any variant where either allele contains bases outside the allowed set."
        )
    )
    eval_group.add_argument("--eval-result", help="Path to save evaluation results")
    
    # General options
    parser.add_argument("--color", choices=["auto", "always", "never"], default=os.environ.get("BIOMATCH_COLOR", "auto"), help="Colorized help output mode")
    parser.add_argument("--list-panels", action="store_true", help="List built-in panel names available to use")
    parser.add_argument("--version", action="version", version=f"BioMatch {__version__}")
    
    return parser

def copy_scripts():
    """Copy analysis scripts to the package directory"""
    source_dir = "/disk227/fengbo/biomatch_software"
    scripts = [
        "00_extractSNPsfromVCF.py",
        "filterRepetiveSNP.py",
        "PLINK_Geno_nonGeno.R",
        "VCF_Geno_nonGeno.R",
        "ref_map_new.py"
    ]
    
    for script in scripts:
        source = os.path.join(source_dir, script)
        dest = os.path.join(ANALYSIS_SCRIPTS_DIR, script)
        if not os.path.exists(dest) and os.path.exists(source):
            shutil.copy2(source, dest)
            if script.endswith(".py"):
                os.chmod(dest, 0o755)  # Make Python scripts executable

def run_biomatch_proc(command: str, check: bool = True, env: dict | None = None, cwd: str | None = None):
    """Run a shell command using bash -lc without process-name masking.

    This avoids altering argv0 (no exec -a biomatch), per user's preference.
    """
    cmd = ["bash", "-lc", command]
    return subprocess.run(cmd, check=check, env=env, cwd=cwd)

def generate_panel(args):
    """Delegate panel generation to biomatch.panels module."""
    from .panels import generate_panel as _generate_panel
    return _generate_panel(args)

def run_counting(args):
    """Delegate counting to biomatch.counting module."""
    from .counting import run_counting as _run_counting
    return _run_counting(args)

def run_ntsm_eval(args):
    """Delegate ntsm evaluation to biomatch.evaluation module."""
    from .evaluation import run_ntsm_eval as _run_ntsm_eval
    return _run_ntsm_eval(args)

def run_deepkin_eval(args):
    """Delegate BioMatch evaluation to biomatch.evaluation module."""
    from .evaluation import run_deepkin_eval as _run_deepkin_eval
    return _run_deepkin_eval(args)

def run_default_eval(args):
    """Alias default eval to ntsm evaluation in biomatch.evaluation."""
    from .evaluation import run_ntsm_eval as _run_ntsm_eval
    return _run_ntsm_eval(args)

def main():
    """Main function to handle command line arguments and run the appropriate mode"""
    # Ensure analysis scripts are available
    copy_scripts()
    
    # Print startup banner
    print_banner()
    
    # Parse command line arguments
    set_color_mode_from_argv()
    parser = setup_parser()
    args = parser.parse_args()
    # Update color mode from parsed args for subsequent outputs
    os.environ["BIOMATCH_COLOR"] = args.color
    globals()["COLOR_MODE"] = args.color
    # Optional: list all available panels and exit
    if getattr(args, "list_panels", False):
        print("Available panel names:")
        try:
            names = []
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
    
    # Determine which mode to run
    if args.gen_panel:
        # Mode 1 or 2: Generate Panel
        success = generate_panel(args)
        if not success:
            return 1
        
        # If count is specified, also run counting and ntsmEval (Mode 2)
        if args.count:
            success = run_counting(args)
            if not success:
                return 1
            
            if args.eval_result:
                success = run_ntsm_eval(args)
                if not success:
                    return 1
    
    elif args.match:
        # Mode 4: BioMatch Eval on VCF/PLINK
        success = run_deepkin_eval(args)
        if not success:
            return 1
    
    elif args.panel_name and args.count:
        # Mode 3: Count & Eval using Existing Panel
        success = run_counting(args)
        if not success:
            return 1
        
        if args.eval_result:
            success = run_ntsm_eval(args)
            if not success:
                return 1
    
    elif args.count_db and args.eval_result:
        # Mode 5: ntsmEval on Count Results
        success = run_ntsm_eval(args)
        if not success:
            return 1
    
    else:
        parser.print_help()
        return 1
    
    return 0

if __name__ == "__main__":
    sys.exit(main())
