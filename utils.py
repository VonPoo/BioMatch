import os
import sys
import shutil
import subprocess

# Version comes from package __init__
try:
    from . import __version__
except Exception:
    __version__ = "0.0.0"

# Package directories
PACKAGE_DIR = os.path.dirname(os.path.abspath(__file__))
ANALYSIS_SCRIPTS_DIR = os.path.join(PACKAGE_DIR, "analysis_scripts")


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


class ColorHelpFormatter:
    # Local lightweight wrapper to avoid argparse dependency here;
    # biomatch.py uses argparse.RawDescriptionHelpFormatter with color.
    pass


def get_colored_banner() -> str:
    if not supports_color():
        return BANNER_TEMPLATE.format(version=__version__)
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


def print_banner():
    print(get_colored_banner())


def run_biomatch_proc(command: str, check: bool = True, env: dict | None = None, cwd: str | None = None):
    cmd = ["bash", "-lc", command]
    return subprocess.run(cmd, check=check, env=env, cwd=cwd)


def safe_copy(src_path: str, dst_path: str):
    try:
        if os.path.abspath(src_path) == os.path.abspath(dst_path):
            print(f"Notice: source equals destination, skipping copy for {dst_path}")
            return
        shutil.copy2(src_path, dst_path)
    except shutil.SameFileError:
        print(f"Notice: SameFileError on copy {src_path} -> {dst_path}, skipping")


def copy_scripts():
    source_dir = "/disk227/fengbo/biomatch_software"
    scripts = [
        "00_extractSNPsfromVCF.py",
        "filterRepetiveSNP.py",
        "PLINK_Geno_nonGeno.R",
        "VCF_Geno_nonGeno.R",
        "ref_map_new.py",
        "special_sequence_method_snp_filter.py",
    ]
    for script in scripts:
        source = os.path.join(source_dir, script)
        dest = os.path.join(ANALYSIS_SCRIPTS_DIR, script)
        if not os.path.exists(dest) and os.path.exists(source):
            shutil.copy2(source, dest)
            if script.endswith(".py"):
                os.chmod(dest, 0o755)


def _normalize(s: str) -> str:
    return s.strip().lower()


# Species autosome mapping
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
}