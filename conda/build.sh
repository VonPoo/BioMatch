#!/bin/bash
set -euo pipefail

# 安装BioMatch
$PYTHON -m pip install --no-deps --ignore-installed .

# 将资源脚本安装到 share 目录，避免与 ntsm 产生路径冲突
RES_DIR="$PREFIX/share/biomatch/resources"
mkdir -p "$RES_DIR"
# Install our updated scripts; keep extractor under 00_ name and post-link will alias
cp "$SRC_DIR/00_extractSNPsfromVCF.py" "$RES_DIR/00_extractSNPsfromVCF.py"
# Use the SAM-based filter script from package analysis_scripts
cp "$SRC_DIR/biomatch/analysis_scripts/filterRepetiveSNP.py" "$RES_DIR/filterRepetiveSNP.py"
cp "$RECIPE_DIR/resources/makefile" "$RES_DIR/makefile.biomatch"
chmod +x "$RES_DIR/00_extractSNPsfromVCF.py" "$RES_DIR/filterRepetiveSNP.py"

# 安装参考面板到 share 目录，供运行时使用
PANELS_DST="$PREFIX/share/biomatch/kmer_ref_panels"
mkdir -p "$PANELS_DST"
if [ -d "$SRC_DIR/biomatch/kmer_ref_panels" ]; then
  cp -r "$SRC_DIR/biomatch/kmer_ref_panels"/* "$PANELS_DST/" || true
fi