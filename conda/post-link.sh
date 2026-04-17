#!/bin/bash
set -euo pipefail

SCRIPTS_DIR="${PREFIX}/bin/ntsm-scripts"
RES_DIR="${PREFIX}/share/biomatch/resources"
if [ -d "${SCRIPTS_DIR}" ]; then
  echo "[post-link] Copying biomatch resources into ntsm-scripts and patching makefile..."
  # Copy our optimized scripts into ntsm-scripts without clobbering existing ones
  cp -f "${RES_DIR}/00_extractSNPsfromVCF.py" "${SCRIPTS_DIR}/00_extractSNPsfromVCF.py"
  cp -f "${RES_DIR}/filterRepetiveSNP.py" "${SCRIPTS_DIR}/filterRepetiveSNP.py"
  # Replace ntsm-scripts makefile with our full pipeline makefile
  cp -f "${RES_DIR}/makefile.biomatch" "${SCRIPTS_DIR}/makefile"
  chmod +x "${SCRIPTS_DIR}/filterRepetiveSNP.py" || true
  chmod +x "${SCRIPTS_DIR}/00_extractSNPsfromVCF.py" || true
  echo "[post-link] ntsm-scripts patched with optimized Python workflow."
else
  echo "[post-link] ntsm-scripts directory not found; skipping patch."
fi