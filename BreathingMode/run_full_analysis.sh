#!/usr/bin/env bash
set -euo pipefail

ROOT_DIR="$(cd "$(dirname "$0")" && pwd)/.."
cd "$ROOT_DIR"

mkdir -p BreathingMode/analysis_results

echo "Starting full extraction: logs -> BreathingMode/analysis_results/extract.log"
./.venv/bin/python3 -u BreathingMode/extract_breathing.py > BreathingMode/analysis_results/extract.log 2>&1

echo "Extraction finished. Starting analysis (nboot=2000): logs -> BreathingMode/analysis_results/analysis.log"
./.venv/bin/python3 -u BreathingMode/analysis_run.py --nboot 2000 > BreathingMode/analysis_results/analysis.log 2>&1

echo "Run complete: $(date)" > BreathingMode/analysis_results/run_complete.flag
echo "All done. Results: BreathingMode/analysis_results/analysis_summary.json"
