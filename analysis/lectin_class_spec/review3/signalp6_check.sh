#!/usr/bin/env bash
# Review 3: re-run SignalP 6 (fast, eukarya) on the 6 AAL-like sequences plus 2 controls (Ecp6 B3VBK9, Emp46 Q12396).
# The module signalp/6 points to the 6.0i package; its binary reports version 6.0h. The binary is called
# directly because `module load` needs conda init on the login shell.
# Usage: signalp6_check.sh WORK_DIR   (WORK_DIR holds aal.faa: the 8 sequences from unilectin_fungal_sequences.faa)
set -euo pipefail
D=${1:?work dir}
S=/opt/linux/rocky/8.x/x86_64/pkgs/signalp/6.0i/env/bin/signalp6
"$S" --version
"$S" --fastafile "$D/aal.faa" --organism eukarya --mode fast --output_dir "$D/sp6" --format txt
cat "$D/sp6/prediction_results.txt"
