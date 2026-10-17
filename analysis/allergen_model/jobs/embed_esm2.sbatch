#!/bin/bash -l
# T3.1: ESM-2 650M embeddings, mean-pooled over residue tokens only (padding, BOS and EOS excluded).
# Uses the repo's function src/surface_glyco/embeddings.py get_esm_embeddings(), with the final layer passed explicitly.
# ESM-2 650M has 33 layers (esm2_t33). The function's default repr_layer is 6, so the layer is set here, not left to the default.
# Environment: the conda env adhesionPred (torch 2.10 with CUDA, fair-esm), at /rhome/jstajich/.conda/envs/adhesionPred.
# Input:  $SEQS (FASTA for one protein set), $ARM, $OUTDIR
# Output: $OUTDIR/embeddings.npz (ids, matrix)
#         $OUTDIR/embed.done  (written last: model, layer, sequence count, truncated count, input hash, torch version, GPU name)
# Submit: sbatch --export=ALL,SEQS=...,OUTDIR=...,ARM=... embed_esm2.sbatch
# Partition: exfab (owner decision O4). One node with 2x ada6000 (48 GB each). The account is in the exfab association.
#SBATCH -p exfab
#SBATCH --gres=gpu:1
#SBATCH -c 8
#SBATCH --mem=64G
#SBATCH --time=4:00:00
#SBATCH -J ap_embed
set -euo pipefail
: "${SEQS:?export SEQS (protein set FASTA) before sbatch}"
: "${OUTDIR:?export OUTDIR before sbatch}"
: "${ARM:?export ARM before sbatch}"
TMP="${SCRATCH:?SCRATCH is not set; run this as a SLURM job}/ap_embed.$$"
trap 'rm -rf "$TMP"' EXIT
mkdir -p "$TMP" "$OUTDIR"
rm -f "$OUTDIR/embed.done"
# Pre-check: the account must be in the exfab association. The field is cut from the sacctmgr pipe-delimited output.
if ! sacctmgr -n -p show assoc user="$USER" format=account | cut -d'|' -f1 | grep -qx "exfab"; then
  echo "FATAL: $USER is not in the exfab account (plan O4); resubmit to another queue or request access" >&2
  exit 1
fi
PY=/rhome/jstajich/.conda/envs/adhesionPred/bin/python
[ -x "$PY" ] || { echo "FATAL: conda env adhesionPred not found at $PY" >&2; exit 1; }
ESM_MODEL="${ESM_MODEL:-esm2_t33_650M_UR50D}"
ESM_LAYER="${ESM_LAYER:-33}"
export PYTHONPATH="/bigdata/stajichlab/jstajich/projects/adhesionPred/src${PYTHONPATH:+:$PYTHONPATH}"
"$PY" - "$SEQS" "$ESM_MODEL" "$ESM_LAYER" "$TMP/embeddings.npz" "$TMP/embed_info.txt" <<'PY'
import sys, hashlib, torch
import numpy as np
from surface_glyco.embeddings import get_esm_embeddings, count_truncated

fasta, model_name, layer, out, info_path = sys.argv[1], sys.argv[2], int(sys.argv[3]), sys.argv[4], sys.argv[5]
if not torch.cuda.is_available():
    raise SystemExit("FATAL: torch sees no GPU; the job would run on CPU for hours")
records, cur_id, parts = [], None, []
with open(fasta) as fh:
    for line in fh:
        line = line.strip()
        if line.startswith(">"):
            if cur_id is not None:
                records.append({"id": cur_id, "sequence": "".join(parts)})
            cur_id, parts = line[1:].split()[0], []
        elif line:
            parts.append(line)
if cur_id is not None:
    records.append({"id": cur_id, "sequence": "".join(parts)})
if not records:
    raise SystemExit(f"FATAL: {fasta} has no sequences")
matrix, ids = get_esm_embeddings(records, model_name=model_name, repr_layer=layer)
if len(ids) != len(records):
    raise SystemExit(f"FATAL: {len(records) - len(ids)} sequences were not embedded; the set is incomplete")
np.savez(out, ids=np.array(ids), matrix=np.asarray(matrix))
trunc = count_truncated(records)
gpu = torch.cuda.get_device_name(0)
sha = hashlib.sha256(open(fasta, "rb").read()).hexdigest()
with open(info_path, "w") as fh:
    fh.write(f"model={model_name}\nrepr_layer={layer}\nsequences={len(ids)}\ntruncated={trunc}\n"
             f"input_sha256={sha}\ntorch={torch.__version__}\ngpu={gpu}\n")
print(f"embedded={len(ids)} truncated={trunc} layer={layer} gpu={gpu}")
PY
[ -s "$TMP/embeddings.npz" ] || { echo "FATAL: no embeddings file written" >&2; exit 1; }
cp "$TMP/embeddings.npz" "$OUTDIR/.tmp.embeddings.npz" && mv "$OUTDIR/.tmp.embeddings.npz" "$OUTDIR/embeddings.npz"
{ cat "$TMP/embed_info.txt"; echo "arm=$ARM"; echo "date=$(date +%F)"; } > "$TMP/embed.done"
cp "$TMP/embed.done" "$OUTDIR/.tmp.embed.done" && mv "$OUTDIR/.tmp.embed.done" "$OUTDIR/embed.done"
