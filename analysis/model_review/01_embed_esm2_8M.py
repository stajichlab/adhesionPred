"""Embed training set + S288C proteome with ESM-2 t6 8M two ways:
  legacy  = exactly what src/adhesion_predict/embeddings.py does (batch 4 on CPU, file order,
            mean over ALL token positions incl. BOS/EOS/padding)
  masked  = mean over residue tokens only (batch-independent)
"""
import gzip, glob, sys, numpy as np, pandas as pd, torch, esm
from Bio import SeqIO

from pathlib import Path
REPO = str(Path(__file__).resolve().parents[2])
torch.set_num_threads(8)

def read(path):
    h = gzip.open(path, "rt") if path.endswith(".gz") else open(path)
    return [(r.id, str(r.seq).replace("J", "L").replace("*", "")) for r in SeqIO.parse(h, "fasta")]

rows = []
# same order as io.load_sequences_from_dir: extension loop, glob order
for lab, d in [(1, "positive"), (0, "negative")]:
    for f in sorted(glob.glob(f"{REPO}/data/{d}/*.pep.gz")):
        src = f.split("/")[-1].split(".")[0]
        for i, s in read(f):
            rows.append(dict(id=i, seq=s, label=lab, source=src))
train = pd.DataFrame(rows)
scer = pd.DataFrame([dict(id=i, seq=s) for i, s in read(sys.argv[1]) if len(s) > 0])

model, alphabet = esm.pretrained.esm2_t6_8M_UR50D()
model.eval()
bc = alphabet.get_batch_converter()

def embed(seqs, bs=4):
    legacy, masked = [], []
    for k in range(0, len(seqs), bs):
        data = [(str(j), s[:1022]) for j, s in enumerate(seqs[k:k + bs])]
        _, _, toks = bc(data)
        with torch.no_grad():
            rep = model(toks, repr_layers=[6])["representations"][6]
        legacy.append(rep.mean(1).numpy())
        m = ((toks != alphabet.padding_idx) & (toks != alphabet.cls_idx) & (toks != alphabet.eos_idx)).unsqueeze(-1)
        masked.append(((rep * m).sum(1) / m.sum(1)).numpy())
        if k % 400 == 0:
            print(k, len(seqs), flush=True)
    return np.vstack(legacy), np.vstack(masked)

L, M = embed(train.seq.tolist())
np.save("train_legacy.npy", L); np.save("train_masked.npy", M)
train.to_csv("train_meta.csv", index=False)
L, M = embed(scer.seq.tolist())
np.save("scer_legacy.npy", L); np.save("scer_masked.npy", M)
scer.to_csv("scer_meta.csv", index=False)
print("done")
