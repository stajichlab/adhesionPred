"""Per-sequence inputs of the Phase C models: features, composition, embeddings. Needs numpy.

One row per unique sequence of Phase B (`row` of features_unique.tsv.gz and of
emb/<model>.nterm.npy). The C-terminal candidates (M8-C, M35-C) read emb/<model>.cterm.npy row
`cterm_row` for proteins longer than 1,022 aa and the N-terminal row for shorter proteins
(Phase C spec 3.3).
"""

import hashlib
from dataclasses import dataclass, field
from pathlib import Path

import evalio
import numpy as np
import rules
import truth_table

AA20 = "ACDEFGHIKLMNPQRSTVWY"
MAX_RESIDUES = 1022
EMBEDDINGS = {
    "M8": ("esm2_t6_8M_UR50D", False),
    "M35": ("esm2_t12_35M_UR50D", False),
    "M8-C": ("esm2_t6_8M_UR50D", True),
    "M35-C": ("esm2_t12_35M_UR50D", True),
}
MODELS = ("esm2_t6_8M_UR50D", "esm2_t12_35M_UR50D")


class ModelError(ValueError):
    """An input of the Phase C models is unusable (re-exported as models.ModelError)."""


@dataclass
class Universe:
    hashes: list[str]
    length: np.ndarray
    sp: np.ndarray
    sp_prob: np.ndarray
    rank: np.ndarray
    gpi_prob: np.ndarray
    st: np.ndarray
    comp: np.ndarray
    cterm_row: np.ndarray
    nterm: dict = field(default_factory=dict)
    cterm: dict = field(default_factory=dict)
    index: dict = field(default_factory=dict)

    def __post_init__(self):
        self.index = {h: i for i, h in enumerate(self.hashes)}
        if len(self.index) != len(self.hashes):
            raise ModelError("the universe lists a hash more than once")
        check_cterm_rows(self.length, self.cterm_row)

    def rows(self, hashes) -> np.ndarray:
        try:
            return np.array([self.index[h] for h in hashes], dtype=np.int64)
        except KeyError as exc:
            raise evalio.StopError(f"hash {exc.args[0]} is not a Phase B unique sequence") from exc


def check_cterm_rows(length, cterm_row) -> None:
    long = np.asarray(length) > MAX_RESIDUES
    has = np.asarray(cterm_row) >= 0
    bad = np.nonzero(long != has)[0]
    if len(bad):
        raise evalio.StopError(
            f"row {bad[0]}: length {int(length[bad[0]])} and cterm_row {int(cterm_row[bad[0]])} "
            f"disagree (a C-terminal row exists exactly for proteins over {MAX_RESIDUES} aa)"
        )


def composition(sequence: str) -> list[float]:
    """Fraction of each of the 20 amino acids. The divisor is the full length, X included,
    so the fractions sum to less than 1 for a sequence with X or other symbols."""
    n = len(sequence)
    if n == 0:
        raise ModelError("composition of an empty sequence is undefined")
    return [sequence.count(a) / n for a in AA20]


def embedding(u: Universe, name: str, idx) -> np.ndarray:
    """Embedding rows of candidate `name` (M8, M35, M8-C, M35-C) for universe rows idx."""
    model, use_cterm = EMBEDDINGS[name]
    X = np.asarray(u.nterm[model][idx], dtype=np.float64)
    if use_cterm:
        crow = u.cterm_row[idx]
        long = crow >= 0
        if long.any():
            X[long] = np.asarray(u.cterm[model][crow[long]], dtype=np.float64)
    return X


def features(u: Universe, candidate: str, idx, h_variant: str | None = None) -> np.ndarray:
    idx = np.asarray(idx, dtype=np.int64)
    loglen = np.log(u.length[idx].astype(np.float64))[:, None]
    if candidate == "B0":
        return loglen
    if candidate == "B1":
        return np.hstack([u.comp[idx], loglen])
    if candidate in EMBEDDINGS:
        return embedding(u, candidate, idx)
    if candidate == "H":
        if h_variant not in EMBEDDINGS:
            raise ModelError(
                f"candidate H needs h_variant in {tuple(EMBEDDINGS)}, got {h_variant!r}"
            )
        extra = np.c_[u.sp_prob[idx], u.gpi_prob[idx], u.st[idx]]
        return np.hstack([embedding(u, h_variant, idx), extra])
    raise ValueError(f"candidate {candidate!r} has no feature matrix")


def array_sha256(arr) -> str:
    return hashlib.sha256(np.ascontiguousarray(arr).tobytes()).hexdigest()


def load(phaseb: Path, verify: bool = True) -> Universe:
    """Read features_unique, unique_sequences and the four embedding matrices of Phase B."""
    phaseb = Path(phaseb)
    fu = truth_table.read_tsv(phaseb / "features_unique.tsv.gz")
    seqs = truth_table.read_tsv(phaseb / "unique_sequences.tsv.gz")
    if [r["seq_sha256"] for r in fu] != [r["seq_sha256"] for r in seqs]:
        raise evalio.StopError("features_unique and unique_sequences list other hashes")
    for i, r in enumerate(fu):
        if r["row"] != str(i):
            raise evalio.StopError(f"features_unique.tsv.gz row {i} has row value {r['row']}")
    num = {}
    for col in ("sp_prob", "gpi_prob", "ser_thr_frac"):
        num[col] = np.array([evalio.float_or_stop(r[col], f"{r['seq_sha256']} {col}") for r in fu])
    run = evalio.read_json(phaseb / "emb" / "embedding_run.json")
    u = Universe(
        hashes=[r["seq_sha256"] for r in fu],
        length=np.array([int(r["length"]) for r in fu], dtype=np.int64),
        sp=np.array([r["sp_prediction"] == "SP" for r in fu], dtype=bool),
        sp_prob=num["sp_prob"],
        rank=rules.gpi_rank([r["gpi_call"] for r in fu]),
        gpi_prob=num["gpi_prob"],
        st=num["ser_thr_frac"],
        comp=np.array([composition(r["sequence"]) for r in seqs], dtype=np.float64),
        cterm_row=np.array([int(r["cterm_row"]) if r["cterm_row"] else -1 for r in fu]),
    )
    n_long = int((u.cterm_row >= 0).sum())
    for model in MODELS:
        for window, store, rows in (("nterm", u.nterm, len(fu)), ("cterm", u.cterm, n_long)):
            arr = np.load(phaseb / "emb" / f"{model}.{window}.npy", mmap_mode="r")
            want = run.get("models", {}).get(model, {}).get(window, {})
            if list(arr.shape) != want.get("shape") or arr.shape[0] != rows:
                raise evalio.StopError(
                    f"emb/{model}.{window}.npy has shape {arr.shape}; embedding_run.json says "
                    f"{want.get('shape')} and the features need {rows} rows"
                )
            if verify and array_sha256(arr) != want.get("array_sha256"):
                raise evalio.StopError(f"emb/{model}.{window}.npy differs from embedding_run.json")
            store[model] = arr
    return u
