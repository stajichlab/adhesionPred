"""Tests for ESM-2 embedding extraction (uses the small esm2_t6_8M model)."""

import numpy as np
import pytest

pytest.importorskip("esm")

from surface_glyco import embeddings as emb  # noqa: E402
from surface_glyco.embeddings import get_esm_embeddings, sanitize_sequence  # noqa: E402

SHORT = "MKTLLVAGLLSSAAFA"
LONG = "MSTTSSTTSTPSSTSA" * 40


def _seqs(*pairs):
    return [{"id": i, "sequence": s} for i, s in pairs]


def test_sanitize_sequence_maps_unknown_characters():
    assert sanitize_sequence("mkJt*a-.Ö") == "MKXTAX"


def test_embedding_does_not_depend_on_batch_mates():
    alone, _ = get_esm_embeddings(_seqs(("a", SHORT)), batch_size=1, device=_cpu())
    batched, ids = get_esm_embeddings(_seqs(("a", SHORT), ("b", LONG)), batch_size=2, device=_cpu())
    assert ids == ["a", "b"]
    np.testing.assert_allclose(alone[0], batched[0], atol=1e-5)


def test_results_are_in_input_order_despite_length_sorting():
    seqs = _seqs(("long", LONG), ("short", SHORT), ("mid", SHORT * 5))
    emb, ids, kept = get_esm_embeddings(seqs, batch_size=2, device=_cpu(), return_indices=True)
    assert ids == ["long", "short", "mid"]
    assert kept == [0, 1, 2]
    single, _ = get_esm_embeddings([seqs[2]], batch_size=1, device=_cpu())
    np.testing.assert_allclose(emb[2], single[0], atol=1e-5)


def test_nonstandard_characters_do_not_drop_sequences():
    seqs = _seqs(("ok", SHORT), ("odd", "mkjt*" + SHORT), ("ok2", SHORT[::-1]))
    emb, ids = get_esm_embeddings(seqs, batch_size=3, device=_cpu())
    assert ids == ["ok", "odd", "ok2"]
    assert emb.shape[0] == 3


def test_repr_layer_out_of_range_raises():
    with pytest.raises(ValueError):
        get_esm_embeddings(_seqs(("a", SHORT)), device=_cpu(), repr_layer=7)


def test_count_truncated_counts_sequences_over_the_limit():
    seqs = _seqs(("a", "M" * 1022), ("b", "M" * 1023), ("c", ""), ("d", "M*" * 1200))
    assert emb.count_truncated(seqs) == 2  # b, and d (1200 residues once '*' is stripped)


def test_failed_sequence_is_skipped_and_indices_stay_aligned(monkeypatch):
    real = emb._embed_batch

    def flaky(model, alphabet, bc, batch, layer, device):
        if any(seq == "BAD" for _, seq in batch):
            raise RuntimeError("boom")
        return real(model, alphabet, bc, batch, layer, device)

    monkeypatch.setattr(emb, "_embed_batch", flaky)
    seqs = _seqs(("a", "MKTAYIAK"), ("bad", "BAD"), ("c", "MKTAYIAKQR"))
    vecs, ids, kept = emb.get_esm_embeddings(seqs, device=_cpu(), batch_size=3, return_indices=True)
    assert ids == ["a", "c"] and kept == [0, 2] and len(vecs) == 2


def _cpu():
    import torch

    return torch.device("cpu")
