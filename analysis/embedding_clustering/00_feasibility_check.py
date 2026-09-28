"""Task 1 feasibility spike: can we load and embed with both ESM2 and ESM C?

This script is a record of the working feasibility check, not a reusable library
module. Findings are written up in task-1-report.md; Task 2 should read that
report before writing embed.py.

IMPORTANT ENVIRONMENT NOTE (see report for details):
- Step 1 (ESM2 via `fair-esm`) runs in the project's normal Python 3.9 environment.
- Step 2 (ESM C via the `esm` PyPI package) requires Python >=3.10, which this
  project's env (requires-python = ">=3.9", currently running 3.9.18) does not
  provide. Step 2 was therefore run in a separate local venv created with:
      uv venv --python 3.11 .venv_esmc
      uv pip install --python .venv_esmc/bin/python esm torch httpx
  This script documents both parts; Step 2's code must be *run* under that venv
  (or an equivalent Python >=3.10 environment), not under the project's normal env.
"""

import sys

import torch

# ---------------------------------------------------------------------------
# Step 1: confirm the classifier's own ESM2 embedding path still works.
# Run this part with the project's normal Python (3.9) environment.
# ---------------------------------------------------------------------------


def step1_confirm_esm2():
    sys.path.insert(0, "src")
    from adhesion_predict.embeddings import get_esm_embeddings

    sequences = [
        {
            "id": "test1",
            "sequence": (
                "MSTNPKPQRKTKRNTNRRPQDVKFPGGGQIVGGVYLLPRRGPRLGVRATRKTSERSQPRGRRQPIPKARRPEG"
                "RTWAQPGYPWPLYGNEGCGWAGWLLSPRGSRPSWGPTDPRRRSRNLGKVIDTLTCGFADLMGYIPLVGAPLGG"
                "AARALAHGVRVLEDGVNYATGNLPGCSFSIFLLALLSCLTVPASA"
            ),
        },
        {
            "id": "test2",
            "sequence": (
                "MKTAYIAKQRQISFVKSHFSRQLEERLGLIEVQAPILSRVGDGTQDNLSGAEKAVQVKVKALPDAQFEVVHSL"
                "AKWKRQTLGQHDFSAGEGLYTHMKALRPDEDRLSPLHSVYVDQWDWELVMGDGERQFSTLKSTVEAIWAGIK"
                "ATEAAVSEEFGLAPFLPDQIHFVHSQELLSRYPDLDAKGRERAIAKDLGAVFLVGIGGKLSDGHRHDVRAPD"
                "YDDWSTPSELGHAGLNGDILVWNPVLEDAFELSSMGIRVDADTLKHQLALTGDEDRLELEWHQALLRGEMPQ"
                "TIGGGIGQSRLTMLLLQLPHIGQVQAGVWPAAVRESVPSLL"
            ),
        },
    ]
    embeddings, ids = get_esm_embeddings(sequences, model_name="esm2_t12_35M_UR50D")
    print(f"ESM2 shape: {embeddings.shape}, expected (2, 480)")
    assert embeddings.shape == (2, 480), f"Unexpected shape: {embeddings.shape}"
    print("Classifier (ESM2) embedding path: OK")
    return embeddings, ids


# ---------------------------------------------------------------------------
# Step 2: ESM C 300M via the `esm` PyPI package (Approach A).
# This must be run with a Python >=3.10 interpreter (see venv note above) --
# it will NOT import under the project's normal Python 3.9 environment.
# ---------------------------------------------------------------------------


def step2_confirm_esmc():
    from esm.models.esmc import ESMC
    from esm.sdk.api import ESMProtein, LogitsConfig

    model = ESMC.from_pretrained("esmc_300m").eval()
    tok = model.tokenizer
    cls_id, eos_id, pad_id = tok.cls_token_id, tok.eos_token_id, tok.pad_token_id

    sequences = [
        "MSTNPKPQRKTKRNTNRRPQDVKFPGGGQIVGG",
        (
            "MKTAYIAKQRQISFVKSHFSRQLEERLGLIEVQAPILSRVGDGTQDNLSGAEKAVQVKVKALPDAQFEVVHSL"
            "AKWKRQTLGQHDFSAGEGLYTHMKALRPDEDRLSPLHSVYVDQWDWELVMGDGERQFSTLKSTVEAIWAGIK"
            "ATEAAVSEEFGLAPFLPDQIHFVHSQELLSRYPDLDAKGRERAIAKDLGAVFLVGIGGKLSDGHRHDVRAPD"
            "YDDWSTPSELGHAGLNGDILVWNPVLEDAFELSSMGIRVDADTLKHQLALTGDEDRLELEWHQALLRGEMPQ"
            "TIGGGIGQSRLTMLLLQLPHIGQVQAGVWPAAVRESVPSLL"
        ),
    ]

    pooled_list = []
    for s in sequences:
        protein = ESMProtein(sequence=s)
        protein_tensor = model.encode(protein)
        out = model.logits(protein_tensor, LogitsConfig(sequence=True, return_embeddings=True))
        emb = out.embeddings[0]  # (L, 960)
        tokens = protein_tensor.sequence  # (L,)

        # Special-token-aware mean pooling: exclude BOS (<cls>), EOS (<eos>),
        # and PAD from the average -- only residue-token positions contribute.
        mask = (tokens != cls_id) & (tokens != eos_id) & (tokens != pad_id)
        pooled = emb[mask].mean(dim=0)
        assert mask.sum().item() == len(s), "pooling mask did not match sequence length"
        pooled_list.append(pooled)

    result = torch.stack(pooled_list)
    print(f"ESM C shape: {result.shape}, expected (2, 960)")
    assert result.shape == (2, 960), f"Unexpected shape: {result.shape}"
    print("ESM C embedding path (Approach A, `esm` package): OK")
    return result.numpy(), sequences


if __name__ == "__main__":
    step1_confirm_esm2()
    # step2_confirm_esmc() must be invoked under a Python >=3.10 interpreter,
    # e.g.: .venv_esmc/bin/python -c "from analysis.embedding_clustering.00_feasibility_check import step2_confirm_esmc; step2_confirm_esmc()"
