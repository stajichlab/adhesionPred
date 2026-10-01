"""ESM-2 embedding extraction for protein sequences."""

import esm
import numpy as np
import torch
from tqdm import tqdm

# Global model cache to avoid reloading models
_MODEL_CACHE = {}

# Layer the shipped models were trained on; see get_esm_embeddings.
DEFAULT_REPR_LAYER = 6
# ESM-2 context is 1024 tokens including BOS and EOS.
MAX_RESIDUES = 1022
POOLING = "residue_mean"
# Residue characters the ESM-2 alphabet accepts; anything else becomes 'X'.
ESM_RESIDUES = set("ACDEFGHIKLMNPQRSTVWYXBUZO")
STRIP_CHARS = "*-."


def get_cached_model(model_name, device=None):
    """Get cached ESM-2 model or load and cache if not present.

    Args:
        model_name: ESM-2 model variant to use.
        device: torch device (cuda or cpu).

    Returns:
        Tuple of (model, alphabet).
    """
    if device is None:
        device = torch.device("cuda" if torch.cuda.is_available() else "cpu")

    cache_key = (model_name, str(device))

    if cache_key not in _MODEL_CACHE:
        print(f"Loading ESM-2 model: {model_name} on {device}...")
        model, alphabet = esm.pretrained.load_model_and_alphabet(model_name)
        model = model.eval()
        model = model.to(device)
        _MODEL_CACHE[cache_key] = (model, alphabet)
        print(f"Model {model_name} cached successfully")
    else:
        print(f"Using cached ESM-2 model: {model_name} on {device}")
        model, alphabet = _MODEL_CACHE[cache_key]

    return model, alphabet


def get_optimal_batch_size(device, model_name):
    """Determine optimal batch size based on GPU memory and model size.

    Args:
        device: torch device.
        model_name: ESM-2 model variant name.

    Returns:
        Optimal batch size for the given configuration.
    """
    if device.type == "cuda":
        try:
            total_memory = torch.cuda.get_device_properties(device).total_memory
            # Estimate batch size based on model size and available memory
            if "t6_8M" in model_name:
                # Smaller model, can handle larger batches
                optimal_size = min(32, max(4, total_memory // (2 * 1024**3)))
            elif "t12_35M" in model_name:
                # Larger model, needs smaller batches
                optimal_size = min(16, max(2, total_memory // (4 * 1024**3)))
            else:
                # Unknown model, use conservative estimate
                optimal_size = min(8, max(2, total_memory // (3 * 1024**3)))

            print(f"Optimized batch size for {model_name} on {device}: {optimal_size}")
            return optimal_size
        except Exception as e:
            print(f"Could not determine optimal batch size: {e}, using default")
            return 8
    else:
        # CPU processing, use smaller batches
        return 4


def sanitize_sequence(sequence):
    """Make a protein sequence safe for the ESM-2 tokenizer.

    Upper-cases, drops stop/gap characters ('*', '-', '.') and maps any other
    character outside ESM-2's residue alphabet (e.g. 'J') to 'X', so that no
    sequence can make a batch fail to tokenize.
    """
    seq = sequence.upper()
    for ch in STRIP_CHARS:
        seq = seq.replace(ch, "")
    return "".join(ch if ch in ESM_RESIDUES else "X" for ch in seq)


def _embed_batch(model, alphabet, batch_converter, batch, repr_layer, device):
    """Embed one batch of (id, sequence) pairs; mean over residue tokens only.

    BOS, EOS and padding positions are excluded from the mean, so a sequence's
    embedding does not depend on the other sequences in its batch.
    """
    _, _, tokens = batch_converter(batch)
    tokens = tokens.to(device)
    with torch.no_grad():
        results = model(tokens, repr_layers=[repr_layer], return_contacts=False)
    reps = results["representations"][repr_layer]
    mask = (
        (tokens != alphabet.padding_idx)
        & (tokens != alphabet.cls_idx)
        & (tokens != alphabet.eos_idx)
    ).unsqueeze(-1)
    pooled = (reps * mask).sum(dim=1) / mask.sum(dim=1).clamp(min=1)
    return pooled.float().cpu().numpy()


def get_esm_embeddings(
    sequences,
    model_name="esm2_t6_8M_UR50D",
    batch_size=None,
    device=None,
    repr_layer=DEFAULT_REPR_LAYER,
    return_indices=False,
):
    """Extract ESM-2 embeddings for protein sequences.

    Sequences are sanitized (see sanitize_sequence), truncated to MAX_RESIDUES,
    batched in length order, and mean-pooled over residue tokens only. Results
    are returned in input order. A batch that fails is retried one sequence at
    a time; a sequence that still fails is skipped with a warning, so callers
    that need labels must align them by the returned ids or indices.

    Args:
        sequences: List of sequence dictionaries with 'id' and 'sequence' keys.
        model_name: ESM-2 model variant to use.
        batch_size: Number of sequences to process at once. If None, will auto-optimize.
        device: torch device (cuda or cpu).
        repr_layer: Transformer layer whose representations are pooled. Defaults to 6,
            the layer the shipped models were trained on (the final layer of the 6-layer
            model, the middle layer of the 12-layer model).
        return_indices: If True, also return the input positions of the embedded sequences.

    Returns:
        Tuple of (embeddings array, sequence ids), plus the list of input indices
        when return_indices is True.
    """
    if device is None:
        device = torch.device("cuda" if torch.cuda.is_available() else "cpu")

    # Use cached model loading
    model, alphabet = get_cached_model(model_name, device)
    if not 0 < repr_layer <= model.num_layers:
        raise ValueError(f"repr_layer={repr_layer} but {model_name} has {model.num_layers} layers")

    # Auto-optimize batch size if not provided
    if batch_size is None:
        batch_size = get_optimal_batch_size(device, model_name)

    print(f"Using batch size: {batch_size}")

    batch_converter = alphabet.get_batch_converter()

    cleaned = [sanitize_sequence(seq["sequence"])[:MAX_RESIDUES] for seq in sequences]
    # Length-sorted batching minimizes padding; results are put back in input order.
    order = sorted(range(len(sequences)), key=lambda k: len(cleaned[k]))
    pooled = {}

    print(f"Extracting embeddings for {len(sequences)} sequences...")
    for i in tqdm(range(0, len(order), batch_size)):
        idx = order[i : i + batch_size]
        batch = [(str(k), cleaned[k]) for k in idx]
        try:
            for k, vec in zip(
                idx,
                _embed_batch(model, alphabet, batch_converter, batch, repr_layer, device),
                strict=False,
            ):
                pooled[k] = vec
        except Exception as e:
            print(f"Warning: batch failed ({e}); retrying its {len(idx)} sequences one at a time")
            for k in idx:
                try:
                    pooled[k] = _embed_batch(
                        model, alphabet, batch_converter, [(str(k), cleaned[k])], repr_layer, device
                    )[0]
                except Exception as e2:
                    print(f"Warning: skipping sequence {sequences[k]['id']}: {e2}")

    kept = sorted(pooled)
    if len(kept) < len(sequences):
        print(
            f"Warning: {len(sequences) - len(kept)} of {len(sequences)} sequences were not embedded"
        )
    embeddings = np.array([pooled[k] for k in kept])
    seq_ids = [sequences[k]["id"] for k in kept]
    if return_indices:
        return embeddings, seq_ids, kept
    return embeddings, seq_ids


ESM2_MODEL_CHOICES = [
    "esm2_t6_8M_UR50D",
    "esm2_t12_35M_UR50D",
]
