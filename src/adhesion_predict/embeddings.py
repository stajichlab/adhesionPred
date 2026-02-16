"""ESM-2 embedding extraction for protein sequences."""

import esm
import numpy as np
import torch
from tqdm import tqdm

# Global model cache to avoid reloading models
_MODEL_CACHE = {}


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


def get_esm_embeddings(sequences, model_name="esm2_t6_8M_UR50D", batch_size=None, device=None):
    """Extract ESM-2 embeddings for protein sequences.

    Args:
        sequences: List of sequence dictionaries with 'id' and 'sequence' keys.
        model_name: ESM-2 model variant to use.
        batch_size: Number of sequences to process at once. If None, will auto-optimize.
        device: torch device (cuda or cpu).

    Returns:
        Tuple of (embeddings array, sequence ids).
    """
    if device is None:
        device = torch.device("cuda" if torch.cuda.is_available() else "cpu")

    # Use cached model loading
    model, alphabet = get_cached_model(model_name, device)

    # Auto-optimize batch size if not provided
    if batch_size is None:
        batch_size = get_optimal_batch_size(device, model_name)

    print(f"Using batch size: {batch_size}")

    batch_converter = alphabet.get_batch_converter()

    embeddings = []
    seq_ids = []

    print(f"Extracting embeddings for {len(sequences)} sequences...")
    for i in tqdm(range(0, len(sequences), batch_size)):
        batch_seqs = sequences[i : i + batch_size]
        batch_data = [(seq["id"], seq["sequence"][:1022]) for seq in batch_seqs]

        try:
            batch_labels, batch_strs, batch_tokens = batch_converter(batch_data)
            batch_tokens = batch_tokens.to(device)

            with torch.no_grad():
                results = model(batch_tokens, repr_layers=[6], return_contacts=False)

            token_representations = results["representations"][6]
            sequence_representations = token_representations.mean(dim=1)

            for j in range(len(batch_seqs)):
                embeddings.append(sequence_representations[j].cpu().numpy())
                seq_ids.append(batch_seqs[j]["id"])

        except Exception as e:
            print(f"Warning: Failed to process batch starting at {i}: {e}")
            continue

    return np.array(embeddings), seq_ids


ESM2_MODEL_CHOICES = [
    "esm2_t6_8M_UR50D",
    "esm2_t12_35M_UR50D",
]
