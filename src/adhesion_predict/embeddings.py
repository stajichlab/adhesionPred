"""ESM-2 embedding extraction for protein sequences."""

import esm
import numpy as np
import torch
from tqdm import tqdm


def get_esm_embeddings(sequences, model_name="esm2_t6_8M_UR50D", batch_size=8, device=None):
    """Extract ESM-2 embeddings for protein sequences.

    Args:
        sequences: List of sequence dictionaries with 'id' and 'sequence' keys.
        model_name: ESM-2 model variant to use.
        batch_size: Number of sequences to process at once.
        device: torch device (cuda or cpu).

    Returns:
        Tuple of (embeddings array, sequence ids).
    """
    if device is None:
        device = torch.device("cuda" if torch.cuda.is_available() else "cpu")

    print(f"Loading ESM-2 model: {model_name}...")
    model, alphabet = esm.pretrained.load_model_and_alphabet(model_name)
    model = model.eval()
    model = model.to(device)

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
