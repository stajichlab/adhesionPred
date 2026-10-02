"""Constants shared by the J2 runner, the chunk store and the assembly. Standard library only.

MODEL_DIM is the width of the pooled embedding of each model. A test compares it with the real
model configuration (needs torch and fair-esm). REPR_LAYER is the layer that the model card
rules use (`surface_glyco.card.DEFAULT_REPR_LAYER`).
"""

from surface_glyco.card import DEFAULT_REPR_LAYER

MODEL_DIM = {"esm2_t6_8M_UR50D": 320, "esm2_t12_35M_UR50D": 480}
MODELS = tuple(MODEL_DIM)
REPR_LAYER = DEFAULT_REPR_LAYER
