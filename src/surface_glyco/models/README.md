# Models

No trained model ships with this package yet. The earlier models (trained on FLO/ALS homologs
versus random proteins, with padding-inclusive pooling) were removed; they remain in git history.
Train one with `surface_glyco_train`, or pass `--model path/to/model.pkl`. A model is accompanied
by a JSON card (`model.json`) that `surface_glyco_predict` reads. Set `SURFACE_GLYCO_MODELS_DIR`
to use a different directory.
