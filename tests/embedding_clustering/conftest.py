"""Make analysis/embedding_clustering importable as plain modules for tests."""

import sys
from pathlib import Path

EMBEDDING_CLUSTERING_DIR = Path(__file__).resolve().parents[2] / "analysis" / "embedding_clustering"
sys.path.insert(0, str(EMBEDDING_CLUSTERING_DIR))
