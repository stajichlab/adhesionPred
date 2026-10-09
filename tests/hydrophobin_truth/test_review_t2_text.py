import importlib.util
from pathlib import Path

SCRIPT = Path(__file__).resolve().parents[2] / "analysis/hydrophobin_truth/review_t2_text.py"
spec = importlib.util.spec_from_file_location("review_t2_text", SCRIPT)
m = importlib.util.module_from_spec(spec)
spec.loader.exec_module(m)


def test_related_texts():
    assert (
        m.classify("Aerial growth, conidiation, and dispersal of filamentous fungi")
        == "hydrophobin_related"
    )
    assert (
        m.classify("Self-assembles to form functional amyloid fibrils called rodlets")
        == "hydrophobin_related"
    )
    assert m.classify("Lowers the surface tension of water; amphipathic") == "hydrophobin_related"
    assert m.classify("Spore wall hydrophobin") == "hydrophobin_related"


def test_other_texts():
    assert m.classify("Transcription factor that regulates development") == "other"
    assert m.classify("") == "other"


def test_entry_text_uses_only_experimental_function_location_subunit_comments():
    entry = {
        "comments": [
            {
                "commentType": "FUNCTION",
                "texts": [
                    {
                        "value": "Forms rodlets",
                        "evidences": [
                            {"evidenceCode": "ECO:0000269", "source": "PubMed", "id": "1"}
                        ],
                    }
                ],
            },
            {
                "commentType": "FUNCTION",
                "texts": [
                    {
                        "value": "Inferred by similarity",
                        "evidences": [{"evidenceCode": "ECO:0000250"}],
                    }
                ],
            },
            {
                "commentType": "INDUCTION",
                "texts": [
                    {
                        "value": "x",
                        "evidences": [
                            {"evidenceCode": "ECO:0000269", "source": "PubMed", "id": "2"}
                        ],
                    }
                ],
            },
        ]
    }
    assert m.experimental_text(entry) == "Forms rodlets"
