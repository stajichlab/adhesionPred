#!/usr/bin/python3.12
"""Fetch the fungal part of two public allergen sources and summarise it (issue #19 scoping).

Sources
  WHO/IUIS Allergen Nomenclature  http://www.allergen.org/csv.php?table=allergen | joint
      free to download; the site asks users to cite the site URL and a recent IUIS publication
  UniProt keyword KW-0020 (Allergen), taxonomy 4751 (Fungi), via the UniProt REST API

Outputs (in --out, default analysis/allergen_scoping):
  iuis_fungal_allergens.tsv       one row per IUIS fungal allergen molecule
  iuis_fungal_isoallergens.tsv    one row per isoallergen, with sequence and accessions
  uniprot_kw0020_fungi.tsv        UniProt fungal entries with the Allergen keyword, with features
  iuis_fungal_uniprot_features.tsv  UniProt features for the IUIS accessions
  summary.txt                     the numbers quoted in docs/reports/2026-10-04-fungal-allergen-scoping.md

Usage: /usr/bin/python3.12 01_fetch_and_summarize.py [--out DIR]
"""

import argparse
import io
import urllib.parse
import urllib.request
from pathlib import Path

import pandas as pd

FIELDS = "accession,organism_name,protein_name,length,ft_signal,ft_lipid,cc_subcellular_location,xref_pfam"


def get(url):
    with urllib.request.urlopen(url, timeout=120) as r:
        return r.read().decode()


def uniprot(query, fields=FIELDS):
    url = f"https://rest.uniprot.org/uniprotkb/stream?format=tsv&fields={fields}&query={urllib.parse.quote(query)}"
    return pd.read_csv(io.StringIO(get(url)), sep="\t").fillna("")


def flags(t):
    t["signal_peptide"] = t["Signal peptide"] != ""
    t["gpi"] = t["Lipidation"].str.contains("GPI", case=False)
    t["secreted_loc"] = t["Subcellular location [CC]"].str.contains("ecreted", na=False)
    t["intracellular_loc"] = t["Subcellular location [CC]"].str.contains(
        "ytoplasm|itochond|eroxisom|ucleus|ndoplasm|acuol", na=False
    )
    return t


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--out", default="analysis/allergen_scoping")
    out = Path(ap.parse_args().out)
    out.mkdir(parents=True, exist_ok=True)
    log = []

    a = pd.read_csv(io.StringIO(get("http://www.allergen.org/csv.php?table=allergen")))
    j = pd.read_csv(
        io.StringIO(get("http://www.allergen.org/csv.php?table=joint")), low_memory=False
    )
    fa = a[a.TaxSource.str.contains("ung", na=False)]
    fj = j[j.TaxSource.str.contains("ung", na=False)]
    fa.to_csv(out / "iuis_fungal_allergens.tsv", sep="\t", index=False)
    fj[
        [
            "TaxSource",
            "TaxOrder",
            "Species",
            "AllergenID",
            "Name",
            "BioNames",
            "IsoName",
            "AccProtein",
            "AccUniProt",
            "AccPDB",
            "Sequence",
            "SeqFeatures",
            "IsoAllergenicity",
            "SequenceRef",
        ]
    ].to_csv(out / "iuis_fungal_isoallergens.tsv", sep="\t", index=False)
    log.append(f"IUIS all allergens {len(a)}; by TaxSource {a.TaxSource.value_counts().to_dict()}")
    log.append(
        f"IUIS fungal molecules {len(fa)}; species {fa.Species.nunique()}; orders {fa.TaxOrder.value_counts().to_dict()}"
    )
    log.append(f"IUIS fungal exposure {fa.Exposure.value_counts().to_dict()}")
    log.append(f"IUIS fungal genera {fa.Species.str.split().str[0].value_counts().to_dict()}")
    log.append(
        f"IUIS fungal isoallergen rows {len(fj)}; with protein sequence {fj.Sequence.notna().sum()}; "
        f"with UniProt accession {fj.AccUniProt.notna().sum()}"
    )

    k = flags(uniprot("keyword:KW-0020 AND taxonomy_id:4751"))
    k.to_csv(out / "uniprot_kw0020_fungi.tsv", sep="\t", index=False)
    log.append(
        f"UniProt KW-0020 fungi {len(k)}; signal peptide {k.signal_peptide.sum()}; GPI {k.gpi.sum()}; "
        f"secreted loc {k.secreted_loc.sum()}; intracellular loc {k.intracellular_loc.sum()}"
    )
    log.append(
        f"UniProt KW-0020 fungi genera {k.Organism.str.split().str[0].value_counts().to_dict()}"
    )

    acc = sorted(set(fj.AccUniProt.dropna().astype(str).str.split(r"[;, ]+").explode()))
    frames = []
    for i in range(0, len(acc), 20):
        q = "(" + " OR ".join(f"accession:{x}" for x in acc[i : i + 20]) + ")"
        frames.append(uniprot(q))
    t = flags(pd.concat(frames))
    t.to_csv(out / "iuis_fungal_uniprot_features.tsv", sep="\t", index=False)
    log.append(
        f"IUIS UniProt accessions {len(acc)}; found {len(t)}; signal peptide {t.signal_peptide.sum()}; GPI {t.gpi.sum()}; "
        f"secreted loc {t.secreted_loc.sum()}; intracellular loc {t.intracellular_loc.sum()}; "
        f"no location annotation {(t['Subcellular location [CC]'] == '').sum()}"
    )
    log.append(
        f"overlap IUIS accessions and KW-0020 fungi set {len(set(acc) & set(k.Entry))}; KW-0020 only {len(set(k.Entry) - set(acc))}"
    )
    (out / "summary.txt").write_text("\n".join(log) + "\n")
    print("\n".join(log))


if __name__ == "__main__":
    main()
