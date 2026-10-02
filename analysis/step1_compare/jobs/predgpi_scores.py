#!/usr/bin/env python
"""PredGPI (module predgpi/202001) with its continuous scores, one TSV row per protein.

The PredGPI CLI (`predgpi.py -m gff3`) writes only a class: a `GPI-anchor` line with score
1.0 (FPR <= 0.0015), 0.70 (<= 0.005) or 0.55 (<= 0.01), else a `Chain` line. This wrapper
calls the same functions of predgpi.py (predGpipe, the HMM and SVM files of the installation)
and also writes the estimated false positive rate and the SVM output, which the hybrid
candidate H needs as a GPI score. The class rules and residue substitutions copy
predgpi.py main(); test_predgpi_wrapper_matches_cli checks the classes against the CLI.

Runs under the module's Python 3.9 (`module load predgpi/202001`; PREDGPI_HOME must be set).
Keep this file Python 3.9 compatible. Columns: id, length, gpi_call, gpi_prob, omega, fpr, svm.
"""

import argparse
import gzip
import os
import sys

COLUMNS = ("id", "length", "gpi_call", "gpi_prob", "omega", "fpr", "svm")


TOO_SHORT_MAX = 40  # predgpi.py main(): a sequence of 40 residues or fewer is not scored


def read_fasta(path):
    """(id, sequence) records. Spaces inside a sequence are removed, as predgpi.py does."""
    opener = gzip.open if path.endswith(".gz") else open
    name, chunks = None, []
    with opener(path, "rt") as handle:
        for line in handle:
            line = line.strip()
            if line.startswith(">"):
                if name is not None:
                    yield name, "".join(chunks).replace(" ", "")
                name, chunks = line[1:].split()[0], []
            elif line:
                chunks.append(line)
    if name is not None:
        yield name, "".join(chunks).replace(" ", "")


def classify(fpr):
    if fpr <= 0.0015:
        return "highly_probable", "1.0"
    if fpr <= 0.005:
        return "probable", "0.70"
    if fpr <= 0.01:
        return "weakly", "0.55"
    return "none", "0"


def score_row(name, seq, predict):
    """One TSV line (without newline) for a protein.

    predict(seq_t) returns (cut, svmout, fpr) for the sequence after the residue
    substitutions of predgpi.py. It is not called for a sequence of TOO_SHORT_MAX residues
    or fewer."""
    if len(seq) <= TOO_SHORT_MAX:
        return "\t".join([name, str(len(seq)), "too_short", "0", "", "", ""])
    seq_t = seq.replace("U", "C").replace("Z", "A").replace("B", "A").replace("X", "A")
    cut, svmout, fpr = predict(seq_t)
    call, prob = classify(fpr)
    omega = str(len(seq) - cut) if call != "none" else ""
    return "\t".join(
        [name, str(len(seq)), call, prob, omega, format(fpr, ".6g"), format(svmout, ".6g")]
    )


def write_scores(fasta, out_path, predict):
    """Write the scores to out_path through out_path + '.tmp'. Return 0, or 2 after STOP.

    On any error the temporary file is removed, `STOP: <reason>` goes to stderr and out_path
    is not written."""
    tmp = out_path + ".tmp"
    seen = set()
    try:
        with open(tmp, "w") as out:
            out.write("\t".join(COLUMNS) + "\n")
            for name, seq in read_fasta(fasta):
                if name in seen:
                    raise ValueError("duplicate id " + name)
                seen.add(name)
                out.write(score_row(name, seq, predict) + "\n")
        os.replace(tmp, out_path)
    except Exception as exc:  # any scoring error must leave no output
        if os.path.exists(tmp):
            os.remove(tmp)
        print(f"STOP: {type(exc).__name__}: {exc}", file=sys.stderr)
        return 2
    return 0


def main(argv=None):
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--fasta", required=True)
    parser.add_argument("--out", required=True)
    args = parser.parse_args(argv)
    home = os.environ.get("PREDGPI_HOME")
    if not home:
        print("STOP: PREDGPI_HOME is not set; run `module load predgpi/202001`", file=sys.stderr)
        return 2
    sys.path.insert(0, home)
    import predgpi  # lazy: only the PredGPI module Python has it
    from predgpilib.hmm import HMM_IO
    from predgpilib.svm import SVMLike

    hmm = HMM_IO.get_hmm(os.path.join(home, "GPIDAT", "PHMM.TOT.ss.mod"))
    svm = SVMLike.getSVMLight(os.path.join(home, "GPIDAT", "MOD"))

    def predict(seq_t):
        _, cut, svmout, fpr = predgpi.predGpipe(seq_t, svm, hmm)
        return cut, svmout, fpr

    return write_scores(args.fasta, args.out, predict)


if __name__ == "__main__":
    sys.exit(main())
