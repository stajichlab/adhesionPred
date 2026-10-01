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


def read_fasta(path):
    opener = gzip.open if path.endswith(".gz") else open
    name, chunks = None, []
    with opener(path, "rt") as handle:
        for line in handle:
            line = line.strip()
            if line.startswith(">"):
                if name is not None:
                    yield name, "".join(chunks)
                name, chunks = line[1:].split()[0], []
            elif line:
                chunks.append(line)
    if name is not None:
        yield name, "".join(chunks)


def classify(fpr):
    if fpr <= 0.0015:
        return "highly_probable", "1.0"
    if fpr <= 0.005:
        return "probable", "0.70"
    if fpr <= 0.01:
        return "weakly", "0.55"
    return "none", "0"


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
    import predgpi
    from predgpilib.hmm import HMM_IO
    from predgpilib.svm import SVMLike

    hmm = HMM_IO.get_hmm(os.path.join(home, "GPIDAT", "PHMM.TOT.ss.mod"))
    svm = SVMLike.getSVMLight(os.path.join(home, "GPIDAT", "MOD"))
    tmp = args.out + ".tmp"
    seen = set()
    with open(tmp, "w") as out:
        out.write("\t".join(COLUMNS) + "\n")
        for name, seq in read_fasta(args.fasta):
            if name in seen:
                print("STOP: duplicate id " + name, file=sys.stderr)
                os.remove(tmp)
                return 2
            seen.add(name)
            if len(seq) <= 40:
                out.write("\t".join([name, str(len(seq)), "too_short", "0", "", "", ""]) + "\n")
                continue
            seq_t = seq.replace("U", "C").replace("Z", "A").replace("B", "A").replace("X", "A")
            _, cut, svmout, fpr = predgpi.predGpipe(seq_t, svm, hmm)
            call, prob = classify(fpr)
            omega = str(len(seq) - cut) if call != "none" else ""
            row = [name, str(len(seq)), call, prob, omega, "%.6g" % fpr, "%.6g" % svmout]
            out.write("\t".join(row) + "\n")
    os.replace(tmp, args.out)
    return 0


if __name__ == "__main__":
    sys.exit(main())
