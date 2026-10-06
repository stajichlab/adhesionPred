"""NCBI taxonomy lineage lookups (from nodes.dmp) for status matching."""

from pathlib import Path


class TaxonError(ValueError):
    """A taxon ID is not in the taxonomy."""


class Lineage:
    def __init__(self, parent):
        self._parent = dict(parent)

    @classmethod
    def from_nodes_dmp(cls, path):
        """Read ``nodes.dmp`` (fields separated by ``\\t|\\t``; field 0 = taxid, field 1 = parent)."""
        parent = {}
        text = Path(path).read_text(encoding="utf-8-sig", errors="replace")
        for n, line in enumerate(text.splitlines(), 1):
            fields = line.split("\t|\t")
            if len(fields) < 2:
                continue
            try:
                parent[int(fields[0])] = int(fields[1].replace("\t|", "").strip())
            except ValueError:
                raise TaxonError(f"{path}:{n}: taxon IDs are not integers") from None
        return cls(parent)

    def ancestors(self, taxon):
        """The taxon first, then its parents up to the root."""
        taxon = int(taxon)
        if taxon not in self._parent:
            raise TaxonError(f"taxon {taxon} is not in the taxonomy")
        chain = [taxon]
        while True:
            up = self._parent[chain[-1]]
            if up == chain[-1] or up in chain:  # root, or a cycle: stop
                return chain
            if up not in self._parent:
                raise TaxonError(f"parent {up} of taxon {chain[-1]} is not in the taxonomy")
            chain.append(up)

    def is_descendant_or_self(self, taxon, ancestor):
        return int(ancestor) in self.ancestors(taxon)

    def depth(self, taxon):
        return len(self.ancestors(taxon)) - 1
