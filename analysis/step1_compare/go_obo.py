"""Parse go-basic.obo: parents over is_a and part_of, obsolete terms (spec 3.1 step 1).

Standard library only. The rules copy /tmp/glyco_spec/d1_count.py: only [Term] stanzas count,
parents are `is_a` and `relationship: part_of`, other relationships are ignored, and alt_id
lines are recorded but not used for ancestors.
"""

from dataclasses import dataclass, field
from pathlib import Path


@dataclass
class Ontology:
    parents: dict[str, set[str]] = field(default_factory=dict)
    names: dict[str, str] = field(default_factory=dict)
    obsolete: set[str] = field(default_factory=set)
    alt_ids: dict[str, str] = field(default_factory=dict)
    data_version: str = ""
    _cache: dict[str, frozenset[str]] = field(default_factory=dict, repr=False)

    def known(self, term: str) -> bool:
        return term in self.names

    def ancestors(self, term: str) -> frozenset[str]:
        """The term itself plus every is_a / part_of ancestor. Unknown terms give {term}."""
        cached = self._cache.get(term)
        if cached is not None:
            return cached
        seen = {term}
        stack = [term]
        while stack:
            for parent in self.parents.get(stack.pop(), ()):
                if parent not in seen:
                    seen.add(parent)
                    stack.append(parent)
        result = frozenset(seen)
        self._cache[term] = result
        return result


def parse_obo(path: str | Path) -> Ontology:
    onto = Ontology()
    current = None
    in_term = False
    with open(path, encoding="utf-8") as handle:
        for raw in handle:
            line = raw.rstrip("\n")
            if line.startswith("["):
                in_term = line == "[Term]"
                current = None
            elif not in_term and current is None and line.startswith("data-version:"):
                onto.data_version = line.split(":", 1)[1].strip()
            elif in_term and line.startswith("id: GO:"):
                current = line.split()[1]
                onto.names.setdefault(current, "")
                onto.parents.setdefault(current, set())
            elif current is None:
                continue
            elif line.startswith("name: "):
                onto.names[current] = line[6:]
            elif line.startswith("is_a: "):
                onto.parents[current].add(line[6:].split()[0])
            elif line.startswith("relationship: part_of "):
                onto.parents[current].add(line.split()[2])
            elif line == "is_obsolete: true":
                onto.obsolete.add(current)
            elif line.startswith("alt_id: "):
                onto.alt_ids[line[8:].split()[0]] = current
    return onto
