def test_ancestors_follow_is_a_and_part_of(mini_ontology):
    assert mini_ontology.ancestors("GO:0009277") == {
        "GO:0009277",
        "GO:0005618",
        "GO:0110165",
        "GO:0005575",
        "GO:0071944",
    }
    assert "GO:0012505" in mini_ontology.ancestors("GO:0005783")  # ER part_of endomembrane
    assert {"GO:0016020", "GO:0005739"} <= mini_ontology.ancestors("GO:0005743")


def test_other_relationships_and_typedefs_are_ignored(mini_ontology):
    assert "GO:0005886" not in mini_ontology.ancestors("GO:0032991")  # has_part ignored
    assert mini_ontology.parents["GO:0008150"] == set()  # Typedef is_a not attached


def test_obsolete_alt_id_and_version(mini_ontology):
    assert mini_ontology.obsolete == {"GO:0031225"}
    assert mini_ontology.alt_ids == {"GO:0099999": "GO:0009277"}
    assert mini_ontology.ancestors("GO:0099999") == {"GO:0099999"}  # as d1_count.py
    assert mini_ontology.data_version == "releases/2026-07-26"
    assert mini_ontology.known("GO:0005576")
    assert not mini_ontology.known("GO:1234567")
