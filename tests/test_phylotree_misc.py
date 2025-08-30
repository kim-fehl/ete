import pytest
from ete4 import PhyloTree


def test_miscelaneus(phylotree_example_newick, sp_name_fn):
    """Miscellaneous PhyloTree features"""
    # Creates a gene phylogeny with several duplication events at
    # different levels.
    t = PhyloTree(phylotree_example_newick, sp_naming_function=sp_name_fn)
    # Create a dictionary with relative ages for the species present in
    # the phylogenetic tree. Note that ages are only relative numbers to
    # define which species are older, and that different species can
    # belong to the same age.
    sp2age = {
        'Hsa': 1,  # Homo sapiens (Hominids)
        'Ptr': 2,  # P. troglodytes (primates)
        'Mmu': 2,  # Macaca mulata (primates)
        'Mms': 3,  # Mus musculus (mammals)
        'Cfa': 3,  # Canis familiaris (mammals)
        'Dme': 4,  # Drosophila melanogaster (metazoa)
    }

    # Check that dup ages are correct
    dup1 = t.common_ancestor(["Hsa_001", "Hsa_004"])
    assert dup1.get_age(sp2age) == 2
    dup2 = t.common_ancestor(["Dme_001", "Dme_002"])
    assert dup2.get_age(sp2age) == 4
    dup3 = t.common_ancestor(["Hsa_001", "Hsa_002"])
    assert dup3.get_age(sp2age) == 3
    dup4 = t.common_ancestor(["Hsa_001", "Hsa_003"])
    assert dup4.get_age(sp2age) == 1

    # Check rooting options
    expected_root = t['Dme_002']
    expected_root.dist += 2.3
    assert t.get_farthest_oldest_leaf(sp2age) == expected_root
    #print t
    #print t.get_farthest_oldest_node(sp2age)

    # Check get species functions
    assert t.get_species() == set(sp2age.keys())
    assert set(sp for sp in t.iter_species()) == set(sp2age.keys())


def test_collapse(sp_name_fn):
    t = PhyloTree('((Dme_001,Dme_002),(((Cfa_001,Mms_001),((((Hsa_001,Hsa_001),Ptr_001),Mmu_001),((Hsa_004,Ptr_004),Mmu_004))),(Ptr_002,(Hsa_002,Mmu_002))));',
                  sp_naming_function=sp_name_fn)
    for n in t.traverse():
        n.dist = 1
        n.support = 1
    collapsed_hsa = '((Dme_001:1,Dme_002:1)1:1,(((Cfa_001:1,Mms_001:1)1:1,(((Ptr_001:1,Hsa_001:1)1:1,Mmu_001:1)1:1,((Hsa_004:1,Ptr_004:1)1:1,Mmu_004:1)1:1)1:1)1:1,(Ptr_002:1,(Hsa_002:1,Mmu_002:1)1:1)1:1)1:1);'
    t2 = t.collapse_lineage_specific_expansions(['Hsa'])
    assert str(collapsed_hsa) == str(t2.write(props=["species"], parser=2))
    with pytest.raises(TypeError):
        print(t.collapse_lineage_specific_expansions('Hsa'))
