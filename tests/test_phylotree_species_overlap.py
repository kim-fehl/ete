from ete4 import PhyloTree

# Tree used by the tests provided by the ``phylotree_example_newick`` fixture.
#   ╭───┬╴Dme_001
#   │   ╰╴Dme_002
#   │       ╭───┬╴Cfa_001
#───┤       │   ╰╴Mms_001
#   │   ╭───┤           ╭───┬╴Hsa_001
#   │   │   │       ╭───┤   ╰╴Hsa_003
#   │   │   │   ╭───┤   ╰╴Ptr_001
#   ╰───┤   ╰───┤   ╰╴Mmu_001
#       │       │   ╭───┬╴Hsa_004
#       │       ╰───┤   ╰╴Ptr_004
#       │           ╰╴Mmu_004
#       ╰───┬╴Ptr_002
#           ╰───┬╴Hsa_002
#               ╰╴Mmu_002

def test_get_sp_overlap_on_all_descendants(phylotree_example_newick, sp_name_fn):
    """Tests orthology prediction using the species overlap algorithm"""
    # Creates a gene phylogeny with several duplication events at
    # different levels.
    t = PhyloTree(phylotree_example_newick, sp_naming_function=sp_name_fn)
    # Scans the tree using the species overlap algorithm and detect all
    # speciation and duplication events
    events = t.get_descendant_evol_events()

    # Check that all duplications are detected
    dup1 = t.common_ancestor(["Hsa_001", "Hsa_004"])
    assert dup1.props.get('evoltype') == "D"

    dup2 = t.common_ancestor(["Dme_001", "Dme_002"])
    assert dup2.props.get('evoltype') == "D"

    dup3 = t.common_ancestor(["Hsa_001", "Hsa_002"])
    assert dup3.props.get('evoltype') == "D"

    dup4 = t.common_ancestor(["Hsa_001", "Hsa_003"])
    assert dup4.props.get('evoltype') == "D"

    # All other nodes should be speciation
    for node in t.traverse():
        if not node.is_leaf and node not in {dup1, dup2, dup3, dup4}:
            assert node.props.get('evoltype') == "S"

    # Check events
    for e in events:
        assert e.node.props.get('evoltype') == e.etype

    # Check orthology/paralogy prediction
    orthologs = set()
    for e in events:
        if e.node == dup1:
            assert e.inparalogs == {'Ptr_001', 'Hsa_001', 'Mmu_001', 'Hsa_003'}
            assert e.outparalogs == {'Mmu_004', 'Ptr_004', 'Hsa_004'}
            assert e.orthologs == set()
            assert e.outparalogs == e.out_seqs
            assert e.inparalogs == e.in_seqs
        elif e.node == dup2:
            assert e.inparalogs == {'Dme_001'}
            assert e.outparalogs == {'Dme_002'}
            assert e.orthologs == set()
            assert e.outparalogs == e.out_seqs
            assert e.inparalogs == e.in_seqs
        elif e.node == dup3:
            assert e.inparalogs == {'Hsa_003', 'Cfa_001', 'Ptr_001', 'Hsa_001', 'Ptr_004', 'Hsa_004', 'Mmu_004', 'Mmu_001', 'Mms_001'}
            assert e.outparalogs == {'Hsa_002', 'Ptr_002', 'Mmu_002'}
            assert e.orthologs == set()
            assert e.outparalogs == e.out_seqs
            assert e.inparalogs == e.in_seqs
        elif e.node == dup4:
            assert e.inparalogs == {'Hsa_001'}
            assert e.outparalogs == {'Hsa_003'}
            assert e.orthologs == set()
            assert e.outparalogs == e.out_seqs
            assert e.inparalogs == e.in_seqs
        else:
            key1 = sorted(e.inparalogs)
            key2 = sorted(e.orthologs)
            orthologs.add(tuple(sorted([tuple(key1), tuple(key2)])))

    orthologies = [
        [
            {'Dme_001', 'Dme_002'},
            {'Ptr_001', 'Cfa_001', 'Hsa_002', 'Hsa_003', 'Ptr_002', 'Hsa_001', 'Ptr_004', 'Hsa_004', 'Mmu_004', 'Mmu_001', 'Mms_001', 'Mmu_002'},
        ],
        [
            {'Mms_001', 'Cfa_001'},
            {'Hsa_003', 'Ptr_001', 'Hsa_001', 'Ptr_004', 'Hsa_004', 'Mmu_004', 'Mmu_001'},
        ],
        [{'Ptr_002'}, {'Hsa_002', 'Mmu_002'}],
        [{'Cfa_001'}, {'Mms_001'}],
        [{'Hsa_002'}, {'Mmu_002'}],
        [{'Hsa_003', 'Hsa_001', 'Ptr_001'}, {'Mmu_001'}],
        [{'Ptr_004', 'Hsa_004'}, {'Mmu_004'}],
        [{'Hsa_003', 'Hsa_001'}, {'Ptr_001'}],
        [{'Hsa_004'}, {'Ptr_004'}],
    ]
    expected_orthologs = set()
    for l1, l2 in orthologies:
        key1 = sorted(l1)
        key2 = sorted(l2)
        expected_orthologs.add(tuple(sorted([tuple(key1), tuple(key2)])))

    # Are all orthologies as expected
    assert expected_orthologs == orthologs

    # Test different sos_thr
    t = PhyloTree('(((SP1_a, SP2_a), (SP3_a, SP1_b)), (SP1_c, SP2_c));', sp_naming_function=sp_name_fn)
    seed = t['SP1_a']
    events = t.get_descendant_evol_events(0.1)
    assert t.common_ancestor([seed, 'SP3_a']).props.get('evoltype') == 'D'
    assert t.common_ancestor([seed, 'SP1_c']).props.get('evoltype') == 'D'

    t = PhyloTree('(((SP1_a, SP2_a), (SP3_a, SP1_b)), (SP1_c, SP2_c));', sp_naming_function=sp_name_fn)
    seed = t['SP1_a']
    events = t.get_descendant_evol_events(0.5)
    assert t.common_ancestor([seed, 'SP3_a']).props.get('evoltype') == 'S'
    assert t.common_ancestor([seed, 'SP1_c']).props.get('evoltype') == 'D'

    t = PhyloTree('(((SP1_a:1, SP2_a:1):1, (SP3_a:1, SP1_b:1):1):1, (SP1_c:1, SP2_c:1):1):0;', sp_naming_function=sp_name_fn)
    seed = t['SP1_a']
    events = seed.get_my_evol_events(0.75)
    assert t.common_ancestor([seed, 'SP3_a']).props.get('evoltype') == 'S'
    assert t.common_ancestor([seed, 'SP1_c']).props.get('evoltype') == 'S'


def test_get_sp_overlap_on_a_seed(phylotree_example_newick, sp_name_fn):
    """Tests orthology prediction using species overlap from a seed node"""
    # Creates a gene phylogeny with several duplication events at
    # different levels.
    t = PhyloTree(phylotree_example_newick, sp_naming_function=sp_name_fn)
    # Scans the tree using the species overlap algorithm
    seed = t['Hsa_001']
    events = seed.get_my_evol_events()

    # Check that duplications are detected
    dup1 = t.common_ancestor(["Hsa_001", "Hsa_004"])
    assert dup1.props.get('evoltype') == "D"

    # This duplication is not in the seed path
    dup2 = t.common_ancestor(["Dme_001", "Dme_002"])
    assert not hasattr(dup2, "evoltype")

    dup3 = t.common_ancestor(["Hsa_001", "Hsa_002"])
    assert dup3.props.get('evoltype') == "D"

    dup4 = t.common_ancestor(["Hsa_001", "Hsa_003"])
    assert dup4.props.get('evoltype') == "D"

    # All other nodes should be speciation
    node = seed
    while node:
        if not node.is_leaf and node not in {dup1, dup2, dup3, dup4}:
            assert node.props.get('evoltype') == "S"
        node = node.up

    # Check events
    for e in events:
        assert e.node.props.get('evoltype') == e.etype

    # Check orthology/paralogy prediction
    orthologs = set()
    for e in events:
        if e.node == dup1:
            assert e.inparalogs == {'Hsa_001', 'Hsa_003'}
            assert e.outparalogs == {'Hsa_004'}
            assert e.orthologs == set()
            assert e.in_seqs == {'Ptr_001', 'Hsa_001', 'Mmu_001', 'Hsa_003'}
            assert e.out_seqs == {'Mmu_004', 'Ptr_004', 'Hsa_004'}
        elif e.node == dup3:
            assert e.inparalogs == {'Hsa_003', 'Hsa_001', 'Hsa_004'}
            assert e.outparalogs == {'Hsa_002'}
            assert e.orthologs == set()
            assert e.in_seqs == {'Hsa_003', 'Cfa_001', 'Ptr_001', 'Hsa_001', 'Ptr_004', 'Hsa_004', 'Mmu_004', 'Mmu_001', 'Mms_001'}
            assert e.out_seqs == {'Hsa_002', 'Ptr_002', 'Mmu_002'}
        elif e.node == dup4:
            assert e.inparalogs == {'Hsa_001'}
            assert e.outparalogs == {'Hsa_003'}
            assert e.orthologs == set()
            assert e.in_seqs == {'Hsa_001'}
            assert e.out_seqs == {'Hsa_003'}
        else:
            key1 = sorted(e.inparalogs)
            key2 = sorted(e.orthologs)
            orthologs.add(tuple(sorted([tuple(key1), tuple(key2)])))

    # Are all orthologies as expected
    orthologies = [
        [{'Dme_001', 'Dme_002'}, {'Hsa_002', 'Hsa_003', 'Hsa_001', 'Hsa_004'}],
        [{'Mms_001', 'Cfa_001'}, {'Hsa_003', 'Hsa_001', 'Hsa_004'}],
        [{'Hsa_003', 'Hsa_001'}, {'Mmu_001'}],
        [{'Hsa_003', 'Hsa_001'}, {'Ptr_001'}],
    ]
    expected_orthologs = set()
    for l1, l2 in orthologies:
        key1 = sorted(l1)
        key2 = sorted(l2)
        expected_orthologs.add(tuple(sorted([tuple(key1), tuple(key2)])))

    assert expected_orthologs == orthologs

    # Test different sos_thr
    t = PhyloTree('(((SP1_a:1, SP2_a:1):1, (SP3_a:1, SP1_b:1):1):1, (SP1_c:1, SP2_c:1):1):0;', sp_naming_function=sp_name_fn)
    seed = t['SP1_a']
    events = seed.get_my_evol_events(0.1)
    assert t.common_ancestor([seed, 'SP3_a']).props.get('evoltype') == 'D'
    assert t.common_ancestor([seed, 'SP1_c']).props.get('evoltype') == 'D'

    t = PhyloTree('(((SP1_a:1, SP2_a:1):1, (SP3_a:1, SP1_b:1):1):1, (SP1_c:1, SP2_c:1):1):0;', sp_naming_function=sp_name_fn)
    seed = t['SP1_a']
    events = seed.get_my_evol_events(0.50)
    assert t.common_ancestor([seed, 'SP3_a']).props.get('evoltype') == 'S'
    assert t.common_ancestor([seed, 'SP1_c']).props.get('evoltype') == 'D'

    t = PhyloTree('(((SP1_a:1, SP2_a:1):1, (SP3_a:1, SP1_b:1):1):1, (SP1_c:1, SP2_c:1):1):0;', sp_naming_function=sp_name_fn)
    seed = t['SP1_a']
    events = seed.get_my_evol_events(0.75)
    assert t.common_ancestor([seed, 'SP3_a']).props.get('evoltype') == 'S'
    assert t.common_ancestor([seed, 'SP1_c']).props.get('evoltype') == 'S'
