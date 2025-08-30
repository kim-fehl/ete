import pytest

from ete4 import Tree, PhyloTree
from ete4.core.tree import TreeError

from .conftest import ROBINSON_FOULDS_CASES, COPHENETIC_MATRIX_CASES


def test_distances():
    # Distances: get_distance, get_farthest_node,
    # get_farthest_descendant, get_midpoint_outgroup
    t = Tree('(((A:0.1, B:0.01):0.001, C:0.0001):1.0[&&NHX:name=I], '
             '(D:0.00001):0.000001[&&NHX:name=J]):2.0[&&NHX:name=root];')
    A = t['A']
    B = t['B']
    C = t['C']
    D = t['D']
    I = t['I']
    J = t['J']
    root = t['root']

    def assertApprox(x, y):
        assert abs(x - y) < 1e-8  # approximately equal numbers

    assert t.common_ancestor([A, I]).name == "I"
    assert t.common_ancestor([A, D]).name == "root"
    assertApprox(t.get_distance(A, I), 0.101)
    assertApprox(t.get_distance(A, B), 0.11)
    assertApprox(t.get_distance(A, A), 0)
    assertApprox(t.get_distance(I, I), 0)
    assertApprox(t.get_distance(A, root), t.get_distance(root, A))

    assertApprox(t.get_distance(A, root), t.get_distance(root, A))
    assertApprox(t.get_distance(root, A), t.get_distance(A, root))

    # Get_farthest_node, get_farthest_leaf
    assert root.get_farthest_leaf() == (A, 1.101)
    assert root.get_farthest_node() == (A, 1.101)
    assert A.get_farthest_leaf() == (A, 0.0)
    assert A.get_farthest_node() == (D, 1.101011)
    assert I.get_farthest_node() == (D, 1.000011)

    # Topology only distances
    t = Tree('(((A:0.5, B:1.0):1.0, C:5.0):1, (D:10.0, F:1.0):2.0):20;')

    assert t.get_closest_leaf() == (t['A'], 2.5)
    assert t.get_farthest_leaf() == (t['D'], 12.0)
    assert t.get_farthest_leaf(topological=True) == (t['A'], 2.0)
    assert t.get_closest_leaf(topological=True) == (t['C'], 1.0)
    assert t.get_distance(t, t) == 0.0
    assert t.get_distance(t, t, topological=True) == 0.0
    assert t.get_distance(t, t['A'], topological=True) == 3.0

    assert (t['F']).get_farthest_node(topological=True) == (t['A'], 3.0)
    assert (t['F']).get_farthest_node(topological=False) == (t['D'], 11.0)


def test_ultrametric(rng):
    EPSILON = 1e-5  # small number for the purposes of comparing distances

    # Convert tree to a ultrametric, in which the distance from
    # leafs to root is always the same.
    t = Tree()
    t.populate(80, dist_fn=rng.random, support_fn=rng.random)
    max_dist = max(t.get_distance(t, n) for n in t)

    t.to_ultrametric()
    assert all(abs(t.get_distance(t, n) - max_dist) < EPSILON for n in t)

    t2 = Tree()
    t2.populate(80, dist_fn=rng.random, support_fn=rng.random)
    max_dist = max(t2.get_distance(t2, n) for n in t2)

    t2.to_ultrametric(topological=True)
    assert all(abs(t2.get_distance(t2, n) - max_dist) < EPSILON for n in t2)
    leaf, _ = t2.get_farthest_leaf(topological=True)
    assert all(abs(node.dist - leaf.dist) < EPSILON
               for node in leaf.ancestors() if not node.is_root)


def test_expand_polytomies_rf():
    gtree = Tree('((a:1, (b:1, (c:1, d:1):1):1), (e:1, (f:1, g:1):1):1);')
    ref1 = Tree('((a:1, (b:1, c:1, d:1):1):1, (e:1, (f:1, g:1):1):1);')
    ref2 = Tree('((a:1, (b:1, c:1, d:1):1):1, (e:1, f:1, g:1):1);')

    for ref in [ref1, ref2]:
        gtree.robinson_foulds(ref, expand_polytomies=True)[0]

    gtree = Tree('((g, h), (a, (b, (c, (d,( e, f))))));')
    ref3 = Tree('((a, b, c, (d, e, f)), (g, h));')
    ref4 = Tree('((a, b, c, d, e, f), (g, h));')
    ref5 = Tree('((a, b, (c, d, (e, f))), (g, h));')

    for ref in [ref3, ref4, ref5]:
        gtree.robinson_foulds(ref, expand_polytomies=True,
                              polytomy_size_limit=8)[0]

    gtree = Tree('((g, h), (a, b, (c, d, (e, f))));')
    ref6 = Tree('((a, b, (c, d, e, f)), (g, h));')
    ref7 = Tree('((a, (b, (c, d, e, f))), (g, h));')
    ref8 = Tree('((a, b, c, (d, e, f)), (g, h));')
    ref9 = Tree('((d, b, c, (a, e, f)), (g, h));')

    for ref in [ref6, ref7, ref8, ref9]:
        gtree.robinson_foulds(ref, expand_polytomies=True)[0]

    gtree = Tree('((g, h), ((a, b), (c, d), (e, f)));')
    ref10 = Tree('((g, h), ((a, c), ((b, d), (e, f))));')

    for ref in [ref10]:
        gtree.robinson_foulds(ref, expand_polytomies=True,
                              polytomy_size_limit=8)[0]


def test_robinson_foulds_exceptions():
    t1 = Tree('(a,b,(c,d,e));')
    t2 = Tree('((a,b),(c,d,e));')
    # testing unrooted trees
    with pytest.raises(TreeError):
        t1.robinson_foulds(t2=t2)

    # expand polytomies and unrooted trees
    with pytest.raises(TreeError):
        t1.robinson_foulds(t2=t2,
                      unrooted_trees=True, expand_polytomies=True)

    # using expand_polytomies and correct_by_size at the same time
    with pytest.raises(TreeError):
        t1.robinson_foulds(t2=t1,
                      unrooted_trees=True, expand_polytomies=True,
                      correct_by_polytomy_size=True)

    # correct by size when polytomies in both sides
    with pytest.raises(TreeError):
        t1.robinson_foulds(t2=t1,
                      unrooted_trees=True, correct_by_polytomy_size=True)

    # polytomy larger than deafult limit
    with pytest.raises(TreeError):
        t2.robinson_foulds(t2=Tree('(a, (b,c,d,e,f,g,h));'),
                      expand_polytomies=True)

    # duplicated items
    t3 = Tree('(a, (b, (c, c)));')
    with pytest.raises(TreeError):
        t3.robinson_foulds(t2=t2)
    with pytest.raises(TreeError):
        t2.robinson_foulds(t2=t3)


@pytest.mark.parametrize("RF, unrooted, nw1, nw2", ROBINSON_FOULDS_CASES)
def test_robinson_foulds_and_more(rng, RF, unrooted, nw1, nw2):
    t1 = Tree(nw1)
    t2 = Tree(nw2)
    rf, rf_max, names, r1, r2, d1, d2 = t1.robinson_foulds(t2, unrooted_trees=unrooted)
    real_max = (20*2) - 4 if not unrooted else (20*2) - 6

    assert len(names) == 20
    assert rf_max == real_max
    assert rf == RF

    comp = t1.compare(t2, unrooted=unrooted)
    assert 20 == comp['effective_tree_size']
    assert rf_max == comp['max_rf']
    assert RF == comp['rf']
    # Let's insert some random nodes, that should be ignored
    for target in rng.sample([n for n in t2.descendants() if not n.is_leaf], 5):
        target.populate(5, dist_fn=rng.random, support_fn=rng.random)
    comp = t1.compare(t2, unrooted=unrooted)
    assert 20 == comp['effective_tree_size']
    assert rf_max == comp['max_rf']
    assert RF == comp['rf']


def test_monophyly():
    """Checks for monophyletic, paraphyletic, and polyphyletic groups."""
    t = Tree("((((((a, e), i), o),h), u), ((f, g), j));")
    #          ╭─┬╴a
    #        ╭─┤ ╰╴e
    #      ╭─┤ ╰╴i
    #    ╭─┤ ╰╴o
    #  ╭─┤ ╰╴h
    # ─┤ ╰╴u
    #  │ ╭─┬╴f
    #  ╰─┤ ╰╴g
    #    ╰╴j

    is_mono, monotype, extra  = t.check_monophyly(values=['a', 'e', 'i', 'o', 'u'])
    assert is_mono is False
    assert monotype == 'paraphyletic'

    is_mono, monotype, extra = t.check_monophyly(values=['a', 'e', 'i', 'o'])
    assert is_mono is True
    assert monotype == 'monophyletic'

    is_mono, monotype, extra =  t.check_monophyly(values=['i', 'o'])
    assert is_mono is False
    assert monotype == 'paraphyletic'

    # Now with unrooted trees, and using species instead of names.
    t = PhyloTree('(aaa1, (aaa3, (aaa4, (bbb1, bbb2))));',
                  sp_naming_function=lambda name: name[:3])

    assert t.check_monophyly(values={'aaa'},
                             prop='species', unrooted=True) == \
        (True, 'monophyletic', set())

    # Variations on that tree.
    t = PhyloTree('(aaa1, (bbb3, (aaa4, (bbb1, bbb2))));',
                  sp_naming_function=lambda name: name[:3])
    is_mono, _, extra = t.check_monophyly(values={'aaa'},
                                          prop='species', unrooted=True)
    assert not is_mono
    assert extra == {t['bbb3']}

    t = PhyloTree('(aaa1, (aaa3, (aaa4, (bbb1, bbb2))));',
                  sp_naming_function=lambda name: name[:3])
    is_mono, _, extra = t.check_monophyly(values={'bbb'},
                                          prop='species', unrooted=True)
    assert is_mono
    assert extra == set()

    t = PhyloTree('(aaa1, (aaa3, (aaa4, (bbb1, ccc2))));',
                  sp_naming_function=lambda name: name[:3])
    is_mono, _, extra = t.check_monophyly(values={'bbb', 'ccc'},
                                          prop='species', unrooted=True)
    assert is_mono
    assert extra == set()

    t = PhyloTree('(aaa1, (aaa3, (bbb4, (bbb1, bbb2))));',
                  sp_naming_function=lambda name: name[:3])
    is_mono, _, extra = t.check_monophyly(values={'bbb4', 'bbb2'},
                                          prop='name', unrooted=True)
    assert not is_mono
    assert extra == {t['bbb1']}

    t = PhyloTree('(aaa1, (aaa3, (bbb4, (bbb1, bbb2))));',
                  sp_naming_function=lambda name: name[:3])
    is_mono, _, extra = t.check_monophyly(values={'bbb1', 'bbb2'},
                                          prop='name', unrooted=True)
    assert is_mono
    assert extra == set()

    t = PhyloTree('(aaa1, aaa3, (aaa4, (bbb1, bbb2)));',
                  sp_naming_function=lambda name: name[:3])
    is_mono, _, extra = t.check_monophyly(values={'aaa'},
                                          prop='species', unrooted=True)
    assert is_mono
    assert extra == set()

    t = PhyloTree('(aaa1, bbb3, (aaa4, (bbb1, bbb2)));',
                  sp_naming_function=lambda name: name[:3])
    is_mono, _, extra = t.check_monophyly(values={'aaa'},
                                          prop='species', unrooted=True)
    assert not is_mono
    assert extra == {t['bbb3']}

    # Testing get_monophyly
    t = Tree('((((((4, e), i)M1, o),h), u), ((3, 4), (i, june))M2);',
             parser=1)
    # we annotate the tree using external data
    colors = {'a': 'red', 'e': 'green', 'i': 'yellow',
              'o': 'black', 'u': 'purple', '4': 'green',
              '3': 'yellow', '1': 'white', '5': 'red',
              'june': 'yellow'}
    for leaf in t:
        leaf.add_props(color=colors.get(leaf.name, 'none'))
    green_yellow_nodes = {t['M1'], t['M2']}
    mono_nodes = t.get_monophyletic(values=['green', 'yellow'],
                                    prop='color')
    assert set(mono_nodes) == green_yellow_nodes


@pytest.mark.parametrize("newick, expected_dists, expected_leaves", COPHENETIC_MATRIX_CASES)
def test_cophenetic_matrix(newick, expected_dists, expected_leaves):
    t = Tree(newick)
    dists, leaves = t.cophenetic_matrix()
    for i in range(len(expected_dists)):
        for j in range(len(expected_dists[i])):
            assert round(abs(expected_dists[i][j] - dists[i][j]), 4) == 0
    assert expected_leaves == leaves

