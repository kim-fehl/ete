import sys
import itertools

import pytest

from ete4 import Tree, PhyloTree
from ete4.core.tree import TreeError

from . import conftest as ds


def test_concat_trees():
    t1 = Tree('((A, B), C);')
    t2 = Tree('((a, b), c);')

    concat_tree = t1 + t2
    concat_tree.sort_descendants()
    assert concat_tree.write(parser=9) == '(((A,B),C),((a,b),c));'

    t3 = PhyloTree('((a, b), c);')

    mixed_types = lambda: t1 + t3
    with pytest.raises(TreeError):
        mixed_types()


def test_tree_manipulation():
    """Test operations which modify the tree topology."""
    nw_tree = "((Hi:1,Turtle:1.3)1:1,(A:0.3,B:2.4)1:0.43);"

    # Manipulate Topologies
    # Adding and removing nodes (add_child, remove_child,
    # add_sister, remove_sister). The resulting newick tree should
    # match the nw_tree defined before.
    t = Tree()

    remove_child_except = lambda: t.remove_child(t)
    add_sister_except = lambda: t.add_sister()
    with pytest.raises(TreeError):
        remove_child_except()
    with pytest.raises(TreeError):
        add_sister_except()

    c1 = t.add_child(dist=1, support=1)
    c2 = t.add_child(dist=0.43, support=1)
    n = Tree({'name': 'Hi', 'dist': 1, 'support': 1})
    _n = c1.add_child(n)
    c3 = _n.add_sister(name="Turtle", dist="1.3")
    c4 = c2.add_child(name="A", dist="0.3")

    c5 = c2.add_child(name="todelete")
    _c5 = c2.remove_child(c5)

    c6 = c2.add_child(name="todelete")
    _c6 = c4.remove_sister()

    c7 = c2.add_child(name="B", dist=2.4)

    assert nw_tree == t.write(parser=0)
    assert _c5 == c5
    assert _c6 == c6
    assert _n == n

    # Delete,
    t = Tree("(((A, B), C)[&&NHX:name=I], (D, F)[&&NHX:name=J])[&&NHX:name=root];")
    D = t['D']
    F = t['F']
    J = t['J']
    root = t['root']
    J.delete()
    assert J.up is None
    assert (J in t) is False
    assert D.up == root
    assert F.up == root

    # Delete preventing non dicotomic
    t = Tree('((((A:1,B:1):1,C:1):1,D:1):1,E:1):0;')
    orig_dist = t.get_distance(t, 'A')
    C = t['C']
    C.delete(preserve_branch_length=True)
    assert orig_dist == t.get_distance(t, 'A')

    t = Tree('((((A:1,B:1):1,C:1):1,D:1):1,E:1):0;')
    orig_dist = t.get_distance(t, 'A')
    C = t['C']
    C.delete(preserve_branch_length=False)
    assert orig_dist == t.get_distance(t, 'A') + 1

    t = Tree('((((A:1,B:1):1,C:1):1,D:1):1,E:1):0;')
    orig_dist = t.get_distance(t, 'A')
    C = t['C']
    C.delete(prevent_nondicotomic=False)
    assert orig_dist == t.get_distance(t, 'A')

    # detach
    t = Tree("(((A, B)[&&NHX:name=H], C)[&&NHX:name=I], (D, F)[&&NHX:name=J])[&&NHX:name=root];")
    D = t["D"]
    F = t["F"]
    J = t["J"]
    root = t["root"]
    J.detach()
    assert J.up is None
    assert (J in t) is False
    assert set([n.name for n in t.descendants()]) == set(["A", "B", "C", "I", "H"])

    # sorting branches
    t1 = Tree('((A,B),(C,D,E,F), (G,H,I));')
    t1.ladderize()
    assert list(t1.leaf_names()) == [_ for _ in 'ABGHICDEF']
    t1.ladderize(reverse=True)
    assert list(t1.leaf_names()) == [_ for _ in 'CDEFGHIAB']
    t1.sort_descendants()
    assert list(t1.leaf_names()) == [_ for _ in 'ABCDEFGHI']

    # prune
    t1 = Tree("(((A, B), C)[&&NHX:name=I], (D, F)[&&NHX:name=J])[&&NHX:name=root];")
    D1 = t1['D']
    t1.prune(["A", "C", D1])
    sys.stdout.flush()
    assert set([n.name for n in t1.descendants()]) == set(["A", "C", "D", "I"])

    t1 = Tree("(((A, B), C)[&&NHX:name=I], (D, F)[&&NHX:name=J])[&&NHX:name=root];")
    D1 = t1['D']
    t1.prune(["A", "B"])
    assert t1.write() == "(A,B);"

    # test prune keeping internal nodes

    t1 = Tree('(((((A,B)C)D,E)F,G)H,(I,J)K)root;', parser=1)
    t1.prune(['A', 'B', 'F', 'H'])
    assert set([n.name for n in t1.traverse()]) == set(['A', 'B', 'F', 'H', 'root'])

    t1 = Tree('(((((A,B)C)D,E)F,G)H,(I,J)K)root;', parser=1)
    t1.prune(['A', 'B'])
    assert set([n.name for n in t1.traverse()]) == set(['A', 'B', 'root'])

    t1 = Tree('(((((A,B)C)D,E)F,G)H,(I,J)K)root;', parser=1)
    t1.prune(['A', 'B', 'C'])
    assert set([n.name for n in t1.traverse()]) == set(['A', 'B', 'C', 'root'])

    t1 = Tree('(((((A,B)C)D,E)F,G)H,(I,J)K)root;', parser=1)
    t1.prune(['A', 'B', 'I'])
    assert set([n.name for n in t1.traverse()]) == set(['A', 'B', 'C', 'I', 'root'])


def test_remove_child(empty_tree, rng):
    """Test removing children."""
    t = empty_tree
    t.populate(20, dist_fn=rng.random, support_fn=rng.random)

    # Removed node loses its parent.
    node1 = t[1]
    assert node1.up == t  # no surprise here

    node1_returned = t.remove_child(t[1])
    assert node1_returned == node1
    assert node1.up is None  # node1 lost its parent

    t.add_child(node1)  # ok, let's put it back

    # A node that was already "stolen" does not lose its parent.
    node0 = t[0]  # first child from t, node that we will be moving around

    t2 = Tree()
    t2.add_child(node0)  # "steals" the first child from t

    assert t[0] == t2[0]  # they should be the same node

    t.remove_child(node0)
    assert t[0] != node0  # node is no longer a child of t

    assert node0.up == t2  # but has not lost track of its parent
    assert t2[0] == node0


def test_pop_child(empty_tree, rng):
    """Test popping children."""
    t = empty_tree
    t.populate(20, dist_fn=rng.random, support_fn=rng.random)

    # Removed node loses its parent.
    node1 = t[1]
    assert node1.up == t  # no surprise here

    node1_returned = t.pop_child()
    assert node1_returned == node1
    assert node1.up is None  # node1 lost its parent

    t.add_child(node1)  # ok, let's put it back

    # A node that was already "stolen" does not lose its parent.
    node0 = t[0]  # first child from t, node that we will be moving around

    t2 = Tree()
    t2.add_child(node0)  # "steals" the first child from t

    assert t[0] == t2[0]  # they should be the same node

    t.pop_child(0)
    assert t[0] != node0  # node is no longer a child of t

    assert node0.up == t2  # but has not lost track of its parent
    assert t2[0] == node0


def test_pruning(rng):
    # test prune preserving distances
    for i in range(3):  # NOTE: each iteration is quite slow
        t = Tree()
        t.populate(40, dist_fn=rng.random, support_fn=rng.random)
        orig_nw = t.write()
        distances = {}
        for a in t.leaves():
            for b in t.leaves():
                distances[(a, b)] = round(t.get_distance(a, b), 6)

        to_keep = set(rng.sample(list(t.leaves()), 6))
        t.prune(to_keep, preserve_branch_length=True)
        for a, b in distances:
            if a in to_keep and b in to_keep:
                assert distances[(a, b)] == round(t.get_distance(a, b), 6)

    # Total number of nodes is correct (no single child nodes)
    for x in range(10):
        t_fuzzy = Tree("(((A,B)1, C)2,(D,E)3)root;", parser=1)
        t_fuzzy.sort_descendants()
        orig_nw = t_fuzzy.write()
        ref_nodes = list(t_fuzzy.descendants())
        t_fuzzy.populate(10, dist_fn=rng.random, support_fn=rng.random)
        t_fuzzy['1'].populate(3, dist_fn=rng.random, support_fn=rng.random)
        t_fuzzy['2'].populate(5, dist_fn=rng.random, support_fn=rng.random)
        t_fuzzy['3'].populate(5, dist_fn=rng.random, support_fn=rng.random)
        t_fuzzy.prune(ref_nodes)
        t_fuzzy.sort_descendants()
        assert orig_nw == t_fuzzy.write()
        assert len(list(t_fuzzy.descendants())) == len(ref_nodes)

    # Total number of nodes is correct (no single child nodes)
    t = Tree()
    sample_size = 5
    t.populate(1000, dist_fn=rng.random, support_fn=rng.random)
    sample = rng.sample(list(t.leaves()), sample_size)
    t.prune(sample)
    assert len(t) == sample_size
    assert len(list(t.descendants())) == (sample_size * 2) - 2

    # Test preserve branch dist when pruning
    t = Tree()
    t.populate(100, dist_fn=rng.random, support_fn=rng.random)
    sample_size = 10  # NOTE: big values make this test very slow
    sample = rng.sample(list(t.leaves()), sample_size)
    matrix1 = ["%f" % t.get_distance(a, b) for (a, b) in itertools.product(sample, sample)]
    t.prune(sample, preserve_branch_length=True)
    matrix2 = ["%f" % t.get_distance(a, b) for (a, b) in itertools.product(sample, sample)]

    assert matrix1 == matrix2
    assert len(list(t.descendants())) == (sample_size * 2) - 2


def test_resolve_polytomy():
    t = Tree('((a,a,a,a),(b,b,b,(c,c,c)));')
    t.resolve_polytomy()
    assert t.write(parser=9) == '((((a,a),a),a),(((b,b),b),((c,c),c)));'

    t = Tree('((((a,a,a,a))),(b,b,b,(c,c,c)));')
    t.standardize()  # calls resolve_polytomy() internally too
    assert t.write(parser=9) == '((((b,b),b),((c,c),c)),(((a,a),a),a));'


def test_common_ancestors():
    # getting nodes, get_childs, get_sisters, root,
    # get_common_ancestor, get_nodes_by_name
    # get_descendants_by_name, is_leaf, is_root
    t = Tree("(((A,B)N1,C)N2[&&NHX:tag=common],D)[&&NHX:tag=root:name=root];", parser=1)
    assert t.get_sisters() == []

    A, B, C = t['A'], t['B'], t['C']
    root = t['root']
    assert "A" == A.name
    test_not_found = lambda: t['notfound']
    with pytest.raises(TreeError):
        test_not_found()

    assert "common" == t.common_ancestor([A, C]).props["tag"]
    assert "common" == t.common_ancestor([C, B]).props["tag"]
    assert root == t.common_ancestor([A, "D"])

    assert "root" == A.root.props["tag"]
    assert "root" == B.root.props["tag"]
    assert "root" == C.root.props["tag"]

    common = t.common_ancestor([C])
    assert "root" == common.root.props["tag"]

    assert common.root.is_root
    assert not A.is_root
    assert A.is_leaf
    assert not A.root.is_leaf
    with pytest.raises(TreeError):
        A.common_ancestor([Tree()])

    # Test multiple target nodes and lineage
    N1, N2 = t['N1'], t['N2']
    common = t.common_ancestor(['A', 'C'])
    assert common == N2

    expected_paths = {A: [A, N1, N2, root], C: [C, N2, root]}

    for node in expected_paths:
        assert list(node.lineage()) == expected_paths[node]

    # Test common ancestor function using self as single argument (issue #398)
    common = A.common_ancestor([A])
    assert common == A
    common = C.common_ancestor(["C"])
    assert common == C


def test_getters_iters(rng):
    # Iter ancestors
    t = Tree("(((((a,b)A,c)B,d)C,e)D,f)root;", parser=1)
    ancestor_names = [n.name for n in (t["a"]).ancestors()]
    assert ancestor_names == ["A", "B", "C", "D", "root"]
    ancestor_names = [n.name for n in (t["B"]).ancestors()]
    assert ancestor_names == ["C", "D", "root"]

    # Tree magic python features
    t = Tree(ds.nw_dflt)
    assert len(t) == 20
    assert "Ddi0002240" in t
    assert t.children[0] in t
    for a in t:
        assert a.name

    # Populate
    t = Tree(ds.nw_full)
    prev_size = len(t)
    t.populate(25, dist_fn=rng.random, support_fn=rng.random)
    assert len(t) == prev_size + 25
    for i in range(10):
        t = Tree()
        t.populate(100, dist_fn=rng.random, support_fn=rng.random)
        # Checks that all names are actually unique
        assert len(set(t.leaf_names())) == 100

    # Adding and removing features
    t = Tree("(((A,B),C)[&&NHX:tag=common],D)[&&NHX:tag=root];")
    A = t['A']

    # Check gettters and itters return the same
    t = Tree(ds.nw2_full)

    assert set([n for n in t.traverse("preorder")]) == set([n for n in t.traverse("postorder")])
    assert t in set([n for n in t.traverse("preorder")])

    # Check order of visiting nodes

    t = Tree("((3,4)2,(6,7)5)1;", parser=1)
    postorder = "3426751"
    preorder = "1234567"
    levelorder = "1253467"

    assert preorder == ''.join(n.name for n in t.traverse("preorder"))

    assert postorder == ''.join(n.name for n in t.traverse("postorder"))

    assert levelorder == ''.join(n.name for n in t.traverse("levelorder"))

    # Swap children.
    n = t.get_children()
    t.reverse_children()
    n.reverse()
    assert n == t.get_children()


def test_tree_navigation(rng):
    t = Tree('(((A,B)H,C)I,(D,F)J)root;', parser=1)
    postorder = [n.name for n in t.traverse("postorder")]
    preorder = [n.name for n in t.traverse("preorder")]
    levelorder = [n.name for n in t.traverse("levelorder")]

    assert postorder == ['A', 'B', 'H', 'C', 'I', 'D', 'F', 'J', 'root']
    assert preorder == ['root', 'I', 'H', 'A', 'B', 'C', 'J', 'D', 'F']
    assert levelorder == ['root', 'I', 'J', 'H', 'C', 'D', 'F', 'A', 'B']
    ancestors = [n.name for n in (t['B']).ancestors()]
    assert ancestors == ['H', 'I', 'root']
    assert list(t.ancestors()) == []

    # add something of is_leaf_fn etc...
    custom_test = lambda x: x.name in 'JCH'
    custom_leaves = t.leaves(is_leaf_fn=custom_test)
    assert {n.name for n in custom_leaves} == {'J', 'H', 'C'}

    # Test cached content
    t = Tree()
    t.populate(20, dist_fn=rng.random, support_fn=rng.random)

    cache_node = t.get_cached_content()
    cache_node_leaves_only_false = t.get_cached_content(leaves_only=False)
    assert cache_node[t] == set(t.leaves())
    assert cache_node_leaves_only_false[t] == set(t.traverse())

    cache_name = t.get_cached_content('name')
    cache_name_leaves_only_false = t.get_cached_content('name', leaves_only=False)
    assert cache_name[t] == set(t.leaf_names())
    assert cache_name_leaves_only_false[t] == set(n.name for n in t.traverse())

    cache_many = {n: {(ni.name, ni.dist, ni.support) for ni in nodes}
                  for n, nodes in t.get_cached_content().items()}
    cache_many_lof = {n: {(ni.name, ni.dist, ni.support) for ni in nodes}
                      for n, nodes in t.get_cached_content(leaves_only=False).items()}
    assert cache_many[t] == set([(leaf.name, leaf.dist, leaf.support) for leaf in t.leaves()])
    assert cache_many_lof[t] == set((n.name, n.dist, n.support) for n in t.traverse())

