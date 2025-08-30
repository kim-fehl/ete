import pytest

from ete4 import Tree

from . import conftest as ds


def strip(text):
    """Return the given text stripping the empty lines and indentation."""
    indent = min(len(line) - len(line.lstrip())
                 for line in text.splitlines() if line.strip())
    return '\n'.join(line[indent:].rstrip()
                     for line in text.splitlines() if line.strip())


def assert_looks(tree, text):
    """Assert that tree looks like the given text (as ascii)."""
    assert tree.to_str(compact=True, props=['name']) == strip(text)


def test_rooting_topology():
    """Test topology changes after rooting"""
    t = Tree('((d,e)b,(f,g)c);', parser=1)
    assert_looks(t, """
               ╭╴b╶┬╴d
            ╴⊗╶┤   ╰╴e
               ╰╴c╶┬╴f
                   ╰╴g
    """)

    t.set_outgroup(t['e'])
    assert_looks(t, """
            ╴⊗╶┬╴e
               ╰╴b╶┬╴d
                   ╰╴c╶┬╴f
                       ╰╴g
    """)

    t.set_outgroup(t['d'])
    assert_looks(t, """
               ╭╴d
            ╴⊗╶┤   ╭╴c╶┬╴f
               ╰╴b╶┤   ╰╴g
                   ╰╴e
    """)

    t.set_outgroup(t['c'])
    assert_looks(t, """
               ╭╴c╶┬╴f
            ╴⊗╶┤   ╰╴g
               ╰╴b╶┬╴e
                   ╰╴d
    """)

    t.set_outgroup(t['b'])
    assert_looks(t, """
               ╭╴b╶┬╴e
            ╴⊗╶┤   ╰╴d
               ╰╴c╶┬╴f
                   ╰╴g
    """)


def test_rooting_distances(rng):
    """Tests that set_outgroup operations never change distance relationships between nodes"""

    # Loads a large tree with different node distance
    t = Tree(ds.nw2_full)

    # records the distance between two random nodes
    YGR028W = t['YGR028W']
    YGR138C = t['YGR138C']
    d1 = t.get_distance(YGR138C, YGR028W)

    # records original sum up distances
    sum_distances = sum([n.dist for n in t.traverse() if n.dist])

    # gets midpoint outgroup and re-root the tree there
    midpoint = t.get_midpoint_outgroup()
    t.set_outgroup(midpoint)

    # saves the two nodes that split the tree by midpoint
    o1, o2 = t.children[0], t.children[1]

    # Test node distances is preserved
    d2 = t.get_distance(YGR138C, YGR028W)
    assert d1 == d2

    # Test sum up distance is intact after rooting
    assert sum_distances == sum(n.dist for n in t.traverse() if n.dist)

    # Let's now test whether can we recover the original state of the tree
    # after many random re-rooting operations

    nodes = list(t.descendants())

    for i in range(100):
        for j in range(100):
            # root at a random place
            n = rng.sample(nodes, 1)[0]
            t.set_outgroup(n)

        # Restore original midpoint outgroup. If everything was ok, sum up
        # distance and the same two root branches should be restored
        midpoint = t.get_midpoint_outgroup()
        t.set_outgroup(midpoint)
        assert set([t.children[0], t.children[1]]) == set([o1, o2])

        d3 = t.get_distance(t["YGR138C"], t["YGR028W"])
        assert d1 == d3

        # Test sum up distance is intact after rooting
        assert abs(
            sum_distances -
            sum(n.dist for n in t.traverse() if n.dist)) < 1e-8

    # Test that the distance of the two root branches after
    # rooting are balanced.
    t = Tree('((A:10,B:1),(C:1,D:1)E:1)root;', parser=1)
    t.set_outgroup(t.get_midpoint_outgroup())
    assert t.children[0].dist == 5.0
    assert t.children[1].dist == 5.0

    # Test that set_outgroup can root an "unrooted" tree (that is, a tree
    # whose root has more than 2 children).
    t = Tree('(A:10,B:1,(C:1,D:1)E:1)root;', parser=1)
    assert t.children[0] == t['A']
    t.set_outgroup(t['A'])

    # Test that the distance of the two root branches
    # after rooting are balanced even for unrooted trees
    assert t.children[0].dist == 5.0
    assert t.children[1].dist == 5.0


def test_unroot():
    # Simple case. We start with a dicotomy from the root.
    t = Tree('((a:0.5,b:0.5):0.5,(c:0.2,d:0.2):0.8);')
    t.unroot()
    assert '(a:0.5,b:0.5,(c:0.2,d:0.2):1.3);' == t.write()

    # If we unroot an unrooted tree, it should stay the same:
    t.unroot()
    assert '(a:0.5,b:0.5,(c:0.2,d:0.2):1.3);' == t.write()

    # Test unrooting when we have more branch properties.
    t = Tree('((a:0.5[&&NHX:color=green],b:0.5[&&NHX:color=red]):0.5'
             '[&&NHX:color=green],(c:0.2[&&NHX:color=green],d:0.2'
             '[&&NHX:color=blue]):0.8[&&NHX:color=green]);')
    t.unroot(bprops=['color'])
    assert t.write(props=None) == \
        '(a:0.5[&&NHX:color=green],b:0.5[&&NHX:color=red],(c:0.2' \
        '[&&NHX:color=green],d:0.2[&&NHX:color=blue]):1.3[&&NHX:color=green]);'

    # If the branch properties are not consistent, we have an exception.
    t = Tree('((a:0.5[&&NHX:color=green],b:0.5[&&NHX:color=red]):0.5'
             '[&&NHX:color=green],(c:0.2[&&NHX:color=green],d:0.2'
             '[&&NHX:color=blue]):0.8[&&NHX:color=red]);')
    with pytest.raises(AssertionError):
        t.unroot(bprops=['color'])


def test_rooting_branch_support(rng):
    """Test that branch support, distances and custom branch properties are correctly handled after re-rooting."""

    # Generate a random tree. Test branch support and distances after rooting.
    t = Tree()
    t.populate(50, dist_fn=rng.random, support_fn=lambda: rng.uniform(0, 100))
    t.unroot()

    # Add a branch property.
    rand_value = rng.random()
    for ch in t.children:
        ch.props['bprop'] = rand_value

    for n in t.descendants():
        if n.up is not t:
            n.props['bprop'] = rng.random()

    # Record the distance and support value of all clades, based on its content
    names = set(t.leaf_names())
    cluster_id2support = {}
    cluster_id2dist = {}
    cluster_id2bprop = {}
    for n in t.descendants():
        cluster_names = set(n.leaf_names())
        cluster_names2 = names - cluster_names
        cluster_id = '_'.join(sorted(cluster_names))
        cluster_id2 = '_'.join(sorted(cluster_names2))
        cluster_id2support[cluster_id] = n.support
        cluster_id2support[cluster_id2] = n.support

        cluster_id2dist[cluster_id] = n.dist
        cluster_id2dist[cluster_id2] = n.dist
        cluster_id2bprop[cluster_id] = n.props["bprop"]
        cluster_id2bprop[cluster_id2] = n.props["bprop"]

    # Root to every single node in the tree and test whether all partitions conserve their properties
    for outgroup in t.descendants():
        t.set_outgroup(outgroup, bprops=["bprop"])

        for n in t.descendants():
            cluster_names = set(n.leaf_names())
            cluster_names2 = names - cluster_names
            cluster_id = '_'.join(sorted(cluster_names))
            cluster_id2 = '_'.join(sorted(cluster_names2))

            assert cluster_id2support.get(cluster_id, None) == n.support
            assert cluster_id2support.get(cluster_id2, None) == n.support

            assert cluster_id2bprop.get(cluster_id, None) == n.props.get("bprop")
            assert cluster_id2bprop.get(cluster_id2, None) == n.props.get("bprop")

            if n.up and n.up.up:
                assert cluster_id2dist.get(cluster_id, None) == n.dist

