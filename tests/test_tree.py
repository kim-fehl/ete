import sys
import random
import itertools
import json
from tempfile import NamedTemporaryFile

from ete4 import Tree, PhyloTree
from ete4.core.tree import TreeError
from ete4.parser.newick import NewickError
from ete4.parser import newick

from . import conftest as ds
from .conftest import CUSTOM_FORMAT_CASES, QUOTED_NAME_CASES, COMPLEX_NAME, COPHENETIC_MATRIX_CASES, ROBINSON_FOULDS_CASES
import pytest

@pytest.fixture
def empty_tree():
    return Tree()


@pytest.fixture
def rng():
    """Deterministic random generator for reproducible tests."""
    return random.Random(0)


def strip(text):
    """Return the given text stripping the empty lines and indentation."""
    # Helps compare tree visualizations.
    indent = min(len(line) - len(line.lstrip())
                 for line in text.splitlines() if line.strip())
    return '\n'.join(line[indent:].rstrip()
        for line in text.splitlines() if line.strip())


def assert_looks(tree, text):
    """Assert that tree looks like the given text (as ascii)."""
    assert tree.to_str(compact=True, props=['name']) == \
                     strip(text)

def test_read_write_exceptions():
    """Test that the right exceptions are risen."""
    def wrong_dist():
        t = Tree()
        t.dist = '1a'

    def wrong_support():
        t = Tree()
        t.support = '1a'

    def wrong_up():
        t = Tree()
        t.up = 'Something'

    def wrong_children():
        t = Tree()
        t.children = 'Something'

    with pytest.raises(ValueError):
        wrong_dist()
    with pytest.raises(ValueError):
        wrong_support()
    with pytest.raises(TypeError):
        wrong_up()
    with pytest.raises(TreeError):
        wrong_children()

def test_add_remove_properties(empty_tree):
    t = empty_tree
    t.add_props(testf1=1, testf2="1", testf3=[1])
    t.add_prop('testf4', set([1]))
    assert t.props["testf1"] == 1
    assert t.props["testf2"] == "1"
    assert t.props["testf3"] == [1]
    assert t.props["testf4"] == set([1])

    t.del_prop('testf4')
    assert 'testf4' not in t.props

def test_basic_properties():
    t = Tree('((a:1,b:2)0.8:3,c:4);')

    assert type(t['a'].dist) == float
    assert type(t[0].support) == float

    t[0].support = None
    assert 'support' not in t[0].props
    assert t[0].support is None

    t[0].dist = None
    assert 'dist' not in t[0].props
    assert t[0].dist is None

    a = t['a']
    a.name = None
    assert 'name' not in a.props
    assert a.name is None

    a.name = 5
    assert a.name == '5'

def test_tree_read_and_write(rng):
    """Test newick support."""
    # Read and write newick tree from/to file.
    with NamedTemporaryFile() as f_tree:  # test reading from file
        f_tree.write(ds.nw_full.encode('utf8'))
        f_tree.flush()
        t = Tree(open(f_tree.name))

    assert ds.nw_full == t.write(props=["flag", "mood"])
    assert ds.nw_topo == t.write(parser=9)
    assert ds.nw_dist == t.write(parser=5)

    with NamedTemporaryFile() as f_writetest:  # test writing to file
        t.write(outfile=f_writetest.name)
        assert Tree(open(f_writetest.name)).write() == t.write()

    # Read and write newick tree from/to string.
    t = Tree(ds.nw_full)

    assert ds.nw_full == t.write(props=["flag", "mood"])
    assert ds.nw_topo == t.write(parser=9)
    assert ds.nw_dist == t.write(parser=5)

    # Read complex newick.
    t = Tree(ds.nw2_full)
    assert ds.nw2_full == t.write()

    # Read weird topologies.
    t = Tree(ds.nw_simple5)
    assert ds.nw_simple5 == t.write(parser=9)

    t = Tree(ds.nw_simple6)
    assert ds.nw_simple6 == t.write(parser=9)

    # Read single node trees.
    assert Tree("hello;").write(parser=9, format_root_node=True) == "hello;"
    assert Tree("(hello);").write(parser=9) == "(hello);"

    # Export root features.
    newick = "(((A[&&NHX:name=A],B[&&NHX:name=B])[&&NHX:name=NoName],C[&&NHX:name=C])[&&NHX:name=I],(D[&&NHX:name=D],F[&&NHX:name=F])[&&NHX:name=J])[&&NHX:name=root];"
    t = Tree(newick)
    assert t.write(parser=9, props=['name'], format_root_node=True) == '(((A,B)[&&NHX:name=NoName],C)[&&NHX:name=I],(D,F)[&&NHX:name=J])[&&NHX:name=root];'

    # Export ordered features.
    t = Tree("((A,B),C);")
    expected_nw = "((A,B[&&NHX:0=0:1=1:2=2:3=3:4=4:5=5:6=6:7=7:8=8:9=9:a=a:b=b:c=c:d=d:e=e:f=f:g=g:h=h:i=i:j=j:k=k:l=l:m=m:n=n:o=o:p=p:q=q:r=r:s=s:t=t:u=u:v=v:w=w]),C);"
    features = list("abcdefghijklmnopqrstuvw0123456789")
    rng.shuffle(features)
    for letter in features:
        t['B'].add_prop(letter, letter)
    assert expected_nw == t.write(props=None)

def test_repr(empty_tree):
    """Test that the Tree representation looks like we expect."""
    t = empty_tree
    r1, r2, r3 = t.__repr__(), repr(t), '%r' % t
    assert r1 == r2 == r3
    assert r1.startswith('<Tree ')
    assert r1.endswith('>')
    assert ' at 0x' in r1

def test_to_str():
    """Test that the ascii representation (to use when printing) works."""
    t = Tree('((a,b)x,(c,d)y);', parser=1)

    assert t.to_str(props=['name']) == strip("""
           ╭╴a
       ╭╴x╶┤
       │   ╰╴b
    ╴⊗╶┤
       │   ╭╴c
       ╰╴y╶┤
           ╰╴d
    """)
    assert t.to_str(compact=True, show_internal=False, props=['name']) == strip("""
     ╭─┬╴a
    ─┤ ╰╴b
     ╰─┬╴c
       ╰╴d
    """)
    assert t.to_str(props=['name'], cascade=True) == strip("""
    ⊗
    ├─┐x
    │ ├─╴a
    │ └─╴b
    └─┐y
      ├─╴c
      └─╴d
    """)
    assert t.to_str(props=['dist', 'support']) == strip("""
               ╭╴⊗,⊗
         ╭╴⊗,⊗╶┤
         │     ╰╴⊗,⊗
    ╴⊗,⊗╶┤
         │     ╭╴⊗,⊗
         ╰╴⊗,⊗╶┤
               ╰╴⊗,⊗
    """)

    t2 = Tree('((a:1,b:2)x:3,(c:4,d:5)y:6);', parser=1)
    assert t2.to_str(props=['dist']) == strip("""
             ╭╴1.0
       ╭╴3.0╶┤
       │     ╰╴2.0
    ╴⊗╶┤
       │     ╭╴4.0
       ╰╴6.0╶┤
             ╰╴5.0
    """)

    assert t2.to_str() == strip("""
                        ╭╴name=a,dist=1.0
      ╭╴name=x,dist=3.0╶┤
      │                 ╰╴name=b,dist=2.0
    ──┤
      │                 ╭╴name=c,dist=4.0
      ╰╴name=y,dist=6.0╶┤
                        ╰╴name=d,dist=5.0
    """)

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

def test_newick_formats(rng):
    """Test different newick subformats."""
    NW_FORMAT = [0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 100]  # predefined parsers ("formats")

    # Let's stress a bit
    for i in range(10):
        t = Tree()
        t.populate(4, dist_fn=rng.random, support_fn=rng.random)
        for n in t.traverse():
            n.name = n.name or 'NoName'
            n.support = n.support or 1
        for f in NW_FORMAT:
            assert t.write(parser=f) == Tree(t.write(parser=f),parser=f).write(parser=f)

    # Format 0 = ((H:1,(G:1,F:1)1:1)1:1,I:1)1:1;
    # Format 1 = ((H:1,(G:1,F:1):1):1,I:1):1;
    # Format 2 = ((H:1,(G:1,F:1)1:1)1:1,I:1)1:1;
    # Format 3 = ((H:1,(G:1,F:1)NoName:1)NoName:1,I:1)NoName:1;
    # Format 4 = ((H:1,(G:1,F:1)),I:1);
    # Format 5 = ((H:1,(G:1,F:1):1):1,I:1):1;
    # Format 6 = ((H,(G,F):1):1,I):1;
    # Format 7 = ((H:1,(G:1,F:1)NoName)NoName,I:1)NoName;
    # Format 8 = ((H,(G,F)NoName)NoName,I)NoName;
    # Format 9 = ((H,(G,F)),I);
    # Format 100 = ((,(,)),);

    t = Tree()
    t.populate(50, dist_fn=rng.random, support_fn=rng.random)
    for n in t.traverse():
        n.name = n.name or 'NoName'
        n.support = n.support if n.support is not None else 1
    t.sort_descendants()
    expected_distances = [round(n.dist, 6) for n in t.traverse('postorder') if n.up]
    expected_leaf_distances = [round(n.dist, 6) for n in t]
    expected_internal_distances = [round(n.dist, 6) for n in t.traverse('postorder') if not n.is_leaf and n.up]
    expected_supports = [round(n.support, 6) for n in t.traverse('postorder') if not n.is_leaf and n.up]
    expected_leaf_names = [n.name for n in t]

    # Check that all formats read names correctly
    for f in [0,1,2,3,5,6,7,8,9]:
        t2 = Tree(t.write(parser=f), parser=f)
        t2.sort_descendants()
        observed_names = [n.name for n in t]
        assert observed_names == expected_leaf_names

    # Check that all formats reading distances, recover original distances
    for f in [0,1,2,3,5]:
        t2 = Tree(t.write(parser=f), parser=f)
        t2.dist, t2.support = 0, 1
        t2.sort_descendants()
        observed_distances = [n.dist for n in t2.traverse('postorder')]
        assert all(abs(x - y) < 1e-6 for x, y in zip(observed_distances, expected_distances))

    # formats reading only leaf distances
    for f in [4,7]:
        t2 = Tree(t.write(parser=f), parser=f)
        t2.dist, t2.support = 0, 1
        t2.sort_descendants()
        observed_distances = [n.dist for n in t2]
        assert all(abs(x - y) < 1e-6 for x, y in zip(observed_distances, expected_leaf_distances))

    # formats reading only leaf distances
    for f in [6]:
        t2 = Tree(t.write(parser=f), parser=f)
        t2.dist, t2.support = 0, 1
        t2.sort_descendants()
        observed_distances = [n.dist for n in t2.traverse('postorder') if not n.is_leaf]
        assert all(abs(x - y) < 1e-6 for x, y in zip(observed_distances, expected_internal_distances))


    # Check that all formats reading supports, recover original distances
    for f in [0,2]:
        t2 = Tree(t.write(parser=f), parser=f)
        t2.dist, t2.support = 0, 1
        t2.sort_descendants()
        observed_supports = [n.support for n in t2.traverse('postorder') if not n.is_leaf]
        assert all(abs(x - y) < 1e-6 for x, y in zip(observed_supports, expected_supports))


   # Check that formats reading supports, do not accept node names
    for f in [0,2]:
        # format 3 forces dumping internal node names, NoName in case is missing
        with pytest.raises(Exception):
            Tree(t.write(parser=3), parser=f)

   # Check that formats reading names, do not load supports
    for f in [1, 3]:
        t2 = Tree(t.write(parser=0), parser=f)
        default_supports = set([n.support for n in t2.traverse()])
        assert {None} == default_supports


    # Check errors reading numbers
    error_nw1 = "((A:0.813705,(E:0.545591,D:0.411772)error:0.137245)1.000000:0.976306,C:0.074268);"
    for f in [0, 2]:
        with pytest.raises(NewickError):
            Tree(error_nw1, parser=f)

    error_nw2 = "((A:0.813705,(E:0.545error,D:0.411772)1.0:0.137245)1.000000:0.976306,C:0.074268);"
    for f in [0, 1, 2]:
        with pytest.raises(NewickError):
            Tree(error_nw2, parser=f)


    error_nw3 = "((A:0.813705,(E:0.545error,D:0.411772)1.0:0.137245)1.000000:0.976306,C:0.074268);"
    for f in [0, 1, 2]:
        with pytest.raises(NewickError):
            Tree(error_nw2, parser=f)

    # Check errors derived from reading names with weird or illegal chars
    base_nw = "((NAME1:0.813705,(NAME2:0.545,NAME3:0.411772)NAME6:0.137245)NAME5:0.976306,NAME4:0.074268);"
#        valid_names = ['[name]', '[name', '"name"', "'name'", "'name", 'name', '[]\'"&%$!*.']
    valid_names = ['name']
#        error_names = ['error)', '(error', "erro()r",  ":error", "error:", "err:or", ",error",  "error,"]
    error_names = ['error)', "erro()r",  ":error", "error:", "err:or"]
    for ename in error_names:
        with pytest.raises(NewickError):
            Tree(base_nw.replace('NAME2', ename), parser=1)
        if not ename.startswith(','):
            with pytest.raises(NewickError):
                Tree(base_nw.replace('NAME6', ename), parser=1)

    for vname in valid_names:
        expected_names = set(['NAME1', vname, 'NAME3', 'NAME4'])
        assert set([n.name for n in Tree(base_nw.replace('NAME2', vname), parser=1)]) == \
                         expected_names

    # invalid NHX format
    with pytest.raises(NewickError):
        Tree("(((A, B), C)[&&NHX:nameI]);")
    # unsupported newick stream
    with pytest.raises(Exception):
        Tree([1,2,3])

def test_newick_multisupport():
    nw = '((a,b)2/3:4,(c,d)5/6:7);'
    t = Tree(nw, parser='multisupport')
    assert t.write(parser='multisupport') == nw

@pytest.mark.parametrize("nw, nw_normalized", QUOTED_NAME_CASES)
def test_quoted_names(nw, nw_normalized):
    with pytest.raises(NewickError):
        Tree(nw, parser=0)
    t = Tree(nw, parser=1)
    assert any(n for n in t if n.name == COMPLEX_NAME)
    # test writing and reloading tree
    nw_back = t.write(parser=1)
    t2 = Tree(nw, parser=1)
    nw_back2 = t2.write(parser=1)
    assert nw_normalized == nw_back
    assert nw_normalized == nw_back2

@pytest.mark.parametrize("flag, result", CUSTOM_FORMAT_CASES)
def test_custom_formatting_formats(flag, result):
    """Test change dist, name and support formatters."""
    t = Tree('((A:1.1111,B:2.2222)C:3.3333[&&NHX:support=1],D:4.4444);',
             parser=1)
    t.sort_descendants()

    parser = newick.make_parser(flag, dist='%0.1f', name='TEST-%s',
                                support='SUP-%0.1f')
    nw = t.write(parser=parser)
    assert nw == result

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
    assert J.up == None
    assert (J in t) == False
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
    assert orig_dist == t.get_distance(t, 'A')+1

    t = Tree('((((A:1,B:1):1,C:1):1,D:1):1,E:1):0;')
    orig_dist = t.get_distance(t, 'A')
    C = t['C']
    C.delete(prevent_nondicotomic=False)
    assert orig_dist == t.get_distance(t, 'A')

    #detach
    t = Tree("(((A, B)[&&NHX:name=H], C)[&&NHX:name=I], (D, F)[&&NHX:name=J])[&&NHX:name=root];")
    D = t["D"]
    F = t["F"]
    J = t["J"]
    root = t["root"]
    J.detach()
    assert J.up == None
    assert (J in t) == False
    assert set([n.name for n in t.descendants()]) == set(["A","B","C","I","H"])

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
    t1.prune(["A","C", D1])
    sys.stdout.flush()
    assert set([n.name for n in t1.descendants()]) == set(["A","C","D","I"])

    t1 = Tree("(((A, B), C)[&&NHX:name=I], (D, F)[&&NHX:name=J])[&&NHX:name=root];")
    D1 = t1['D']
    t1.prune(["A","B"])
    assert t1.write() == "(A,B);"

    # test prune keeping internal nodes

    t1 = Tree('(((((A,B)C)D,E)F,G)H,(I,J)K)root;', parser=1)
    t1.prune(['A', 'B', 'F', 'H'])
    assert set([n.name for n in t1.traverse()]) == \
                     set(['A', 'B', 'F', 'H', 'root'])

    t1 = Tree('(((((A,B)C)D,E)F,G)H,(I,J)K)root;', parser=1)
    t1.prune(['A', 'B'])
    assert set([n.name for n in t1.traverse()]) == \
                     set(['A', 'B', 'root'])

    t1 = Tree('(((((A,B)C)D,E)F,G)H,(I,J)K)root;', parser=1)
    t1.prune(['A', 'B', 'C'])
    assert set([n.name for n in t1.traverse()]) == \
                     set(['A', 'B', 'C', 'root'])

    t1 = Tree('(((((A,B)C)D,E)F,G)H,(I,J)K)root;', parser=1)
    t1.prune(['A', 'B', 'I'])
    assert set([n.name for n in t1.traverse()]) == \
                     set(['A', 'B', 'C', 'I', 'root'])

def test_remove_child(empty_tree, rng):
    """Test removing children."""
    t = empty_tree
    t.populate(20, dist_fn=rng.random, support_fn=rng.random)

    # Removed node loses its parent.
    node1 = t[1]
    assert node1.up == t  # no surprise here

    node1_returned = t.remove_child(t[1])
    assert node1_returned == node1
    assert node1.up == None  # node1 lost its parent

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
    assert node1.up == None  # node1 lost its parent

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
                distances[(a,b)] = round(t.get_distance(a, b), 6)

        to_keep = set(rng.sample(list(t.leaves()), 6))
        t.prune(to_keep, preserve_branch_length=True)
        for a,b in distances:
            if a in to_keep and b in to_keep:
                assert distances[(a,b)] == round(t.get_distance(a, b), 6)

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
    assert len(list(t.descendants())) == (sample_size*2)-2

    # Test preserve branch dist when pruning
    t = Tree()
    t.populate(100, dist_fn=rng.random, support_fn=rng.random)
    sample_size = 10  # NOTE: big values make this test very slow
    sample = rng.sample(list(t.leaves()), sample_size)
    matrix1 = ["%f" % t.get_distance(a, b) for (a,b) in itertools.product(sample, sample)]
    t.prune(sample, preserve_branch_length=True)
    matrix2 = ["%f" % t.get_distance(a, b) for (a,b) in itertools.product(sample, sample)]

    assert matrix1 == matrix2
    assert len(list(t.descendants())) == (sample_size*2)-2

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
    t = Tree("(((A,B)N1,C)N2[&&NHX:tag=common],D)[&&NHX:tag=root:name=root];",
             parser=1)
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
    prev_size= len(t)
    t.populate(25, dist_fn=rng.random, support_fn=rng.random)
    assert len(t) == prev_size+25
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

    assert set([n for n in t.traverse("preorder")]) == \
                         set([n for n in t.traverse("postorder")])
    assert t in set([n for n in t.traverse("preorder")])

    # Check order or visiting nodes

    t = Tree("((3,4)2,(6,7)5)1;", parser=1)
    postorder = "3426751"
    preorder = "1234567"
    levelorder = "1253467"

    assert preorder == \
                     ''.join(n.name for n in t.traverse("preorder"))

    assert postorder == \
                     ''.join(n.name for n in t.traverse("postorder"))

    assert levelorder == \
                     ''.join(n.name for n in t.traverse("levelorder"))

    # Swap children.
    n = t.get_children()
    t.reverse_children()
    n.reverse()
    assert n == t.get_children()

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
    assert root.get_farthest_leaf() == (A,1.101)
    assert root.get_farthest_node() == (A,1.101)
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
    """ Tests that set_outgroup operations never changes distance relationships between nodes"""

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
    assert sum_distances == \
                     sum(n.dist for n in t.traverse() if n.dist)

    # Let's now test whether can we recover the original state of the tree
    # after many random re-rooting operations

    nodes = list(t.descendants())

    for i in range(100):
        for j in range(100):
            # root at a random place
            n = rng.sample(nodes, 1)[0]
            t.set_outgroup(n)

        # Restore original midpoint outgroup. If everything was ok, sum up
        # distance and the same two root branches should be restoredœ
        midpoint = t.get_midpoint_outgroup()
        t.set_outgroup(midpoint)
        assert set([t.children[0], t.children[1]]) == set([o1, o2])

        d3 = t.get_distance(t["YGR138C"], t["YGR028W"])
        assert d1 == d3

        # Test sum up distance is intact after rooting
        assert abs(
            sum_distances -
            sum(n.dist for n in t.traverse() if n.dist)) \
                        < 1e-8

    # Test that the distance of the two root branches after
    # rooting are balanced.
    t = Tree('((A:10,B:1),(C:1,D:1)E:1)root;', parser=1);
    t.set_outgroup(t.get_midpoint_outgroup())
    assert t.children[0].dist == 5.0
    assert t.children[1].dist == 5.0

    # Test that set_outgroup can root an "unrooted" tree (that is, a tree
    # whose root has more than 2 children).
    t = Tree('(A:10,B:1,(C:1,D:1)E:1)root;', parser=1);
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


    #self.assertEqual(cache_name_lof[t], [t.name])

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

def test_describe():
    assert Tree().describe() == \
                     'Number of leaf nodes: 1\n' \
                     'Total number of nodes: 1\n' \
                     'Rooted: No children\n' \
                     'Most distant node: \n' \
                     'Max. distance: 0'
    assert Tree('(a,b,c);').describe() == \
                     'Number of leaf nodes: 3\n' \
                     'Total number of nodes: 4\n' \
                     'Rooted: No\n' \
                     'Most distant node: a\n' \
                     'Max. distance: 1'
    assert Tree('(a,(b,c));').describe() == \
                     'Number of leaf nodes: 3\n' \
                     'Total number of nodes: 5\n' \
                     'Rooted: Yes\n' \
                     'Most distant node: b\n' \
                     'Max. distance: 2'

def test_treeid(rng):
    t = Tree()
    t.populate(50, dist_fn=rng.random, support_fn=rng.random)
    orig_id = t.get_topology_id()
    nodes = list(t.descendants())
    for i in range(20):
        for n in rng.sample(nodes, 10):
            n.reverse_children()
            assert t.get_topology_id() == orig_id

def test_node_id():
    """Test the node_id corresponding to a node inside a tree."""
    t = Tree('((a,b)x,(c,d)y);', parser=1)
    #     ╭╴x (0,)╶┬╴a (0,0)
    # ╴()╶┤        ╰╴b (0,1)
    #     ╰╴y (1,)╶┬╴c (1,0)
    #              ╰╴d (1,1)

    assert t.id == ()
    assert t['x'].id == (0,)
    assert t['y'].id == (1,)
    assert t['a'].id == (0,0)
    assert t['d'].id == (1,1)

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

# TODO: Fix the compare() function and this test.
# def test_tree_compare():
#     def _astuple(d):
#         keynames = ["norm_rf", "rf", "max_rf", "ref_edges_in_source",
#                     "source_edges_in_ref", "effective_tree_size",
#                     "source_subtrees", "treeko_dist"]
#         return tuple([d[v] for v in keynames])

#     ref1 = Tree('((((A, B)0.91, (C, D))0.9, (E,F)0.96), (G, H));', parser=0)
#     ref2 = Tree('(((A, B)0.91, (C, D))0.9, (E,F)0.96);', parser=0)
#     s1 = Tree('(((A, B)0.9, (C, D))0.9, (E,F)0.9);', parser=0)

#     small = Tree("((A, B), C);")
#     # RF unrooted in too small trees for rf, but with at least one internal node
#     self.assertEqual(_astuple(small.compare(ref1, unrooted=True)),
#                      ("NA", "NA", 0.0, 1.0, 1.0, 3, 1, "NA"))

#     small = Tree("(A, B);")
#     # RF unrooted in too small trees
#     self.assertEqual(_astuple(small.compare(ref1, unrooted=True)),
#                      ("NA", "NA", 0.0, "NA", "NA", 2, 1, "NA"))

#     small = Tree("(A, B);")
#     # RF unrooted in too small trees
#     self.assertEqual(_astuple(small.compare(ref1, unrooted=False)),
#                      ("NA", "NA", 0.0, "NA", "NA", 2, 1, "NA"))

#     # identical trees, 8 rooted partitions in total (4 an 4), and 6 unrooted
#     self.assertEqual(_astuple(s1.compare(ref1)),
#                      (0.0, 0.0, 8, 1.0, 1.0, 6, 1, "NA"))

#     self.assertEqual(_astuple(s1.compare(ref1, unrooted=True)),
#                      (0.0, 0.0, 6, 1.0, 1.0, 6, 1, "NA"))

#     # The same stats should be return discarding branches, as the topology
#     # is still identical, but branches used should be different
#     self.assertEqual(_astuple(s1.compare(ref1, min_support_source=0.99, min_support_ref=.99)),
#                      (0.0, 0.0, 2, 1.0, 1.0, 6, 1, "NA"))

#     self.assertEqual(_astuple(s1.compare(ref1, min_support_source=0.99, min_support_ref=.99, unrooted=True)),
#                      (0.0, 0.0, 2, 1.0, 1.0, 6, 1, "NA"))


#     self.assertEqual(_astuple(s1.compare(ref1, min_support_source=0.99)),
#                      (0.0, 0.0, 5, 1/4., 1.0, 6, 1, "NA"))


#     self.assertEqual(_astuple(s1.compare(ref1, min_support_source=0.99, unrooted=True)),
#                      (0.0, 0.0, 4, 6/8., 1.0, 6, 1, "NA"))


#     self.assertEqual(_astuple(s1.compare(ref1, min_support_ref=0.99)),
#                      (0.0, 0.0, 5, 1.0, 1/4., 6, 1, "NA"))


#     self.assertEqual(_astuple(s1.compare(ref1, min_support_ref=0.99, unrooted=True)),
#                      (0.0, 0.0, 4, 1.0, 6/8., 6, 1, "NA"))


#     # Three partitions different
#     s2 = Tree('(((A, E)0.9, (C, D))0.98, (B,F)0.95);')
#     self.assertEqual(_astuple(s2.compare(ref1)),
#                      (6/8., 6, 8, 1/4., 1/4., 6, 1, "NA"))

#     self.assertEqual(_astuple(s2.compare(ref1, unrooted=True)),
#                      (4/6., 4, 6, 6/8., 6/8., 6, 1, "NA"))

#     # lets discard one branch from source tree.  there are 4 valid edges in
#     # ref, 3 in source there is only 2 edges in common, CD and root (which
#     # should be discounted for % of found branches)
#     self.assertEqual(_astuple(s2.compare(ref1, min_support_source=0.95)),
#                      (5/7., 5, 7, 1/4., 1/3., 6, 1, "NA"))

#     # similar in unrooted, but we don not need to discount root edges
#     self.assertEqual(_astuple(s2.compare(ref1, min_support_source=0.95, unrooted=True)),
#                      (3/5., 3, 5, 6/8., 6/7., 6, 1, "NA"))


#     # totally different trees
#     s3 = Tree('(((A, C)0.9, (E, D))0.98, (B,F)0.95);')
#     self.assertEqual(_astuple(s3.compare(ref1)),
#                      (1.0, 8, 8, 0.0, 0.0, 6, 1, "NA"))

# TODO: Merge this function with the previous one? (test_tree_compare())
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

    # usisng expand_polytomies and correct_by_size at the same time
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
    assert is_mono == False
    assert monotype == 'paraphyletic'

    is_mono, monotype, extra= t.check_monophyly(values=['a', 'e', 'i', 'o'])
    assert is_mono == True
    assert monotype == 'monophyletic'

    is_mono, monotype, extra =  t.check_monophyly(values=['i', 'o'])
    assert is_mono == False
    assert monotype == 'paraphyletic'

    # Now with unrooted trees, and using species instead of names.
    t = PhyloTree('(aaa1, (aaa3, (aaa4, (bbb1, bbb2))));',
                  sp_naming_function=lambda name: name[:3])
    # ─┬╴aaa1        # the species will be 'aaa' for this node, etc.
    #  ╰─┬╴aaa3
    #    ╰─┬╴aaa4
    #      ╰─┬╴bbb1
    #        ╰╴bbb2

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

    # # Check monophyly randomization test
    # t = PhyloTree(,
    # t.populate(100, dist_fn=rng.random, support_fn=rng.random)
    # ancestor = t.common_ancestor(['aaaaaaaaaa', 'aaaaaaaaab', 'aaaaaaaaac'])
    # all_nodes = list(t.descendants())
    # # I test every possible node as root for the tree. The content of ancestor
    # # should allways be detected as monophyletic
    # results = set()
    # for x in all_nodes:
    #     mono, part, extra = t.check_monophyly(values=set(ancestor.leaf_names()), prop='name', unrooted=True)
    #     results.add(mono)
    #     t.set_outgroup(x)
    # self.assertEqual(list(results), [True])

    # TODO: The previous test looks like it is wrong. Review with Jaime.

    # Testing get_monophyly
    t = Tree('((((((4, e), i)M1, o),h), u), ((3, 4), (i, june))M2);',
             parser=1)
    # we annotate the tree using external data
    colors = {'a': 'red', 'e': 'green', 'i': 'yellow',
              'o': 'black', 'u': 'purple', '4':'green',
              '3': 'yellow', '1': 'white', '5': 'red',
              'june': 'yellow'}
    for leaf in t:
        leaf.add_props(color=colors.get(leaf.name, 'none'))
    green_yellow_nodes = {t['M1'], t['M2']}
    mono_nodes = t.get_monophyletic(values=['green', 'yellow'],
                                    prop='color')
    assert set(mono_nodes) == green_yellow_nodes

def test_copy():
    t = Tree("((A, B)Internal_1:0.7, (C, D)Internal_2:0.5)root:1.3;",
             parser=1)
    # we add a custom annotation to the node named A
    t['A'].add_props(label="custom Value")
    # we add a complex feature to the A node, consisting of a list of lists
    t['A'].add_props(complex=[[0,1], [2,3], [1,11], [1,0]])

    nw2 = t.write(props=None, format_root_node=True, parser=1)

    t_nw  = t.copy("newick")
    t_nwx = t.copy("newick-extended")
    t_pkl = t.copy("cpickle")
    t['A'].props['testfn'] = lambda: "YES"
    t_deep = t.copy("deepcopy")

    assert (t_nw["root"]).name == "root"
    assert (t_nwx["A"]).props['label'] == "custom Value"
    assert (t_pkl["A"]).props['complex'][0] == [0,1]
    assert (t_deep["A"]).props['testfn']() == "YES"

@pytest.mark.parametrize("newick, expected_dists, expected_leaves", COPHENETIC_MATRIX_CASES)
def test_cophenetic_matrix(newick, expected_dists, expected_leaves):
    t = Tree(newick)
    dists, leaves = t.cophenetic_matrix()
    for i in range(len(expected_dists)):
        for j in range(len(expected_dists[i])):
            assert round(abs(expected_dists[i][j] - dists[i][j]), 4) == 0
    assert expected_leaves == leaves
