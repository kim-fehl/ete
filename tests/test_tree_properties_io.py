from tempfile import NamedTemporaryFile

import pytest

from ete4 import Tree
from ete4.core.tree import TreeError
from ete4.parser.newick import NewickError
from ete4.parser import newick

from . import conftest as ds
from .conftest import CUSTOM_FORMAT_CASES, QUOTED_NAME_CASES, COMPLEX_NAME


def strip(text):
    """Return the given text stripping the empty lines and indentation."""
    indent = min(len(line) - len(line.lstrip())
                 for line in text.splitlines() if line.strip())
    return '\n'.join(line[indent:].rstrip()
                     for line in text.splitlines() if line.strip())


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
            assert t.write(parser=f) == Tree(t.write(parser=f), parser=f).write(parser=f)

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
    for f in [0, 1, 2, 3, 5, 6, 7, 8, 9]:
        t2 = Tree(t.write(parser=f), parser=f)
        t2.sort_descendants()
        observed_names = [n.name for n in t]
        assert observed_names == expected_leaf_names

    # Check that all formats reading distances, recover original distances
    for f in [0, 1, 2, 3, 5]:
        t2 = Tree(t.write(parser=f), parser=f)
        t2.dist, t2.support = 0, 1
        t2.sort_descendants()
        observed_distances = [n.dist for n in t2.traverse('postorder')]
        assert all(abs(x - y) < 1e-6 for x, y in zip(observed_distances, expected_distances))

    # formats reading only leaf distances
    for f in [4, 7]:
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
    for f in [0, 2]:
        t2 = Tree(t.write(parser=f), parser=f)
        t2.dist, t2.support = 0, 1
        t2.sort_descendants()
        observed_supports = [n.support for n in t2.traverse('postorder') if not n.is_leaf]
        assert all(abs(x - y) < 1e-6 for x, y in zip(observed_supports, expected_supports))

    # Check that formats reading supports, do not accept node names
    for f in [0, 2]:
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
        assert set([n.name for n in Tree(base_nw.replace('NAME2', vname), parser=1)]) == expected_names

    # invalid NHX format
    with pytest.raises(NewickError):
        Tree("(((A, B), C)[&&NHX:nameI]);")
    # unsupported newick stream
    with pytest.raises(Exception):
        Tree([1, 2, 3])


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
    t = Tree('((A:1.1111,B:2.2222)C:3.3333[&&NHX:support=1],D:4.4444);', parser=1)
    t.sort_descendants()

    parser = newick.make_parser(flag, dist='%0.1f', name='TEST-%s',
                                support='SUP-%0.1f')
    nw = t.write(parser=parser)
    assert nw == result


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
    assert t['a'].id == (0, 0)
    assert t['d'].id == (1, 1)


def test_copy():
    t = Tree("((A, B)Internal_1:0.7, (C, D)Internal_2:0.5)root:1.3;", parser=1)
    # we add a custom annotation to the node named A
    t['A'].add_props(label="custom Value")
    # we add a complex feature to the A node, consisting of a list of lists
    t['A'].add_props(complex=[[0, 1], [2, 3], [1, 11], [1, 0]])

    t_nw = t.copy("newick")
    t_nwx = t.copy("newick-extended")
    t_pkl = t.copy("cpickle")
    t['A'].props['testfn'] = lambda: "YES"
    t_deep = t.copy("deepcopy")

    assert (t_nw["root"]).name == "root"
    assert (t_nwx["A"]).props['label'] == "custom Value"
    assert (t_pkl["A"]).props['complex'][0] == [0, 1]
    assert (t_deep["A"]).props['testfn']() == "YES"

