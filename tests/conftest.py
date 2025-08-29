import os
import pytest
from ete4 import SeqGroup
from ete4.core.tree import Tree
from ete4 import PhyloTree


TREEMATCHER_CASES = [
    ('((hello:1,(1:1,2:1,3:1)xx:1)accept:1, NODE):0;', ['accept']),
    ('((hello:1,(1:1,2:1,3:1)xx:1)accept:0.4, NODE):0;', []),
    ('(hello:1,(1:1,2:1,3:1)xx:1)accept:1;', ['accept']),
    ('((bye:1,(1:1,2:1,3:1)xx:1)none:1, NODE):0;', []),
    ('((bye:1,(1:1,2:1,3:1)xx:1)y:1, NODE):0;', ['y']),
    ('((bye,(,,))x:1,((,,),bye)y:1):0;', ['x', 'y'])
]

COMPLEX_NAME = "((A:0.0001[&&NHX:hello=true],B:0.011)90:0.01[&&NHX:hello=true],(C:0.01, D:0.001)hello:0.01);"

QUOTED_NAME_CASES = [
    (
        '(("A:0.1":1,"%s":2)"C:0.00":3,"D":4);' % COMPLEX_NAME,
        "(('A:0.1':1,'%s':2)'C:0.00':3,D:4);" % COMPLEX_NAME,
    ),
    (
        '''(("A:\\"0.1\\"":1,"%s":2)"C:'0.00'":3,"D'sd'x":4);'''
        % COMPLEX_NAME,
        '''(('A:\\"0.1\\"':1,'%s':2)'C:''0.00''\':3,'D''sd''x':4);'''
        % COMPLEX_NAME,
    ),
]

CUSTOM_FORMAT_CASES = [
    (0, '((TEST-A:1.1,TEST-B:2.2)SUP-1.0:3.3,TEST-D:4.4);'),
    (1, '((TEST-A:1.1,TEST-B:2.2)TEST-C:3.3,TEST-D:4.4);'),
    (2, '((TEST-A:1.1,TEST-B:2.2)SUP-1.0:3.3,TEST-D:4.4);'),
    (3, '((TEST-A:1.1,TEST-B:2.2)TEST-C:3.3,TEST-D:4.4);'),
    (4, '((TEST-A:1.1,TEST-B:2.2),TEST-D:4.4);'),
    (5, '((TEST-A:1.1,TEST-B:2.2):3.3,TEST-D:4.4);'),
    (6, '((TEST-A,TEST-B):3.3,TEST-D);'),
    (7, '((TEST-A:1.1,TEST-B:2.2)TEST-C,TEST-D:4.4);'),
    (8, '((TEST-A,TEST-B)TEST-C,TEST-D);'),
    (9, '((TEST-A,TEST-B),TEST-D);'),
]


@pytest.fixture
def example1_newick():
    return '(((ao:1,(ap:1,aq:1)1:1)1:1,(ar:1,(as:1,at:1)1:1)1:1)1:1,((aa:1,ab:1)1:1,((ac:1,(ad:1,(ae:1,(af:1,(ag:1,ah:1)1:1)1:1)1:1)1:1)1:1,((ai:1,(aj:1,(ak:1,al:1)1:1)1:1)1:1,(am:1,an:1)1:1)1:1)1:1)1:1);'


@pytest.fixture
def example2_newick():
    return '(((2ao:1,(2ap:1,2aq:1)1:1)1:1,(2ar:1,(2as:1,2at:1)1:1)1:1)1:1,((2aa:1,2ab:1)1:1,((2ac:1,(2ad:1,(2ae:1,(2af:1,(2ag:1,2ah:1)1:1)1:1)1:1)1:1)1:1,((2ai:1,(2aj:1,(2ak:1,2al:1)1:1)1:1)1:1,(2am:1,2an:1)1:1)1:1)1:1)1:1);'


@pytest.fixture
def example3_newick():
    return '(((ao:1,(ap:1,2aq:1)1:1)1:1,(ar:1,(as:1,2at:1)1:1)1:1)1:1,((aa:1,ab:1)1:1,((2ac:1,(2ad:1,(2ae:1,(2af:1,(ag:1,2ah:1)1:1)1:1)1:1)1:1)1:1,((2ai:1,(aj:1,(2ak:1,al:1)1:1)1:1)1:1,(2am:1,an:1)1:1)1:1)1:1)1:1);'


@pytest.fixture
def example4_newick():
    return '(((2ao:1,(ap:1,aq:1)1:1)1:1,(2ar:1,(2as:1,at:1)1:1)1:1)1:1,((aa:1,2ab:1)1:1,((ac:1,(ad:1,(ae:1,(2af:1,(2ag:1,2ah:1)1:1)1:1)1:1)1:1)1:1,((ai:1,(aj:1,(2ak:1,al:1)1:1)1:1)1:1,(2am:1,an:1)1:1)1:1)1:1)1:1);'


@pytest.fixture
def example1_tree(example1_newick):
    return Tree(example1_newick)


@pytest.fixture
def example2_tree(example2_newick):
    return Tree(example2_newick)


@pytest.fixture
def example3_tree(example3_newick):
    return Tree(example3_newick)


@pytest.fixture
def example4_tree(example4_newick):
    return Tree(example4_newick)


@pytest.fixture
def phylotree_example_newick():
    return '((Dme_001:1,Dme_002:1):1,(((Cfa_001:1,Mms_001:1):1,((((Hsa_001:1,Hsa_003:1):1,Ptr_001:1):1,Mmu_001:1):1,((Hsa_004:1,Ptr_004:1):1,Mmu_004:1):1):1):1,(Ptr_002:1,(Hsa_002:1,Mmu_002:1):1):1):1):0;'


@pytest.fixture
def phylotree_example(phylotree_example_newick):
    return PhyloTree(phylotree_example_newick)


@pytest.fixture
def expected_seqgroup():
    path = os.path.join(os.path.dirname(__file__), "data", "expected_sequences.fa")
    return SeqGroup(path)
