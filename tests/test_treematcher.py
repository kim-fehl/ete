"""
Tests related to the treematcher module.
"""

from ete4 import Tree, PhyloTree
import ete4.treematcher as tm

import pytest
from .conftest import TREEMATCHER_CASES


def strip(text):
    """Return the given text stripping the empty lines and indentation."""
    # Helps compare tree visualizations.
    indent = min(len(line) - len(line.lstrip())
                 for line in text.splitlines() if line.strip())
    return '\n'.join(line[indent:].rstrip()
        for line in text.splitlines() if line.strip())


def test_str():
    # See if the representation of a pattern as a text is the expected one.
    pattern = tm.TreePattern("""
    (
      "len(ch) > 2",
      "name in ['hello', 'bye']"
    )
    "(len(name) < 3 or name == 'accept') and d >= 0.5"
    """)
    assert str(pattern) == strip("""
                                                  ╭╴len(ch) > 2
╴(len(name) < 3 or name == 'accept') and d >= 0.5╶┤
                                                  ╰╴name in ['hello', 'bye']
    """)

    # See if we can use quotes in a different way, and format more widely.
    pattern2 = tm.TreePattern("""
    (
    '  len(ch) > 2  '  ,
    '  name in ["hello", "bye"]'
    )
    '(len(name) < 3 or name == "accept") and d >= 0.5  '
    """)
    assert str(pattern) == str(pattern2).replace('"', "'")


@pytest.mark.parametrize("newick, expected", TREEMATCHER_CASES)
def test_search(newick, expected):
    pattern = tm.TreePattern("""
    (
      "len(ch) > 2",
      "name in ['hello', 'bye']"
    )
    "(len(name) < 3 or name == 'accept') and d >= 0.5"
    """)

    tree = Tree(newick, parser=1)

    assert ([n.name for n in tm.search(pattern, tree)] ==
            [n.name for n in pattern.search(tree)] ==
            expected)


def test_safer():
    t = PhyloTree('(a,(b,c));', sp_naming_function=lambda name: name)

    tp_unsafe = tm.TreePattern('("node.get_species()=={\'c\'}",'
                               '  node.species=="b")')
    assert list(tp_unsafe.search(t)) == [t.common_ancestor(['b', 'c'])]

    tp_safer = tm.TreePattern('("node.get_species()=={\'c\'}",'
                              '  node.species=="b")', safer=True)
    with pytest.raises(SyntaxError):
        list(tp_safer.search(t))  # asked for unknown function get_species()
