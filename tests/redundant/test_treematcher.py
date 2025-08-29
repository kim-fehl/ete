"""
Tests related to the treematcher module.
"""

import unittest

from ete4 import Tree
import ete4.treematcher as tm
from tests.conftest import TREEMATCHER_CASES



def strip(text):
    """Return the given text stripping the empty lines and indentation."""
    # Helps compare tree visualizations.
    indent = min(len(line) - len(line.lstrip())
                 for line in text.splitlines() if line.strip())
    return '\n'.join(line[indent:].rstrip()
        for line in text.splitlines() if line.strip())


class TestTreematcher(unittest.TestCase):

    def test_str(self):

        # See if the representation of a pattern as a text is the expected one.
        pattern = tm.TreePattern("""
        (
          "len(ch) > 2",
          "name in ['hello', 'bye']"
        )
        "(len(name) < 3 or name == 'accept') and d >= 0.5"
        """)
        self.assertEqual(str(pattern), strip("""
                                                  ╭╴len(ch) > 2
╴(len(name) < 3 or name == 'accept') and d >= 0.5╶┤
                                                  ╰╴name in ['hello', 'bye']
        """))

        # See if we can use quotes in a different way, and format more widely.
        pattern2 = tm.TreePattern("""
        (
        '  len(ch) > 2  '  ,
        '  name in ["hello", "bye"]'
        )
        '(len(name) < 3 or name == "accept") and d >= 0.5  '
        """)
        self.assertEqual(str(pattern), str(pattern2).replace('"', "'"))

    def test_search(self):
        pattern = tm.TreePattern("""
        (
          "len(ch) > 2",
          "name in ['hello', 'bye']"
        )
        "(len(name) < 3 or name == 'accept') and d >= 0.5"
        """)

        for newick, expected_result in TREEMATCHER_CASES:
            tree = Tree(newick, parser=1)

            self.assertEqual([n.name for n in tm.search(pattern, tree)],
                             expected_result)
            self.assertEqual([n.name for n in pattern.search(tree)],
                             expected_result)
