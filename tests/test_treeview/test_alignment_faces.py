import os
import re

import pytest

# Ensure Qt can operate in headless environments
os.environ.setdefault("QT_QPA_PLATFORM", "offscreen")

from PyQt6.QtCore import QIODevice
from PyQt6.QtWidgets import QApplication

# Qt6 removed WriteOnly/ReadOnly attributes; alias them for compatibility
if not hasattr(QIODevice, "WriteOnly"):
    QIODevice.WriteOnly = QIODevice.OpenModeFlag.WriteOnly
    QIODevice.ReadOnly = QIODevice.OpenModeFlag.ReadOnly

from ete4 import PhyloTree
from ete4.treeview import TreeStyle, faces
from ete4.treeview.faces import SequenceFace


@pytest.fixture(scope="module")
def qapp():
    app = QApplication([])
    yield app


def test_alignment_face_colors_and_order(qapp):
    alignment = ">A\nACGT\n>B\nTGCA\n"
    tree = PhyloTree("(A,B);")
    tree.link_to_alignment(alignment, alg_format="fasta")

    ts = TreeStyle()

    def layout(node):
        seq = node.props.get("sequence")
        if node.is_leaf and seq:
            seqface = SequenceFace(seq, seqtype="nt", fsize=10)
            faces.add_face_to_node(seqface, node, column=1, aligned=True)

    ts.layout_fn = layout

    data = tree.render("%%return.svg", tree_style=ts)[0]
    svg = eval(data).decode().lower()

    # Colors for first sequence ACGT
    assert re.search(
        r"#a0a0ff.*?#ff8c4b.*?#ff7070.*?#a0ffa0", svg, re.DOTALL
    )
    # Colors for second sequence TGCA
    assert re.search(
        r"#a0ffa0.*?#ff7070.*?#ff8c4b.*?#a0a0ff", svg, re.DOTALL
    )
