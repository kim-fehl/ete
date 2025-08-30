import pytest

pytest.importorskip("PyQt6.QtGui")
from PyQt6.QtCore import QIODevice

# Qt6 removed WriteOnly/ReadOnly attributes; alias them for compatibility
if not hasattr(QIODevice, "WriteOnly"):
    QIODevice.WriteOnly = QIODevice.OpenModeFlag.WriteOnly
    QIODevice.ReadOnly = QIODevice.OpenModeFlag.ReadOnly

from pathlib import Path
from ete4 import PhyloTree
from ete4.treeview import TreeStyle, faces
from ete4.treeview.faces import SequenceFace

from .conftest import normalize_svg

def test_alignment_face_colors_and_order(qapp, file_regression):
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
    normalized = normalize_svg(svg)

    baseline = (
        Path(__file__).parent
        / "data"
        / "expected"
        / "test_treeview_alignment_face_colors_and_order.svg"
    )
    file_regression.check(normalized, extension=".svg", fullpath=baseline)
