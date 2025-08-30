import pytest

pytest.importorskip("PyQt6.QtGui")
from ete4 import Tree
from ete4.treeview import TreeStyle
from PyQt6.QtCore import QIODevice
from pathlib import Path

from .conftest import normalize_svg

# PyQt6 changed enum locations; provide backwards-compatible aliases
if not hasattr(QIODevice, "WriteOnly"):
    QIODevice.WriteOnly = QIODevice.OpenModeFlag.WriteOnly
    QIODevice.ReadOnly = QIODevice.OpenModeFlag.ReadOnly


def _render_svg(tree, ts):
    """Render *tree* using *ts* and return the raw SVG string."""
    try:
        # Newer API
        svg = tree.render(tree_style=ts, format="svg", render_string=True)
    except TypeError:
        # Older API uses special filename token to return a string
        svg, _ = tree.render("%%return.svg", tree_style=ts)
    if not isinstance(svg, str):
        svg = svg[0]
    if svg.startswith("b'"):
        import ast
        svg = ast.literal_eval(svg).decode("utf-8")
    return svg

def test_svg_regressions(qapp, file_regression):
    tree = Tree("((A,B),C);")
    ts = TreeStyle()
    ts.show_leaf_name = True

    svg = _render_svg(tree, ts)
    normalized = normalize_svg(svg)
    baseline = (
        Path(__file__).parent
        / "data"
        / "expected"
        / "test_treeview_svg_regressions.svg"
    )
    file_regression.check(normalized, extension=".svg", fullpath=baseline)
