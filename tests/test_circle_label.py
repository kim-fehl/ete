import math

import pytest

from pathlib import Path

pytest.importorskip("PyQt6.QtGui")
from PyQt6.QtCore import QIODevice

# Qt6 removed WriteOnly/ReadOnly attributes; alias them for compatibility
if not hasattr(QIODevice, "WriteOnly"):
    QIODevice.WriteOnly = QIODevice.OpenModeFlag.WriteOnly
    QIODevice.ReadOnly = QIODevice.OpenModeFlag.ReadOnly

pytestmark = pytest.mark.skip("Circle label rendering segfaults headlessly")

from ete4 import Tree, CircleFace, RectFace, TextFace, TreeStyle, add_face_to_node
from .conftest import normalize_svg

NEWLINE = "\n"
TAB = "\t"
UNDERBAR = "_"
SPACE = " "
COMMA = ","

_circle_tested = False


def node_layout(annotation):
    def _node_layout(node):
        global _circle_tested
        node.img_style["size"] = 0

        if node.is_leaf:
            if node.name not in annotation:
                print(f'Got unknown leaf "{node.name}"')
                return
            size = annotation[node.name]
        else:
            children = COMMA.join(sorted(leaf.name for leaf in node))
            if children not in annotation:
                print("Got unknown node " + children)
                return
            size = annotation[children]
        dimension = 4 * math.sqrt(size)
        labeld = {
            "text": f"{size:d}",
            "color": "white",
            "font": "Helvetica",
            "size": 8,
        }
        if size % 2:
            clabel = dict(labeld)
            if size == 35:
                if not _circle_tested:
                    _circle_tested = True
                    clabel["text"] = "Ly"
                clabel["size"] = 12
                thisFace = CircleFace(
                    dimension / 2, "steelblue", "circle", label=clabel
                )
            elif size == 43:
                clabel["size"] = 6
                del clabel["color"]
                thisFace = CircleFace(
                    dimension / 2, "steelblue", "sphere", label=clabel
                )
            else:
                thisFace = CircleFace(
                    dimension / 2, "steelblue", "sphere", label=f"{size:d}"
                )
        else:
            thisFace = RectFace(dimension, dimension, "green", "blue", label=labeld)
        thisFace.opacity = 0.7
        add_face_to_node(thisFace, node, column=0, position="float")
        textF = TextFace(str(size), fsize=12, fgcolor="steelblue")
        add_face_to_node(textF, node, column=0, position="aligned")

    return _node_layout


TREE = '''
(S032           :0.68575,(S032 473      :0.29028,(S032 476      :0.19027,
(S032 477       :0.21930,(S032 478      :0.15731,(S032 475      :0.08076,
S032 474        :0.08914):0.04624):0.04277):0.01821):0.07911):0.02825,S032 202  :0.31425);
'''

SETS = '''
  43 S032_477
  42 S032_202
  35 S032_476
  35 S032_474
  33 S032_475
  32 S032_478
  29 S032_473
  25 S032_474,S032_475
  20 S032_476,S032_477
  20 S032_475,S032_478
  20 S032_475,S032_477
  20 S032_474,S032_478
  20 S032_474,S032_476
  19 S032_474,S032_475,S032_478
  19 S032_474,S032_475,S032_477
  19 S032_475,S032_476
  18 S032_474,S032_475,S032_476
  18 S032_477,S032_478
  17 S032_475,S032_477,S032_478
  17 S032_476,S032_478
  16 S032_474,S032_475,S032_477,S032_478
  16 S032_474,S032_476,S032_478
  15 S032_474,S032_475,S032_476,S032_478
  15 S032_474,S032_475,S032_476,S032_477
  15 S032_473,S032_477
  14 S032_473,S032_476
  13 S032_474,S032_475,S032_476,S032_477,S032_478
  12 S032_202,S032_475
  11 S032_473,S032_475,S032_476
  11 S032_202,S032_474,S032_475
  10 S032_473,S032_474,S032_475,S032_476,S032_478
  10 S032_202,S032_474,S032_475,S032_477,S032_478
   9 S032_473,S032_476,S032_477
   8 S032_473,S032_474,S032_475,S032_476,S032_477,S032_478
   8 S032_202,S032_473,S032_475,S032_476
   7 S032_202,S032_473,S032_474,S032_475,S032_476,S032_477,S032_478
'''

treesets = [s for s in SETS.strip().split(NEWLINE) if s.strip()]
thesets = {}
for s in treesets:
    if TAB in s:
        num, namelist = s.strip().split(TAB)
    else:
        num, namelist = s.strip().split(None, 1)
    thesets[namelist.replace(UNDERBAR, SPACE)] = int(num)


def test_circle_label(qapp, file_regression):
    global _circle_tested
    _circle_tested = False
    t = Tree(TREE)
    ts = TreeStyle()
    ts.layout_fn = node_layout(thesets)
    ts.show_leaf_name = True
    ts.show_branch_length = True

    data = t.render("%%return.svg", tree_style=ts)[0]
    svg = eval(data).decode().lower()
    normalized = normalize_svg(svg)

    baseline = (
        Path(__file__).parent / "data" / "expected" / "test_circle_label.svg"
    )
    file_regression.check(normalized, extension=".svg", fullpath=baseline)

