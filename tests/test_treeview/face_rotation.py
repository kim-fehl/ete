import random

from ete4 import Tree
from ete4.treeview import TreeStyle, add_face_to_node, TextFace


def rotation_layout(node, rng):
    if node.is_leaf:
        F = TextFace(node.name, tight_text=True)
        F.rotation = rng.randint(0, 360)
        add_face_to_node(TextFace("third"), node, column=8, position="branch-right")
        add_face_to_node(TextFace("second"), node, column=2, position="branch-right")
        add_face_to_node(F, node, column=0, position="branch-right")

        F.border.width = 1
        F.inner_border.width = 1

def get_example_tree(rng=None):
    if rng is None:
        rng = random
    t = Tree()
    t.populate(10, dist_fn=rng.random, support_fn=rng.random)
    ts = TreeStyle()
    ts.rotation = 45
    ts.show_leaf_name = False
    ts.layout_fn = lambda node: rotation_layout(node, rng)

    return t, ts

if __name__ == "__main__":
    t, ts = get_example_tree()
    t.show(tree_style=ts)
