"""Minimal script to reproduce potential Qt segmentation faults.

This script mirrors the tree rendering steps from the test suite. It
runs the rendering in an offscreen QApplication and exits immediately.
"""
import os

# Force headless Qt backend
os.environ.setdefault("QT_QPA_PLATFORM", "offscreen")
os.environ.setdefault("PYTHONFAULTHANDLER", "1")

from PyQt6.QtWidgets import QApplication
from PyQt6.QtCore import QObject

from ete4 import Tree
from ete4.treeview import TreeStyle


def main() -> None:
    app = QApplication([])
    t = Tree("(A,B);")
    ts = TreeStyle()
    # Render to a temporary file to exercise the Qt drawing pipeline
    t.render("/tmp/minimal.svg", tree_style=ts)
    # Print current object tree to spot lingering widgets
    QObject.dumpObjectTree(app)
    app.quit()


if __name__ == "__main__":
    main()
