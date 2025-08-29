import os
import sys
import pytest

APPSPATH = os.path.join(sys.prefix, "bin", "ete4_apps", "bin")
if not os.path.isfile(os.path.join(APPSPATH, "statal")):
    pytest.skip("external ete4 apps not installed", allow_module_level=True)

from .__main__ import *
