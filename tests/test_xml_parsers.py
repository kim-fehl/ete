from pathlib import Path
import time
import pytest
from ete4 import phyloxml

EXAMPLE_PATH = Path(__file__).resolve().parent.parent / "examples" / "phyloxml"
XML_FILES = sorted(EXAMPLE_PATH.glob("*.xml"))


class _Sink:
    def write(self, _):
        pass


@pytest.mark.parametrize("xml_file", XML_FILES)
def test_phyloxml_parser(xml_file, tmp_path):
    print(xml_file.name, "...", end=" ")
    p = phyloxml.Phyloxml()
    t1 = time.time()
    p.build_from_file(str(xml_file))
    etime = time.time() - t1
    print(f"{etime:0.1f} secs")
    p.export(outfile=_Sink())



