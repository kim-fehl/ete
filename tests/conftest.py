import os
import pytest
from ete4 import SeqGroup


@pytest.fixture
def expected_seqgroup():
    path = os.path.join(os.path.dirname(__file__), "data", "expected_sequences.fa")
    return SeqGroup(path)
