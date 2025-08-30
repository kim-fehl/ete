import pytest
from ete4 import PhyloTree, SeqGroup


def test_link_alignments(sp_name_fn):
    """PhyloTree can be linked to SeqGroup objects"""
    fasta = """
         >seqA
         MAEIPDETIQQFMALT---HNIAVQYLSEFGDLNEALNSYYASQTDDIKDRREEAH
         >seqB
         MAEIPDATIQQFMALTNVSHNIAVQY--EFGDLNEALNSYYAYQTDDQKDRREEAH
         >seqC
         MAEIPDATIQ---ALTNVSHNIAVQYLSEFGDLNEALNSYYASQTDDQPDRREEAH
         >seqD
         MAEAPDETIQQFMALTNVSHNIAVQYLSEFGDLNEAL--------------REEAH
    """
    # Caution with iphylip string. blank spaces in the beginning are important
    iphylip = """
         4 76
      seqA   MAEIPDETIQ QFMALT---H NIAVQYLSEF GDLNEALNSY YASQTDDIKD RREEAHQFMA
      seqB   MAEIPDATIQ QFMALTNVSH NIAVQY--EF GDLNEALNSY YAYQTDDQKD RREEAHQFMA
      seqC   MAEIPDATIQ ---ALTNVSH NIAVQYLSEF GDLNEALNSY YASQTDDQPD RREEAHQFMA
      seqD   MAEAPDETIQ QFMALTNVSH NIAVQYLSEF GDLNEAL--- ---------- -REEAHQ---

             LTNVSHQFMA LTNVSH
             LTNVSH---- ------
             LTNVSH---- ------
             -------FMA LTNVSH
    """
    # Loads a tree and link it to an alignment. As usual, 'alignment' can be
    # the path to a file or the data themselves in text string format
    alg1 = SeqGroup(fasta)
    alg2 = SeqGroup(iphylip, format="iphylip")
    t = PhyloTree("(((seqA,seqB),seqC),seqD);", alignment=fasta, alg_format="fasta",
                  sp_naming_function=sp_name_fn)
    for leaf in t.leaves():
        assert leaf.props.get('sequence') == alg1.get_seq(leaf.name)
    # The associated alignment can be changed at any time
    t.link_to_alignment(alignment=alg2, alg_format="iphylip")
    for leaf in t.leaves():
        assert leaf.props.get('sequence') == alg2.get_seq(leaf.name)
