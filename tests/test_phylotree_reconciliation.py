from ete4 import PhyloTree


def test_reconciliation(sp_name_fn):
    """Tests orthology prediction based on species reconciliation"""
    gene_tree_nw = '((Dme_001,Dme_002),(((Cfa_001,Mms_001),((Hsa_001,Ptr_001),Mmu_001)),(Ptr_002,(Hsa_002,Mmu_002))));'
    species_tree_nw = "((((Hsa, Ptr), Mmu), (Mms, Cfa)), Dme);"

    genetree = PhyloTree(gene_tree_nw, sp_naming_function=sp_name_fn)
    sptree = PhyloTree(species_tree_nw, sp_naming_function=sp_name_fn)

    recon_tree, events = genetree.reconcile(sptree)
    # Check that reconcilied tree nodes have the correct lables:
    # gene loss, duplication, etc.
    expected_recon = "((Dme_001:1,Dme_002:1)1:1[&&NHX:evoltype=D],(((Cfa_001:1,Mms_001:1)1:1[&&NHX:evoltype=S],((Hsa_001:1,Ptr_001:1)1:1[&&NHX:evoltype=S],Mmu_001:1)1:1[&&NHX:evoltype=S])1:1[&&NHX:evoltype=S],((Mms:1[&&NHX:evoltype=L],Cfa:1[&&NHX:evoltype=L])1:1[&&NHX:evoltype=L],(((Hsa:1[&&NHX:evoltype=L],Ptr_002:1)1:1[&&NHX:evoltype=L],Mmu:1[&&NHX:evoltype=L])1:1[&&NHX:evoltype=L],((Ptr:1[&&NHX:evoltype=L],Hsa_002:1)1:1[&&NHX:evoltype=L],Mmu_002:1)1:1[&&NHX:evoltype=S])1:1[&&NHX:evoltype=D])1:1[&&NHX:evoltype=L])1:1[&&NHX:evoltype=D])[&&NHX:evoltype=S];"
    assert recon_tree.write(props=["evoltype"], parser=9) == PhyloTree(expected_recon).write(props=["evoltype"], parser=9)
