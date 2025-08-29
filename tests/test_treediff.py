
from ete4.tools import ete_diff as ediff


def almost_equal(x, y, precision=1e-6):
    return abs(x - y) / max(abs(x), abs(y)) < precision


class Test_Treediff:
    """Test specific methods for trees linked to treediff."""

    def test_treediff_basic(self, example1_tree):
        """Test tree-diff basic functionality."""
        t1 = example1_tree
        t2 = example1_tree

        difftable = ediff.treediff(t1, t2, prop1='name', prop2='name',
                                   dist_fn=ediff.EUCL_DIST, support=False,
                                   reduce_matrix=False, extended=None,
                                   jobs=1, parallel=None)

        assert isinstance(difftable, list)
        assert isinstance(difftable[0], list)
        assert len(difftable[0]) == 7
        assert len(difftable) == 39

    def test_treediff_EUCL_DIST_1(self, example1_tree, example2_tree):
        """Test tree-diff EUCL_DIST distance."""
        t1 = example1_tree
        t2 = example2_tree

        difftable = ediff.treediff(t1, t2, prop1='name', prop2='name',
                                   dist_fn=ediff.EUCL_DIST, support=False,
                                   reduce_matrix=False, extended=None,
                                   jobs=1, parallel=None)

        assert sum(i[0] for i in difftable) == 39

    def test_treediff_EUCL_DIST_2(self, example1_tree, example3_tree):
        """ Tests tree-diff EUCL_DIST distance."""
        t1 = example1_tree
        t2 = example3_tree

        difftable = ediff.treediff(t1, t2, prop1='name', prop2='name',
                                   dist_fn=ediff.EUCL_DIST, support=False,
                                   reduce_matrix=False, extended=None,
                                   jobs=1, parallel=None)

        assert almost_equal(sum(i[0] for i in difftable), 19.621428)

    def test_treediff_EUCL_DIST_3(self, example1_tree, example3_tree):
        """ Tests tree-diff EUCL_DIST diffs"""
        t1 = example1_tree
        t2 = example3_tree

        difftable = ediff.treediff(t1, t2, prop1='name', prop2='name',
                                   dist_fn=ediff.EUCL_DIST, support=False,
                                   reduce_matrix=False, extended=None,
                                   jobs=1, parallel=None)

        assert sorted(i[4] for i in difftable) == DIFFS

    def test_treediff_RF_DIST_1(self, example1_tree, example2_tree):
        """ Tests tree-diff RF_DIST distance"""
        t1 = example1_tree
        t2 = example2_tree

        difftable = ediff.treediff(t1, t2, prop1='name', prop2='name',
                                   dist_fn=ediff.RF_DIST, support=False,
                                   reduce_matrix=False, extended=None,
                                   jobs=1, parallel=None)

        assert sum(i[0] for i in difftable) == 39.0

    def test_treediff_RF_DIST_2(self, example1_tree, example3_tree):
        """ Tests tree-diff RF_DIST distance"""
        t1 = example1_tree
        t2 = example3_tree

        difftable = ediff.treediff(t1, t2, prop1='name', prop2='name',
                                   dist_fn=ediff.RF_DIST, support=False,
                                   reduce_matrix=False, extended=None,
                                   jobs=1, parallel=None)

        assert sum(i[0] for i in difftable) == 10.0

    def test_treediff_extendend_cc(self, example1_tree, example3_tree):
        """ Tests tree-diff Extended distance. Cophenetic Compared"""
        t1 = example1_tree
        t2 = example3_tree

        difftable = ediff.treediff(t1, t2, prop1='name', prop2='name',
                                   dist_fn=ediff.RF_DIST, support=False,
                                   reduce_matrix=False,
                                   extended=ediff.cc_distance,
                                   jobs=1, parallel=None)

        assert sum(i[1] for i in difftable) == 863.9737020175473

    def test_treediff_extendend_be(self, example1_tree, example3_tree):
        """ Tests tree-diff  Extended distance. Branch Extended"""
        t1 = example1_tree
        t2 = example3_tree

        difftable = ediff.treediff(t1, t2, prop1='name', prop2='name',
                                   dist_fn=ediff.RF_DIST, support=False,
                                   reduce_matrix=False,
                                   extended=ediff.be_distance,
                                   jobs=1, parallel=None)

        assert sum(i[1] for i in difftable) == 616.0

    def test_treediff_reports(self, example1_tree, example4_tree):
        """ Tests tree-diff Reports"""
        t1 = example1_tree
        t2 = example4_tree

        difftable = ediff.treediff(t1, t2, prop1='name', prop2='name',
                                   dist_fn=ediff.EUCL_DIST, support=False,
                                   reduce_matrix=False,
                                   extended=ediff.cc_distance,
                                   jobs=1, parallel=None)

        rf , rf_max = t1.robinson_foulds(t2)[:2]

        # NOTE: This is not a test. Calling the functions below just print
        #       stuff on the screen, but there are no checks. Commenting out.
        # ediff.show_difftable_summary(difftable, rf=rf, rf_max=rf_max, extended=ediff.cc_distance))
        # ediff.show_difftable(difftable, extended=ediff.cc_distance)
        # ediff.show_difftable_tab(difftable, extended=ediff.cc_distance)
        # ediff.show_difftable_topo(difftable, 'name', 'name', usecolor=False, extended=ediff.cc_distance)



DIFFS = [
    set(),
    set(),
    set(),
    set(),
    set(),
    set(),
    set(),
    set(),
    set(),
    set(),
    set(),
    {'2aq', 'aq'},
    {'2ai', 'at'},
    {'2ah', 'ac'},
    {'2af', 'ad'},
    {'2ae', 'ae'},
    {'2ad', 'af'},
    {'2ac', 'ah'},
    {'2at', 'ai'},
    {'2am', 'ak'},
    {'2ak', 'am'},
    {'2aq', 'aq'},
    {'2at', 'at'},
    {'2ah', 'ah'},
    {'2ak', 'ak'},
    {'2am', 'am'},
    {'2aq', 'aq'},
    {'2at', 'at'},
    {'2af', '2ah', 'af', 'ah'},
    {'2ak', 'ak'},
    {'2ae', '2af', '2ah', 'ae', 'af', 'ah'},
    {'2ai', '2ak', 'ai', 'ak'},
    {'2ad', '2ae', '2af', '2ah', 'ad', 'ae', 'af', 'ah'},
    {'2aq', '2at', 'aq', 'at'},
    {'2ac', '2ad', '2ae', '2af', '2ah', 'ac', 'ad', 'ae', 'af', 'ah'},
    {'2ai', '2ak', '2am', 'ai', 'ak', 'am'},
    {'2ac', '2ad', '2ae', '2af', '2ah', '2ai', '2ak', '2am', 'ac', 'ad',
     'ae', 'af', 'ah', 'ai', 'ak', 'am'},
    {'2ac', '2ad', '2ae', '2af', '2ah', '2ai', '2ak', '2am', 'ac', 'ad',
     'ae', 'af', 'ah', 'ai', 'ak', 'am'},
    {'2ac', '2ad', '2ae', '2af', '2ah', '2ai', '2ak', '2am', '2aq', '2at',
     'ac', 'ad', 'ae', 'af', 'ah', 'ai', 'ak', 'am', 'aq', 'at'}]



