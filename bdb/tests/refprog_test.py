from nose.tools import eq_, ok_, raises

from bdb.refprog import (is_bdb_includable_refprog, except_refprog_warn,
                         filter_progs, last_used, one_of_the_two,
                         parse_refprog)


def test_is_bdb_includable_refprog_none():
    refprog = None
    result = is_bdb_includable_refprog(refprog)
    eq_(result, True)


def test_is_bdb_includable_refprog_cns():
    refprog = "CNS"
    result = is_bdb_includable_refprog(refprog)
    eq_(result, True)


def test_is_bdb_includable_refprog_xplor():
    refprog = "X-PLOR"
    result = is_bdb_includable_refprog(refprog)
    eq_(result, True)


def test_is_bdb_includable_refprog_buster():
    refprog = "BUSTER"
    result = is_bdb_includable_refprog(refprog)
    eq_(result, True)


def test_is_bdb_includable_refprog_refmac():
    refprog = "REFMAC"
    result = is_bdb_includable_refprog(refprog)
    eq_(result, True)


def test_is_bdb_includable_refprog_phenix():
    refprog = "PHENIX"
    result = is_bdb_includable_refprog(refprog)
    eq_(result, True)


def test_except_refprog_warn():
    pdb_id = "test"
    eq_(except_refprog_warn(pdb_id), None)


def test_filter_progs_empty():
    pin = []
    pv = []
    result = filter_progs(pin, pv)
    eq_(result, ([],[]))


@raises(AssertionError)
def test_filter_progs_type_error():
    pin = "program"
    pv = []
    filter_progs(pin, pv)

    pin = []
    pv = "version"
    filter_progs(pin, pv)


@raises(AssertionError)
def test_filter_progs_different_length():
    pin = []
    pv = ["version"]
    filter_progs(pin, pv)

    pin = ["program"]
    pv = []
    filter_progs(pin, pv)


def test_filter_progs_no_filter():
    pin = ["REFMAC", "SHELX"]
    pv = ["5.6", "L"]
    result = filter_progs(pin, pv)
    eq_(result, (pin, pv))


def test_filter_progs_filtered():
    pin = ["REFMAC", "COOT"]
    pv = ["5.7", "-"]
    result = filter_progs(pin, pv)
    expected = (["REFMAC"], ["5.7"])
    eq_(result, expected)

    pin = ["SOLVE/RESOLVE", "TOM/FRODO", "XTALVIEW"]
    pv = ["-", "-", "-"]
    result = filter_progs(pin, pv)
    expected = ([], [])
    eq_(result, expected)

    pop_progs = ["ARP", "ARP/WARP", "CHAIN", "COOT", "DM", "FRODO", "HKL-3000",
                 "LAFIRE", "MOLPROBITY", "O", "OOPS", "PIKSOL", "PRODRG",
                 "PROTEIN", "QUANTA", "SCWRL", "SFALL", "SOLVE/RESOLVE", "TOM",
                 "TOM/FRODO", "XFIT", "XPLEO", "XTALVIEW"]
    pv = ["-"]
    for p in pop_progs:
        result = filter_progs([p], pv)
        eq_(result, ([], []))



def test_last_used_empty():
    pin = []
    pv = []
    result = last_used(pin, pv)
    eq_(result, [])


@raises(AssertionError)
def test_last_used_type_error():
    pin = "program"
    pv = []
    last_used(pin, pv)

    pin = []
    pv = "version"
    last_used(pin, pv)


@raises(AssertionError)
def test_last_used_different_length():
    pin = []
    pv = ["version"]
    last_used(pin, pv)

    pin = ["program"]
    pv = []
    last_used(pin, pv)


def test_last_used_filtered():
    pin = ["REFMAC", "COOT"]
    pv = ["5.7", "-"]
    result = last_used(pin, pv)
    expected = ["REFMAC"]
    eq_(result, expected)

    pin = ["SOLVE/RESOLVE", "TOM/FRODO", "XTALVIEW"]
    pv = ["-", "-", "-"]
    result = last_used(pin, pv)
    expected = []
    eq_(result, expected)


def test_last_used_shelx():
    pin = ["REFMAC", "SHELX"]
    pv = ["5.8", "L"]
    result = last_used(pin, pv)
    expected = ["SHELX"]
    eq_(result, expected)

    pin = ["REFMAC", "SHELX"]
    pv = ["5.7", "H"]
    result = last_used(pin, pv)
    expected = ["SHELX"]
    eq_(result, expected)

    pin = ["REFMAC", "SHELX"]
    pv = ["5.6", "-"]
    result = last_used(pin, pv)
    expected = ["REFMAC"]
    eq_(result, expected)

    pin = ["REFMAC", "SHELX"]
    pv = ["5.5", "np"]
    result = last_used(pin, pv)
    expected = ["REFMAC"]
    eq_(result, expected)

    pin = ["SHELX"]
    pv = ["-"]
    result = last_used(pin, pv)
    eq_(result, pin)

    pin = ["SHELX"]
    pv = ["L"]
    result = last_used(pin, pv)
    eq_(result, pin)


def test_last_used_same_refprogs():
    # e.g. 1c4k
    pin = ["X-PLOR", "X-PLOR"]
    pv = ["3.1", "3.8"]
    result = last_used(pin, pv)
    expected = ["X-PLOR"]
    eq_(result, expected)


def test_last_used_corels():
    pin = ["REFMAC", "CORELS"]
    pv = ["5.8", "-"]
    result = last_used(pin, pv)
    expected = ["REFMAC"]
    eq_(result, expected)

    pin = ["CORELS"]
    pv = ["-"]
    result = last_used(pin, pv)
    eq_(result, pin)


def test_last_used_known_combis():
    # order doesn't matter
    pin = ["CNS", "CCP4"]
    pv = ["-", "-"]
    result = last_used(pin, pv)
    expected = ["CNS"]
    eq_(result, expected)

    known_combis = [("CCP4", "CNS"), ("CCP4", "X-PLOR"), ("CNS", "BUSTER"),
                    ("CNS", "PHENIX.REFINE"), ("CNS", "REFMAC"),
                    ("CNS", "TNT"), ("CNX", "BUSTER"),
                    ("CNX", "PHENIX.REFINE"), ("CNX", "REFMAC"),
                    ("CNX", "TNT"), ("EREF", "PROLSQ"), ("EREF", "X-PLOR"),
                    ("PROLSQ", "TNT"), ("PROLSQ", "REFMAC"),
                    ("REFMAC", "MAIN"), ("RESTRAIN", "X-PLOR"),
                    ("TWIN_LSQ", "CNS"), ("TNT", "REFMAC"),
                    ("X-PLOR", "BUSTER"), ("X-PLOR", "CNS"),
                    ("X-PLOR", "PHENIX.REFINE"), ("X-PLOR", "PROLSQ"),
                    ("X-PLOR", "REFMAC"), ("X-PLOR", "TNT"), ]
    for c in known_combis:
        result = last_used([c[0], c[1]], pv)
        expected = [c[1]]
        eq_(result, expected)



def test_last_used_unknown_combis():
    pin = ["X-PLOR", "CNS", "REFMAC", "PHENIX"]
    pv = ["3.85", "1.3", "5.8", "-"]
    result = last_used(pin, pv)
    expected = ["PHENIX", "REFMAC"]
    eq_(result, expected)


def test_last_used_sort():
    pin = ["B", "A"]
    pv = ["-", "-"]
    result = last_used(pin, pv)
    expected = ["A", "B"]
    eq_(result, expected)


def test_one_of_the_two():
    lst = ["A", "B"]
    loser = "A"
    winner = "B"
    result = one_of_the_two(lst, loser, winner)
    expected = ["B"]
    eq_(result, expected)


@raises(AssertionError)
def test_one_of_the_two_type_error():
    lst = "A"
    loser = "A"
    winner = "B"
    one_of_the_two(lst, loser, winner)


def test_one_of_the_two_fail():
    lst = ["A", "B"]
    loser = "A"
    winner = "C"
    result = one_of_the_two(lst, loser, winner)
    expected = ["A", "B"]
    eq_(result, expected)

    lst = ["A", "B"]
    loser = "C"
    winner = "B"
    result = one_of_the_two(lst, loser, winner)
    expected = ["A", "B"]
    eq_(result, expected)

    lst = ["A", "B"]
    loser = "B"
    winner = "B"
    result = one_of_the_two(lst, loser, winner)
    expected = ["A", "B"]
    eq_(result, expected)


def test_parse_refprog_null():
    refprog = "NULL"
    result = parse_refprog(refprog, "test")
    expected = ([None], [None], [None])
    eq_(result, expected)


def test_parse_refprog_predefined_exceptions():
    refprog = "X-PLOR 3.1 AND 3.85"
    result = parse_refprog(refprog, "test")
    expected =  (["X-PLOR 3.85"], ["X-PLOR"], ["3.85"])
    eq_(result, expected)

    refprog = "X-PLOR 3.1, 3.816"
    result = parse_refprog(refprog, "test")
    expected =  (["X-PLOR 3.816"], ["X-PLOR"], ["3.816"])
    eq_(result, expected)

    refprog = "CNS 1.1 & 1.3"
    result = parse_refprog(refprog, "test")
    expected =  (["CNS 1.3"], ["CNS"], ["1.3"])
    eq_(result, expected)

    refprog = "CNS 0.4, O, OOPS"
    result = parse_refprog(refprog, "test")
    expected =  (["CNS 0.4", "O", "OOPS"],
                 ["CNS", "O", "OOPS"],
                 ["0.4", "-", "-"])
    eq_(result, expected)

    refprog = "CNS 0.1-0.4"
    result = parse_refprog(refprog, "test")
    expected =  (["CNS 0.4"], ["CNS"], ["0.4"])
    eq_(result, expected)

    refprog = "CNS 0.9,1.0,1.1"
    result = parse_refprog(refprog, "test")
    expected =  (["CNS 1.1"], ["CNS"], ["1.1"])
    eq_(result, expected)

    refprog = "CNS 1.3 WITH DEN REFINEMENT"
    result = parse_refprog(refprog, "test")
    expected =  (["CNS 1.3"], ["CNS"], ["1.3"])
    eq_(result, expected)

    refprog = "CNS 1.2 (USING XTAL_TWIN UTILITIES)"
    result = parse_refprog(refprog, "test")
    expected =  (["CNS 1.2"], ["CNS"], ["1.2"])
    eq_(result, expected)

    refprog = "PHENIX.REFINE_REFMAC 5.5.0070"
    result = parse_refprog(refprog, "test")
    expected =  (["PHENIX.REFINE", "REFMAC 5.5.0070"],
                 ["PHENIX.REFINE", "REFMAC"],
                 ["-", "5.5.0070"])
    eq_(result, expected)

    refprog = "PHENIX (CCI APPS 2007_04_06_1210)"
    result = parse_refprog(refprog, "test")
    expected =  (["PHENIX (PHENIX.REFINE: 2007_04_06_1210)"],
                 ["PHENIX.REFINE"],
                 ["2007_04_06_1210"])
    eq_(result, expected)

    refprog = "PHENIX VERSION 1.8_1069 (PHENIX.REFINE)"
    result = parse_refprog(refprog, "test")
    expected =  (["PHENIX (PHENIX.REFINE: 1.8_1069)"], ["PHENIX.REFINE"],
                 ["1.8_1069"])
    eq_(result, expected)

    refprog = "PHENIX 1.6.2_432 - REFINE"
    result = parse_refprog(refprog, "test")
    expected =  (["PHENIX (PHENIX.REFINE: 1.6.2_432)"], ["PHENIX.REFINE"],
                 ["1.6.2_432"])
    eq_(result, expected)

    refprog = "PHENIX REFINE"
    result = parse_refprog(refprog, "test")
    expected = (["PHENIX.REFINE"], ["PHENIX.REFINE"], ["-"])
    eq_(result, expected)

    refprog = "PHENIX, REFINE"
    result = parse_refprog(refprog, "test")
    expected = (["PHENIX.REFINE"], ["PHENIX.REFINE"], ["-"])
    eq_(result, expected)

    refprog = "PHENIX AUTOREFINE"
    result = parse_refprog(refprog, "test")
    expected = (["PHENIX.REFINE"], ["PHENIX.REFINE"], ["-"])
    eq_(result, expected)

    refprog = "REFMAC 5.1.24/TLS"
    result = parse_refprog(refprog, "test")
    expected =  (["REFMAC 5.1.24"], ["REFMAC"], ["5.1.24"])
    eq_(result, expected)

    refprog = "REFMAC 5.2.0005 24/04/2001"
    result = parse_refprog(refprog, "test")
    expected =  (["REFMAC 5.2.0005"], ["REFMAC"], ["5.2.0005"])
    eq_(result, expected)

    refprog = "REFMAC 5.2.0019 24/04/2001"
    result = parse_refprog(refprog, "test")
    expected =  (["REFMAC 5.2.0019"], ["REFMAC"], ["5.2.0019"])
    eq_(result, expected)

    refprog = "REFMAC X-PLOR 3.843"
    result = parse_refprog(refprog, "test")
    expected =  (["REFMAC", "X-PLOR 3.843"],
                 ["REFMAC", "X-PLOR"],
                 ["-", "3.843"])
    eq_(result, expected)

    refprog = "REFMAC5 5.2.0019"
    result = parse_refprog(refprog, "test")
    expected =  (["REFMAC 5.2.0019"], ["REFMAC"], ["5.2.0019"])
    eq_(result, expected)

    refprog = "REFMAC 5.5.0109 (AND PHENIX)"
    result = parse_refprog(refprog, "test")
    expected =  (["REFMAC 5.5.0109", "PHENIX.REFINE"],
                 ["REFMAC", "PHENIX.REFINE"],
                 ["5.5.0109", "-"])
    eq_(result, expected)

    refprog = "BUSTER, BETA VERSION"
    result = parse_refprog(refprog, "test")
    expected =  (["BUSTER BETA"], ["BUSTER"], ["BETA"])
    eq_(result, expected)

    refprog = "TNT BUSTER/TNT"
    result = parse_refprog(refprog, "test")
    expected =  (["TNT", "BUSTER/TNT"], ["TNT", "BUSTER"], ["-", "-"])
    eq_(result, expected)

    refprog = "O, VERSION 9.0.7"
    result = parse_refprog(refprog, "test")
    expected =  (["O 9.0.7"], ["O"], ["9.0.7"])
    eq_(result, expected)


def test_parse_refprog_nosplit():
    refprog = "PHENIX (PHENIX.REFINE: DEV_1549+SVN)      "
    result = parse_refprog(refprog, "4ow3")
    expected =  (["PHENIX (PHENIX.REFINE: DEV_1549+SVN)"], ["PHENIX.REFINE"],
                 ["np"])
    eq_(result, expected)

    refprog = "BUSTER/TNT"
    result = parse_refprog(refprog, "test")
    expected =  (["BUSTER/TNT"], ["BUSTER"], ["-"])
    eq_(result, expected)

    refprog = "BUSTER/RESOLVE"
    result = parse_refprog(refprog, "test")
    expected =  (["BUSTER/RESOLVE"], ["BUSTER"], ["np"])
    eq_(result, expected)

    refprog = "ARP/WARP"
    result = parse_refprog(refprog, "test")
    expected =  (["ARP/WARP"], ["ARP/WARP"], ["-"])
    eq_(result, expected)

    refprog = "WARP/ARP"
    result = parse_refprog(refprog, "test")
    expected =  (["WARP/ARP"], ["OTHER"], ["np"])
    eq_(result, expected)


def test_parse_refprog_split():
    refprog = "TNT/BUSTER"
    result = parse_refprog(refprog, "test")
    expected =  (["TNT", "BUSTER"], ["TNT", "BUSTER"], ["-", "-"])
    eq_(result, expected)

    refprog = "TNT, BUSTER"
    result = parse_refprog(refprog, "test")
    expected =  (["TNT", "BUSTER"], ["TNT", "BUSTER"], ["-", "-"])
    eq_(result, expected)

    refprog = "TNT&BUSTER"
    result = parse_refprog(refprog, "test")
    expected =  (["TNT", "BUSTER"], ["TNT", "BUSTER"], ["-", "-"])
    eq_(result, expected)

    refprog = "TNT ; BUSTER"
    result = parse_refprog(refprog, "test")
    expected =  (["TNT", "BUSTER"], ["TNT", "BUSTER"], ["-", "-"])
    eq_(result, expected)

    refprog = "TNT +, BUSTER"
    result = parse_refprog(refprog, "test")
    expected =  (["TNT", "BUSTER"], ["TNT", "BUSTER"], ["-", "-"])
    eq_(result, expected)

    refprog = "TNT AND BUSTER,"
    result = parse_refprog(refprog, "test")
    expected =  (["TNT", "BUSTER"], ["TNT", "BUSTER"], ["-", "-"])
    eq_(result, expected)

    refprog = "/REFMAC5.7/////COOT//"
    result = parse_refprog(refprog, "test")
    expected =  (["REFMAC5.7", "COOT"], ["REFMAC", "COOT"], ["5.7", "-"])
    eq_(result, expected)


def test_parse_refprog_no_versions():
    refprogs = ["ARP", "BILDER", "CCP4", "CEDAR", "CHAIN", "CORELS", "DM",
                "DYNAMIX", "EREF", "FFX", "FMLS/VP", "FRODO", "HIPHOP",
                "IMPLOR", "LAFIRE", "LALS", "MAIN", "MOLPROBITY", "MOLLY",
                "MOPRO", "NCNS", "NCNS-TINKER", "NMREF", "PIKSOL", "PHASER",
                "PMB", "POLYVISION", "PRIMEX", "PRODRG", "PROTEIN", "QUANTA",
                "RESTRAIN", "SCWRL", "SHARP", "SFALL", "SOLVE",
                "SOLVE/RESOLVE", "TIBBITTS", "TOM", "TOM/FRODO", "XFIT",
                "XPLEO", "XTALVIEW"]
    for refprog in refprogs:
        result = parse_refprog(refprog, "test")
        expected = ([refprog], [refprog], ["-"])
        eq_(result, expected)


def test_parse_refprog_no_versions_with_version():
    refprog = "CCP4 VERSION"
    result = parse_refprog(refprog, "test")
    expected = ([refprog], ["CCP4"], ["np"])
    eq_(result, expected)


def test_parse_refprog_refmac():
    refprog = "refamc"
    result = parse_refprog(refprog, "test")
    expected = (["REFAMC"], ["REFMAC"], ["-"])
    eq_(result, expected)

    refprog = "REFMEC"
    result = parse_refprog(refprog, "test")
    expected = (["REFMEC"], ["REFMAC"], ["-"])
    eq_(result, expected)

    # Theoretically possible...
    refprog = "CCCCCC"
    result = parse_refprog(refprog, "test")
    expected = (["CCCCCC"], ["REFMAC"], ["-"])
    eq_(result, expected)

    refprog = "REFMAC"
    result = parse_refprog(refprog, "test")
    expected = (["REFMAC"], ["REFMAC"], ["-"])
    eq_(result, expected)

    refprog = "REFMAC 5.7"
    result = parse_refprog(refprog, "test")
    expected = (["REFMAC 5.7"], ["REFMAC"], ["5.7"])
    eq_(result, expected)

    refprog = "REFMAC5.7"
    result = parse_refprog(refprog, "test")
    expected = (["REFMAC5.7"], ["REFMAC"], ["5.7"])
    eq_(result, expected)

    refprog = "REFMAC V 5.7"
    result = parse_refprog(refprog, "test")
    expected = (["REFMAC V 5.7"], ["REFMAC"], ["5.7"])
    eq_(result, expected)

    refprog = "REFMAC V. 5.7"
    result = parse_refprog(refprog, "test")
    expected = (["REFMAC V. 5.7"], ["REFMAC"], ["np"])
    eq_(result, expected)

    refprog = "REFMAC NONE"
    result = parse_refprog(refprog, "test")
    expected = (["REFMAC NONE"], ["REFMAC"], ["NONE"])
    eq_(result, expected)

    refprog = "REFMAC                   5.8                  "
    result = parse_refprog(refprog, "test")
    expected = (["REFMAC                   5.8"], ["REFMAC"], ["np"])
    eq_(result, expected)

    refprog = "REFMAC 5.8                  "
    result = parse_refprog(refprog, "test")
    expected = (["REFMAC 5.8"], ["REFMAC"], ["5.8"])
    eq_(result, expected)


def test_parse_refprog_cns():
    refprog = "CNS"
    result = parse_refprog(refprog, "test")
    expected = (["CNS"], ["CNS"], ["-"])
    eq_(result, expected)

    refprog = "CNS          "
    result = parse_refprog(refprog, "test")
    expected = (["CNS"], ["CNS"], ["-"])
    eq_(result, expected)

    refprog = "CNS (2)"
    result = parse_refprog(refprog, "test")
    expected = (["CNS (2)"], ["CNS"], ["2"])
    eq_(result, expected)

    refprog = "CNS ( 2 )"
    result = parse_refprog(refprog, "test")
    expected = (["CNS ( 2 )"], ["CNS"], ["np"])
    eq_(result, expected)

    refprog = "CNS 2"
    result = parse_refprog(refprog, "test")
    expected = (["CNS 2"], ["CNS"], ["2"])
    eq_(result, expected)

    refprog = "CNS2.0"
    result = parse_refprog(refprog, "test")
    expected = (["CNS2.0"], ["CNS"], ["2.0"])
    eq_(result, expected)

    refprog = "CNS V. ABC123.456"
    result = parse_refprog(refprog, "test")
    expected = (["CNS V. ABC123.456"], ["CNS"], ["ABC123.456"])
    eq_(result, expected)

    refprog = "CNS-3"
    result = parse_refprog(refprog, "test")
    expected = (["CNS-3"], ["CNS"], ["3"])
    eq_(result, expected)

    # Theoretically possible...
    refprog = "CNS-V.3"
    result = parse_refprog(refprog, "test")
    expected = (["CNS-V.3"], ["CNS"], ["V.3"])
    eq_(result, expected)

    # Theoretically possible...
    refprog = "CNS (V.3)"
    result = parse_refprog(refprog, "test")
    expected = (["CNS (V.3)"], ["CNS"], ["V.3"])
    eq_(result, expected)

    # Theoretically possible...
    refprog = "CNS-"
    result = parse_refprog(refprog, "test")
    expected = (["CNS-"], ["CNS"], ["np"])
    eq_(result, expected)


def test_parse_refprog_twinlsq():
    refprog = "TWIN_LSQ"
    result = parse_refprog(refprog, "test")
    expected = (["TWIN_LSQ"], ["CNS"], ["TWIN_LSQ"])
    eq_(result, expected)

    refprog = "TWIN_LSQ 4"
    result = parse_refprog(refprog, "test")
    expected = (["TWIN_LSQ 4"], ["CNS"], ["np"])
    eq_(result, expected)


def test_parse_refprog_cnx():
    refprog = "CNX"
    result = parse_refprog(refprog, "test")
    expected = (["CNX"], ["CNX"], ["-"])
    eq_(result, expected)

    refprog = "CNX          "
    result = parse_refprog(refprog, "test")
    expected = (["CNX"], ["CNX"], ["-"])
    eq_(result, expected)

    refprog = "CNX (ACCELRYS)"
    result = parse_refprog(refprog, "test")
    expected = (["CNX (ACCELRYS)"], ["CNX"], ["-"])
    eq_(result, expected)

    refprog = "CNX V0"
    result = parse_refprog(refprog, "test")
    expected = (["CNX V0"], ["CNX"], ["np"])
    eq_(result, expected)

    refprog = "CNX 0"
    result = parse_refprog(refprog, "test")
    expected = (["CNX 0"], ["CNX"], ["0"])
    eq_(result, expected)

    refprog = "CNX 0.9-9"
    result = parse_refprog(refprog, "test")
    expected = (["CNX 0.9-9"], ["CNX"], ["0.9-9"])
    eq_(result, expected)

    refprog = "CNX0.9-9"
    result = parse_refprog(refprog, "test")
    expected = (["CNX0.9-9"], ["CNX"], ["0.9-9"])
    eq_(result, expected)


def test_parse_refprog_xplor():
    refprog = "X-PLOR"
    result = parse_refprog(refprog, "test")
    expected = (["X-PLOR"], ["X-PLOR"], ["-"])
    eq_(result, expected)

    refprog = "XPLOR"
    result = parse_refprog(refprog, "test")
    expected = (["XPLOR"], ["X-PLOR"], ["-"])
    eq_(result, expected)

    refprog = "X-PLOR (ONLINE)"
    result = parse_refprog(refprog, "test")
    expected = (["X-PLOR (ONLINE)"], ["X-PLOR"], ["-"])
    eq_(result, expected)

    refprog = "X-PLOR V5.4321"
    result = parse_refprog(refprog, "test")
    expected = (["X-PLOR V5.4321"], ["X-PLOR"], ["V5.4321"])
    eq_(result, expected)

    refprog = "X-PLOR V 5"
    result = parse_refprog(refprog, "test")
    expected = (["X-PLOR V 5"], ["X-PLOR"], ["np"])
    eq_(result, expected)


def test_parse_refprog_phenix_ensemble_refinement():
    refprog = "PHENIX.ENSEMBLE_REFINEMENT"
    result = parse_refprog(refprog, "test")
    expected = (["PHENIX.ENSEMBLE_REFINEMENT"], ["PHENIX.ENSEMBLE_REFINEMENT"],
                ["np"])
    eq_(result, expected)

    refprog = "PHENIX (PHENIX.ENSEMBLE_REFINEMENT: dev-9999)"
    result = parse_refprog(refprog, "test")
    expected = (["PHENIX (PHENIX.ENSEMBLE_REFINEMENT: DEV-9999)"],
                ["PHENIX.ENSEMBLE_REFINEMENT"], ["DEV-9999"])
    eq_(result, expected)


def test_parse_refprog_phenix():
    refprog = "PHENIX"
    result = parse_refprog(refprog, "test")
    expected = (["PHENIX"], ["PHENIX.REFINE"], ["-"])
    eq_(result, expected)

    refprog = "PHENIX.REFINE"
    result = parse_refprog(refprog, "test")
    expected = (["PHENIX.REFINE"], ["PHENIX.REFINE"], ["-"])
    eq_(result, expected)

    refprog = "PHENIX.REFINEMENT"
    result = parse_refprog(refprog, "test")
    expected = (["PHENIX.REFINEMENT"], ["PHENIX.REFINE"], ["-"])
    eq_(result, expected)

    refprog = "PHENIX (PHENIX.REFINE)"
    result = parse_refprog(refprog, "test")
    expected = (["PHENIX (PHENIX.REFINE)"], ["PHENIX.REFINE"], ["-"])
    eq_(result, expected)

    refprog = "PHENIX (PHENIX)"
    result = parse_refprog(refprog, "test")
    expected = (["PHENIX (PHENIX)"], ["PHENIX.REFINE"], ["-"])
    eq_(result, expected)

    refprog = "PHENIX.REFINE (PHENIX.REFINE)"
    result = parse_refprog(refprog, "test")
    expected = (["PHENIX.REFINE (PHENIX.REFINE)"], ["PHENIX.REFINE"], ["-"])
    eq_(result, expected)

    refprog = "PHENIX.REFINEMENT (PHENIX.REFINE)"
    result = parse_refprog(refprog, "test")
    expected = (["PHENIX.REFINEMENT (PHENIX.REFINE)"], ["PHENIX.REFINE"],
                ["-"])
    eq_(result, expected)

    refprog = "PHENIX.REFINEMENT (PHENIX)"
    result = parse_refprog(refprog, "test")
    expected = (["PHENIX.REFINEMENT (PHENIX)"], ["PHENIX.REFINE"], ["-"])
    eq_(result, expected)

    refprog = "PHENIX 1.8"
    result = parse_refprog(refprog, "test")
    expected = (["PHENIX 1.8"], ["PHENIX.REFINE"], ["1.8"])
    eq_(result, expected)

    refprog = "PHENIX.REFINE 1.8AB"
    result = parse_refprog(refprog, "test")
    expected = (["PHENIX.REFINE 1.8AB"], ["PHENIX.REFINE"], ["np"])
    eq_(result, expected)

    refprog = "PHENIX.REFINE: 1.8AB"
    result = parse_refprog(refprog, "test")
    expected = (["PHENIX.REFINE: 1.8AB"], ["PHENIX.REFINE"], ["1.8AB"])
    eq_(result, expected)

    refprog = "PHENIX (1.8_1060)"
    result = parse_refprog(refprog, "test")
    expected = (["PHENIX (1.8_1060)"], ["PHENIX.REFINE"], ["1.8_1060"])
    eq_(result, expected)

    refprog = "PHENIX (PHENIX.REFINE: 1.8)"
    result = parse_refprog(refprog, "test")
    expected = (["PHENIX (PHENIX.REFINE: 1.8)"], ["PHENIX.REFINE"], ["1.8"])
    eq_(result, expected)

    refprog = "AND PHENIX.REFINE"
    result = parse_refprog(refprog, "test")
    expected = (["PHENIX.REFINE"], ["PHENIX.REFINE"], ["-"])
    eq_(result, expected)


def test_parse_refprog_buster():
    refprog = "BUSTER"
    result = parse_refprog(refprog, "test")
    expected = (["BUSTER"], ["BUSTER"], ["-"])
    eq_(result, expected)

    refprog = "AUTOBUSTER 4"
    result = parse_refprog(refprog, "test")
    expected = (["AUTOBUSTER 4"], ["BUSTER"], ["4"])
    eq_(result, expected)

    refprog = "BUSTER 2.10.0"
    result = parse_refprog(refprog, "test")
    expected = (["BUSTER 2.10.0"], ["BUSTER"], ["2.10.0"])
    eq_(result, expected)

    refprog = "BUSTER-TNT 2.11.2"
    result = parse_refprog(refprog, "test")
    expected = (["BUSTER-TNT 2.11.2"], ["BUSTER"], ["2.11.2"])
    eq_(result, expected)

    refprog = "BUSTER-TNT 2.X"
    result = parse_refprog(refprog, "test")
    expected = (["BUSTER-TNT 2.X"], ["BUSTER"], ["2.X"])
    eq_(result, expected)

    refprog = "BUSTER-TNT V. 1.1.1"
    result = parse_refprog(refprog, "test")
    expected = (["BUSTER-TNT V. 1.1.1"], ["BUSTER"], ["1.1.1"])
    eq_(result, expected)

    refprog = "BUSTER-TNT BUSTER 2.8.0"
    result = parse_refprog(refprog, "test")
    expected = (["BUSTER-TNT BUSTER 2.8.0"], ["BUSTER"], ["2.8.0"])
    eq_(result, expected)

    refprog = "BUSTER TNT"
    result = parse_refprog(refprog, "test")
    expected = (["BUSTER TNT"], ["BUSTER"], ["-"])
    eq_(result, expected)

    refprog = "BUSTER/TNT"
    result = parse_refprog(refprog, "test")
    expected = (["BUSTER/TNT"], ["BUSTER"], ["-"])
    eq_(result, expected)

    refprog = "AUTOBUSTER"
    result = parse_refprog(refprog, "test")
    expected = (["AUTOBUSTER"], ["BUSTER"], ["-"])
    eq_(result, expected)


def test_parse_refprog_tnt():
    refprog = "TNT"
    result = parse_refprog(refprog, "test")
    expected = (["TNT"], ["TNT"], ["-"])
    eq_(result, expected)

    refprog = "TNT V. 5-E"
    result = parse_refprog(refprog, "test")
    expected = (["TNT V. 5-E"], ["TNT"], ["5-E"])
    eq_(result, expected)

    refprog = "TNT 5.6.1"
    result = parse_refprog(refprog, "test")
    expected = (["TNT 5.6.1"], ["TNT"], ["5.6.1"])
    eq_(result, expected)

    refprog = "TNT 5F"
    result = parse_refprog(refprog, "test")
    expected = (["TNT 5F"], ["TNT"], ["5F"])
    eq_(result, expected)

    refprog = "TNT V. 5-F PRERELEASE"
    result = parse_refprog(refprog, "test")
    expected = (["TNT V. 5-F PRERELEASE"], ["TNT"], ["5-F PRERELEASE"])
    eq_(result, expected)

    refprog = "TNT 5-F PRERELEASE"
    result = parse_refprog(refprog, "test")
    expected = (["TNT 5-F PRERELEASE"], ["TNT"], ["5-F PRERELEASE"])
    eq_(result, expected)


def test_parse_refprog_shelx():
    refprog = "SHELX"
    result = parse_refprog(refprog, "test")
    expected = (["SHELX"], ["SHELX"], ["-"])
    eq_(result, expected)

    refprog = "SHELX97"
    result = parse_refprog(refprog, "test")
    expected = (["SHELX97"], ["SHELX"], ["97"])
    eq_(result, expected)

    refprog = "SHELXH"
    result = parse_refprog(refprog, "test")
    expected = (["SHELXH"], ["SHELX"], ["H"])
    eq_(result, expected)

    refprog = "SHELXS"
    result = parse_refprog(refprog, "test")
    expected = (["SHELXS"], ["SHELX"], ["S"])
    eq_(result, expected)

    refprog = "SHELXL"
    result = parse_refprog(refprog, "test")
    expected = (["SHELXL"], ["SHELX"], ["L"])
    eq_(result, expected)

    refprog = "SHELX-97"
    result = parse_refprog(refprog, "test")
    expected = (["SHELX-97"], ["SHELX"], ["97"])
    eq_(result, expected)

    refprog = "SHELXL-97"
    result = parse_refprog(refprog, "test")
    expected = (["SHELXL-97"], ["SHELX"], ["L-97"])
    eq_(result, expected)

    refprog = "SHELXH-97"
    result = parse_refprog(refprog, "test")
    expected = (["SHELXH-97"], ["SHELX"], ["H-97"])
    eq_(result, expected)

    refprog = "SHELX-97-1"
    result = parse_refprog(refprog, "test")
    expected = (["SHELX-97-1"], ["SHELX"], ["97-1"])
    eq_(result, expected)

    refprog = "SHELX 97-1"
    result = parse_refprog(refprog, "test")
    expected = (["SHELX 97-1"], ["SHELX"], ["97-1"])
    eq_(result, expected)

    refprog = "SHELXL 97"
    result = parse_refprog(refprog, "test")
    expected = (["SHELXL 97"], ["SHELX"], ["L 97"])
    eq_(result, expected)

    refprog = "SHELXL97"
    result = parse_refprog(refprog, "test")
    expected = (["SHELXL97"], ["SHELX"], ["L97"])
    eq_(result, expected)

    refprog = "SHELX-L"
    result = parse_refprog(refprog, "test")
    expected = (["SHELX-L"], ["SHELX"], ["L"])
    eq_(result, expected)


def test_parse_refprog_coot():
    refprog = "COOT"
    result = parse_refprog(refprog, "test")
    expected = (["COOT"], ["COOT"], ["-"])
    eq_(result, expected)

    refprog = "COOT 0.1-3-PRE-1"
    result = parse_refprog(refprog, "test")
    expected = (["COOT 0.1-3-PRE-1"], ["COOT"], ["0.1-3-PRE-1"])
    eq_(result, expected)

    refprog = "COOT 0.0.33"
    result = parse_refprog(refprog, "test")
    expected = (["COOT 0.0.33"], ["COOT"], ["0.0.33"])
    eq_(result, expected)

    refprog = "COOT V. 0.0.33"
    result = parse_refprog(refprog, "test")
    expected = (["COOT V. 0.0.33"], ["COOT"], ["0.0.33"])
    eq_(result, expected)

    refprog = "COOT $"
    result = parse_refprog(refprog, "test")
    expected = (["COOT $"], ["COOT"], ["np"])
    eq_(result, expected)


def test_parse_refprog_o():
    refprog = "O"
    result = parse_refprog(refprog, "test")
    expected = (["O"], ["O"], ["-"])
    eq_(result, expected)

    refprog = "O 9.0.7"
    result = parse_refprog(refprog, "test")
    expected = (["O 9.0.7"], ["O"], ["np"])
    eq_(result, expected)

    refprog = "O9.0.7"
    result = parse_refprog(refprog, "test")
    expected = (["O9.0.7"], ["O"], ["9.0.7"])
    eq_(result, expected)


def test_parse_refprog_profft():
    refprog = "PROFFT"
    result = parse_refprog(refprog, "test")
    expected = (["PROFFT"], ["PROLSQ"], ["PROFFT"])
    eq_(result, expected)

    refprog = "PROFFT (MODIFIED BY Z.OTWINOWSKI)"
    result = parse_refprog(refprog, "test")
    expected = ([refprog], ["PROLSQ"], ["PROFFT MODIFIED BY Z.OTWINOWSKI"])
    eq_(result, expected)

    refprog = "PROFFT anything else"
    result = parse_refprog(refprog, "test")
    expected = (["PROFFT ANYTHING ELSE"], ["PROLSQ"], ["np"])
    eq_(result, expected)


def test_parse_refprog_prolsq():
    refprog = "PROLSQ"
    result = parse_refprog(refprog, "test")
    expected = (["PROLSQ"], ["PROLSQ"], ["PROLSQ"])
    eq_(result, expected)

    refprog = "PROLSQ (MODIFIED BY G.J.QUIGLEY)"
    result = parse_refprog(refprog, "test")
    expected = ([refprog], ["PROLSQ"], ["PROLSQ MODIFIED BY G.J.QUIGLEY"])
    eq_(result, expected)

    refprog = "PROLSQ anything else"
    result = parse_refprog(refprog, "test")
    expected = (["PROLSQ ANYTHING ELSE"], ["PROLSQ"], ["np"])
    eq_(result, expected)


def test_parse_refprog_protin():
    refprog = "PROTIN"
    result = parse_refprog(refprog, "test")
    expected = (["PROTIN"], ["PROLSQ"], ["PROLSQ"])
    eq_(result, expected)

    refprog = "PROTIN anything else"
    result = parse_refprog(refprog, "test")
    expected = (["PROTIN ANYTHING ELSE"], ["PROLSQ"], ["np"])
    eq_(result, expected)


def test_parse_refprog_nuclsq():
    refprog = "NUCLSQ"
    result = parse_refprog(refprog, "test")
    expected = (["NUCLSQ"], ["PROLSQ"], ["NUCLSQ"])
    eq_(result, expected)

    refprog = "NUCLSQ (MODIFIED BY G.J.QUIGLEY)"
    result = parse_refprog(refprog, "test")
    expected = ([refprog], ["PROLSQ"], ["NUCLSQ MODIFIED BY G.J.QUIGLEY"])
    eq_(result, expected)

    refprog = "NUCLSQ anything else"
    result = parse_refprog(refprog, "test")
    expected = (["NUCLSQ ANYTHING ELSE"], ["PROLSQ"], ["np"])
    eq_(result, expected)


def test_parse_refprog_nuclin():
    refprog = "NUCLIN"
    result = parse_refprog(refprog, "test")
    expected = (["NUCLIN"], ["PROLSQ"], ["NUCLSQ"])
    eq_(result, expected)

    refprog = "PROTIN anything else"
    result = parse_refprog(refprog, "test")
    expected = (["PROTIN ANYTHING ELSE"], ["PROLSQ"], ["np"])
    eq_(result, expected)


def test_parse_refprog_gprlsa():
    refprog = "GPRLSA"
    result = parse_refprog(refprog, "test")
    expected = (["GPRLSA"], ["PROLSQ"], ["GPRLSA"])
    eq_(result, expected)

    refprog = "GPRLSA anything else"
    result = parse_refprog(refprog, "test")
    expected = (["GPRLSA ANYTHING ELSE"], ["PROLSQ"], ["np"])
    eq_(result, expected)


def test_parse_refprog_deriv():
    refprog = "DERIV"
    result = parse_refprog(refprog, "test")
    expected = (["DERIV"], ["PROLSQ"], ["DERIV"])
    eq_(result, expected)

    refprog = "DERIV anything else"
    result = parse_refprog(refprog, "test")
    expected = (["DERIV ANYTHING ELSE"], ["PROLSQ"], ["np"])
    eq_(result, expected)


def test_parse_refprog_cerius():
    refprog = "CERIUS"
    result = parse_refprog(refprog, "test")
    expected = (["CERIUS"], ["CERIUS"], ["-"])
    eq_(result, expected)

    refprog = "CERIUS 2"
    result = parse_refprog(refprog, "test")
    expected = (["CERIUS 2"], ["CERIUS"], ["2"])
    eq_(result, expected)

    refprog = "CERIUS2"
    result = parse_refprog(refprog, "test")
    expected = (["CERIUS2"], ["CERIUS"], ["2"])
    eq_(result, expected)

    refprog = "CERIUS anything else"
    result = parse_refprog(refprog, "test")
    expected = (["CERIUS ANYTHING ELSE"], ["CERIUS"], ["np"])
    eq_(result, expected)


def test_parse_refprog_hkl3000():
    refprog = "HKL-3000"
    result = parse_refprog(refprog, "test")
    expected = (["HKL-3000"], ["HKL-3000"], ["-"])
    eq_(result, expected)

    refprog = "HKL3000"
    result = parse_refprog(refprog, "test")
    expected = (["HKL3000"], ["HKL-3000"], ["-"])
    eq_(result, expected)

    refprog = "HKL-3000 anything else"
    result = parse_refprog(refprog, "test")
    expected = (["HKL-3000 ANYTHING ELSE"], ["HKL-3000"], ["np"])
    eq_(result, expected)


def test_parse_refprog_gromos():
    refprog = "GROMOS"
    result = parse_refprog(refprog, "test")
    expected = (["GROMOS"], ["GROMOS"], ["-"])
    eq_(result, expected)

    refprog = "GROMOS87"
    result = parse_refprog(refprog, "test")
    expected = (["GROMOS87"], ["GROMOS"], ["87"])
    eq_(result, expected)

    refprog = "GROMOS anything else"
    result = parse_refprog(refprog, "test")
    expected = (["GROMOS ANYTHING ELSE"], ["GROMOS"], ["np"])
    eq_(result, expected)


def test_parse_refprog_other():
    refprogs = ["AGARWAL FAST-FOURIER TRANSFORM LEAST-SQUARES",
                "HENDRICKSON-KONNERT LEAST-SQUARES REFINEMENT",
                "* OF T. A. JONES.  THE R VALUE IS 0.180.",
                "RESTRAINED RECIPROCAL-SPACE LEAST-SQUARES",
                "SIMULATED ANNEALING METHOD",
                "CONSTRAINED RECIPROCAL-SPACE LEAST-SQUARES",
                "FAST-FOURIER LEAST-SQUARES REFINEMENT",
                "JACK-LEVITT",
                "REAL-SPACE REFINEMENT",
                "SOME OTHER REFPROG",]
    for p in refprogs:
        result = parse_refprog(p, "test")
        expected = ([p], ["OTHER"], ["np"])
