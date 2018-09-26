from gripql import eq
"""
test gdc data in graph
V is a test fixture greated in conftest.py
"""


def test_project(V):
    q = (V.
         where(eq("_label", "Project")).
         where(eq("project_id", "TCGA-BRCA")))

    res = q.execute()
    assert len(res) == 1
    proj = res[0]
    assert proj.gid == "Project:TCGA-BRCA"
    assert proj.label == "Project"
    assert proj.data.project_id == "TCGA-BRCA"


def test_project_biosamples(V):
    q = (V.
         where(eq("_label", "Project")).
         where(eq("project_id", "TCGA-BRCA")).
         in_("InProject").
         in_("BiosampleFor"))

    res = q.execute()
    assert len(res) > 0
    first = res[0]
    print(first)


def test_biosample_gid(V):
    q = (V.
         where(eq("_gid", "Biosample:2ba23626-b0f4-449b-9e4f-a4163c1cf474")))
    res = q.execute()
    assert len(res) == 1


def test_project_aliquots(V):
    q = (V.
         where(eq("_label", "Project")).
         where(eq("project_id", "TCGA-BRCA")).
         in_("InProject").
         in_("BiosampleFor").
         in_("AliquotFor"))

    res = q.execute()
    assert len(res) > 0
    first = res[0]
    print(first)


def test_aliquot_gid(V):
    res = V.where(eq("_gid", "Aliquot:TCGA-AR-A0U1-01A-11W-A12T-09")).execute()
    assert len(res) == 1
