import gene_enricher


def test_simple():
    """ straightforward """
    tp53 = gene_enricher.get_gene('TP53')
    assert(tp53), 'Should exist'
    assert tp53 == [{'symbol': u'TP53', 'entrez_id': u'7157',
                    'ensembl_gene_id': u'ENSG00000141510'}]


def test_ambiguous():
    """ "ABC1" can point to both "ABCA1" and "HEATR6", """
    try:
        gene_enricher.get_gene('ABC1')
        assert 'Should have raised value error'
    except ValueError:
        pass


def test_ZUFSP():
    """ "ZUFSP" is a previous synonym """
    zufsp = gene_enricher.get_gene('ZUFSP')
    assert(zufsp), 'Should exist'


def test_FOOBAR():
    try:
        gene_enricher.get_gene('FOOBAR')
        assert 'Should have raised value error'
    except ValueError:
        pass
