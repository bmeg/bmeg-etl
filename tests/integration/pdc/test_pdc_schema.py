from transform.pdc.pdc_schemas import PDCSchema


def test_all():
    pdc_schema = PDCSchema()

    query_names = list(pdc_schema.queries_dict().keys())
    assert len(query_names) > 0, 'Schema should have at least one query'
    assert len(pdc_schema.query_fields(query_names[0])) > 0, 'Query should return at least one field'

    entity_names = list(pdc_schema.entities_dict().keys())
    assert len(entity_names) > 0, 'Schema should have at least one entity'
    assert len(pdc_schema.entity_fields(entity_names[0])) > 0, 'Entity should return at least one field'

    studyExperimentalDesign = pdc_schema.make_query('studyExperimentalDesign')
    assert '$study_id' in studyExperimentalDesign, 'Should have parameter'
    assert '$study_submitter_id' in studyExperimentalDesign, 'Should have parameter'
