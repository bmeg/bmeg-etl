import orjson
from fhir.resources.patient import Patient
from fhir.resources.specimen import Specimen
from fhir.resources.reference import Reference
from fhir.resources.identifier import Identifier
from fhir.resources.codeableconcept import CodeableConcept


def build_patient(row):
    # TODO: replace and confirm BMEG id def
    patient_identifier = Identifier(**{"value": row['SUBJID'], "system": "".join(["https://gtexportal.org/", "subject"])})
    patient = Patient(**{"id": "-".join([row['SUBJID'], "gtex"]), "identifier": [patient_identifier], "birthDate": row['fhir_birthDate'], "gender": row['SEX_string']})
    return orjson.loads(patient.json())


def build_specimen(row):
    # TODO: replace and confirm BMEG id def
    cc = CodeableConcept(**{"coding": [{"system": "http://snomed.info/sct", "display": row['SMTSD'], "code": row['sctid']}]})
    subject_ref = Reference(**{"reference": "/".join(["Patient", "-".join([row['SUBJID'], "gtex"])])})
    specimen = Specimen(**{"id": "-".join([row['SAMPID'], "gtex"]), "type": cc, "subject": subject_ref})
    return orjson.loads(specimen.json())

