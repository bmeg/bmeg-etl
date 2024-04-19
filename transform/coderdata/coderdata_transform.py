#!/usr/bin/env python

import pandas as pd
import orjson
import uuid
from fhir.resources.substance import Substance
from fhir.resources.reference import Reference
from fhir.resources.identifier import Identifier
from fhir.resources.codeablereference import CodeableReference
from fhir.resources.codeableconcept import CodeableConcept
from fhir.resources.substancedefinition import SubstanceDefinition, SubstanceDefinitionStructure, \
    SubstanceDefinitionStructureRepresentation
from bmegDrugResponse import DrugResponse
from bmegPatient import Patient


def build_substance(row):
    smile_rep = row['canSMILES']
    InChI_rep = row['InChIKey']
    codeable_concept_smile = CodeableConcept(
        **{"coding": [{"system": "http://hl7.org/fhir/substance-representation-format",
                       "code": "SMILES",
                       "display": "SMILES"}]})
    sdfr_smile = SubstanceDefinitionStructureRepresentation(
        **{"representation": smile_rep, "format": codeable_concept_smile})

    codeable_concept_InChI = CodeableConcept(
        **{"coding": [{"system": "http://hl7.org/fhir/substance-representation-format",
                       "code": "InChI",
                       "display": "InChI"}]})
    sdfr_InChI = SubstanceDefinitionStructureRepresentation(
        **{"representation": InChI_rep, "format": codeable_concept_InChI})

    sds = SubstanceDefinitionStructure(**{"representation": [sdfr_smile, sdfr_InChI]})
    id = row['pubchem_id']
    substanceDef_ident = Identifier(**{"value": id, "system": "https://pubchem.ncbi.nlm.nih.gov/"})
    sd = SubstanceDefinition(**{"id": id, "identifier": [substanceDef_ident], "structure": sds})

    cr = CodeableReference(**{"reference": Reference(**{"reference": "/".join(["SubstanceDefinition", id])})})
    substance = Substance(**{"instance": True, "code": cr})
    return orjson.loads(substance.json()) 


def build_response(row):
    id = row['id']
    acc = row['acc']
    auc = row['auc']
    dss1 = row['dss']
    ec50 = row['ec50']
    einf = row['einf']
    hs = row['hs']
    ic50 = row['ic50']

    substance_row = row[['pubchem_id', 'canSMILES', 'InChI']]
    substance = build_substance(substance_row)

    # drugResponse BMEG BaseModel (TODO: replace and confirm BMEG id def)
    dr = DrugResponse(**{"id": id, "projectId": id, "submitterId": id, "acc": acc, "auc":auc,
                    "dss1": dss1, "ec50": ec50, "einf": einf, "hs": hs, "ic50": ic50,
                    "substances": [substance]})
    return orjson.loads(substance.json())


def build_patient(row):
    # patient BMEG BaseModel (TODO: substanceDef schema) 
    # Specimen has - sample.type info ex. cancer / normal

    patient_id = str(uuid.uuid4()) # TODO: replace and confirm BMEG id def
    patient_identifier = Identifier(**{"value": row['improve_sample_id'], "system": "https://pnnl-compbio.github.io/coderdata/"})

    substance_row = row[['pubchem_id', 'canSMILES', 'InChI']]
    substance = build_substance(substance_row)

    patient = Patient(**{"id": patient_id, "identifier": patient_identifier, "substances": [substance]})
    return orjson.loads(patient.json())

