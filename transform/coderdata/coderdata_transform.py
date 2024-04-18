#!/usr/bin/env python

import pandas as pd
from fhir.resources.substance import Substance
from fhir.resources.reference import Reference
from fhir.resources.identifier import Identifier
from fhir.resources.codeablereference import CodeableReference
from fhir.resources.codeableconcept import CodeableConcept
from fhir.resources.substancedefinition import SubstanceDefinition, SubstanceDefinitionStructure, \
    SubstanceDefinitionStructureRepresentation
from bmegDrugResponse import DrugResponse
from bmegPatient import Patient

# build plan
# substance
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
    substanceDef_ident = Identifier(**{"value": "PLACE-HOLDER-1234", "system": "https://pubchem.ncbi.nlm.nih.gov/"})
    sd = SubstanceDefinition(**{"id": id, "identifier": [substanceDef_ident], "structure": sds})

    cr = CodeableReference(**{"reference": Reference(**{"reference": "/".join(["SubstanceDefinition", id])})})
    substance = Substance(**{"instance": True, "code": cr})
    return substance


# -----------------------------------------------------------------
# patient BMEG BaseModel (TODO: substanceDef schema) 
# patient_id = "PLACE-HOLDER-3456"
# patient_identifier = Identifier(**{"value": "PLACE-HOLDER-1234", "system": "https://pnnl-compbio.github.io/coderdata/"})
# patient = Patient(**{"id": patient_id, "substances": [substance]})

# drugResponse BMEG BaseModel (TODO: confirm BMEG id def)
# DrugResponse(**{"id": "BeatAML:synapse:4819:SMI_1006", "projectId": "123456", "submitterId": "1234", "substances": [substance]})

