#!/usr/bin/env python

import pandas as pd
import coderdata as cd
from fhir.resources.substance import Substance
from fhir.resources.reference import Reference
from fhir.resources.identifier import Identifier
from fhir.resources.codeablereference import CodeableReference
from fhir.resources.codeableconcept import CodeableConcept
from fhir.resources.substancedefinition import SubstanceDefinition, SubstanceDefinitionStructure, \
    SubstanceDefinitionStructureRepresentation
from bmegDrugResponse import DrugResponse
from bmegPatient import Patient


beataml = cd.DatasetLoader("beataml")
mpnst = cd.DatasetLoader("mpnst")
depmap = cd.DatasetLoader("broad_sanger")  # handled via pharmacodb pipeline (can be handled here)


bm = cd.join_datasets(depmap, beataml, mpnst)
bm.experiments.loc[:, 'study'] = bm.experiments['study'].str.replace(' ', '_', regex=False)
bm.experiments.loc[:, 'source'] = bm.experiments['source'].str.replace(' ', '_', regex=False)
bm.experiments.loc[:, 'improve_sample_id'] = pd.to_numeric(bm.experiments['improve_sample_id'], errors='coerce').astype(
    'Int64').astype('str')

bm.drugs['pubchem_id'] = pd.to_numeric(bm.drugs['pubchem_id'], errors='coerce').astype('Int64').astype('str')

bm.experiments["id"] = bm.experiments["study"] + ":" + bm.experiments["source"] + ":" + bm.experiments[
    "improve_sample_id"] + ":" + bm.experiments["improve_drug_id"]

experiments = bm.experiments[bm.experiments.dose_response_metric.isin(
    ['fit_auc', 'fit_ic50', 'fit_ec50', 'fit_r2', 'fit_ec50se', 'fit_einf', 'fit_hs', 'aac', 'dss'])]

experiments.loc[:, 'dose_response_metric'] = experiments.dose_response_metric.str.replace('fit_auc', 'auc', regex=False)
experiments.loc[:, 'dose_response_metric'] = experiments.dose_response_metric.str.replace('fit_ic50', 'ic50',
                                                                                          regex=False)
experiments.loc[:, 'dose_response_metric'] = experiments.dose_response_metric.str.replace('fit_ec50', 'ec50',
                                                                                          regex=False)
experiments.loc[:, 'dose_response_metric'] = experiments.dose_response_metric.str.replace('fit_einf', 'einf',
                                                                                          regex=False)
experiments.loc[:, 'dose_response_metric'] = experiments.dose_response_metric.str.replace('fit_hs', 'hs', regex=False)

# save data wrangled 
bm.experiments.to_csv("../../source/coderdata/bm_experiments.csv", index=False)
bm.drugs.to_csv("../../source/coderdata/bm_drugs.csv", index=False)
bm.samples.to_csv("../../source/coderdata/bm_samples.csv", index=False)


# ----------------------------------------------------------------
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
# patient BMEG BaseModel 
# patient_id = "PLACE-HOLDER-3456"
# patient_identifier = Identifier(**{"value": "PLACE-HOLDER-1234", "system": "https://pnnl-compbio.github.io/coderdata/"})
# patient = Patient(**{"id": patient_id, "substances": [substance]})

# drugResponse BaseModel
# DrugResponse(**{"id": "BeatAML:synapse:4819:SMI_1006", "projectId": "123456", "submitterId": "1234", "substances": [substance]})

