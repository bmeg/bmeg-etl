#!/usr/bin/env python

import json
import datetime
import pandas as pd 

def assign_sex(x):
    if x == 2: 
        return "female"
    elif x == 1:
        return "male"
    else: 
        pass

def get_snomed_code(x, gtex_tissue_type):
    for item in gtex_tissue_type:
        if item["value"] == x:
            return item["sctid"]

sample_pheno = pd.read_csv("../../source/gtex/GTEx_Analysis_v8_Annotations_SampleAttributesDS.txt", sep="\t")
patient_pheno = pd.read_csv("../../source/gtex/GTEx_Analysis_v8_Annotations_SubjectPhenotypesDS.txt", sep="\t")


patient_pheno["SEX_string"] = patient_pheno.SEX.apply(lambda x: assign_sex(x))
patient_pheno['fhir_birthDate'] = patient_pheno.AGE.apply(lambda x: datetime.datetime(datetime.datetime.now().year - int(x.split('-')[1]), 1, 1)) # needs to be reverted back to int for analysis 


sample_pheno['SUBJID'] = sample_pheno.SAMPID.apply(lambda x: "-".join(x.split('-', 2)[0:2]))
sample_pheno = sample_pheno.merge(patient_pheno, on='SUBJID', how='left')

tissue_list = list(sample_pheno.SMTSD.unique())

d = []
for item in tissue_list:
    d.append({"value": item, "sctid": ""})

# snomed was annotated via SNOMED CT Browser
# with open("../../source/gtex/gtex_tissue_type_empty.json", "w") as f:
#    json.dump(d, f)

with open("../../source/gtex/gtex_tissue_type.json", encoding='utf-8') as f:
   gtex_tissue_type  = json.load(f)


sample_pheno['sctid'] = sample_pheno.SMTSD.apply(lambda x: get_snomed_code(x, gtex_tissue_type))

# save data wrangled 
sample_pheno.to_csv("../../source/gtex/sample_pheno.csv", index=False)
patient_pheno.to_csv("../../source/gtex/patient_pheno.csv", index=False)
