#!/usr/bin/env python

import pandas as pd
import coderdata as cd


beataml = cd.DatasetLoader("beataml")
mpnst = cd.DatasetLoader("mpnst")
depmap = cd.DatasetLoader("broad_sanger")  # handled via pharmacodb pipeline (can be handled here)


bm = cd.join_datasets(depmap, beataml, mpnst)
bm.experiments.loc[:, 'study'] = bm.experiments['study'].str.replace(' ', '_', regex=False)
bm.experiments.loc[:, 'source'] = bm.experiments['source'].str.replace(' ', '_', regex=False)
bm.experiments.loc[:, 'improve_sample_id'] = pd.to_numeric(bm.experiments['improve_sample_id'], errors='coerce').astype(
    'Int64').astype('str')

bm.drugs['pubchem_id'] = pd.to_numeric(bm.drugs['pubchem_id'], errors='coerce').astype('Int64').astype('str')
bm.drugs = bm.drugs.drop_duplicates(subset=['canSMILES', 'InChIKey', 'improve_drug_id'])
bm.samples = bm.samples.drop_duplicates(subset=['cancer_type', 'model_type', 'improve_sample_id'])


bm.samples['improve_sample_id'] = pd.to_numeric(bm.samples['improve_sample_id'], errors='coerce').astype('Int64').astype('str')

bm.experiments["id"] = bm.experiments["study"] + ":" + bm.experiments["source"] + ":" + bm.experiments[
    "improve_sample_id"] + ":" + bm.experiments["improve_drug_id"]

experiments = bm.experiments[bm.experiments.dose_response_metric.isin(
    ['fit_auc', 'fit_ic50', 'fit_ec50', 'fit_einf', 'fit_hs', 'aac', 'dss'])]

experiments.loc[:, 'dose_response_metric'] = experiments.dose_response_metric.str.replace('fit_auc', 'auc', regex=False)
experiments.loc[:, 'dose_response_metric'] = experiments.dose_response_metric.str.replace('fit_ic50', 'ic50',
                                                                                          regex=False)
experiments.loc[:, 'dose_response_metric'] = experiments.dose_response_metric.str.replace('fit_ec50', 'ec50',
                                                                                          regex=False)
experiments.loc[:, 'dose_response_metric'] = experiments.dose_response_metric.str.replace('fit_einf', 'einf',
                                                                                          regex=False)
experiments.loc[:, 'dose_response_metric'] = experiments.dose_response_metric.str.replace('fit_hs', 'hs', regex=False)

experiments_pivot = experiments.pivot(index="id", columns=['dose_response_metric'], values="dose_response_value")
experiments_pivot.fillna(0, inplace=True)
experiments_pivot.replace('nan', 0, inplace=True)
experiments_pivot['id'] = experiments_pivot.index

df = experiments_pivot['id'].str.split(':', expand=True)
df.columns = ['study', 'source', 'improve_sample_id', 'improve_drug_id']
df2 = pd.concat([experiments_pivot, df], axis=1)

# create experiments / dugs table for metadata required to build DrugResponse
df3 = df2.merge(bm.drugs, on='improve_drug_id', how='left')
df3 = df3.drop_duplicates(subset=['id'])

# create samples / dugs table for metadata required to build Patient
df4 = df3.merge(bm.samples, on='improve_sample_id', how='left')
df4 = df4.drop_duplicates(subset=['id'])

# save data wrangled 
experiments_pivot.to_csv("../../source/coderdata/cleaned_pivoted_experiments_drugs_samples.csv", index=True)
# bm.experiments.to_csv("../../source/coderdata/bm_experiments.csv", index=False)
# experiments.to_csv("../../source/coderdata/cleaned_experiments.csv", index=False)
# experiments_pivot.to_csv("../../source/coderdata/cleaned_pivoted_experiments.csv", index=True)
# bm.drugs.to_csv("../../source/coderdata/bm_drugs.csv", index=False)
# bm.samples.to_csv("../../source/coderdata/bm_samples.csv", index=False)
