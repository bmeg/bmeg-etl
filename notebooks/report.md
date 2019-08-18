# Overview

Method:

* We downloaded the 3 target datasets, creating a simple graph from each
* Apply simple transform to harmonize [Patient, Case] to 'Subject'
* Introspect the Subject nodes, identifying the 'same_as' edge
* Illustrate the diagnosis vocabulary 

## TCIA

![image](figures/tcia_summary.png)

## PDC

![image](figures/pdc_summary.png)

## GDC

![image](figures/gdc_summary.png)

## Composed CDA 

![image](figures/cda_summary.png)


## Shared cases distribution
* gdc 2053
* tcia 1791
* pdc 262

![image](figures/shared_cases_details.png)

## Diagnoses

* diagnoses in pdc & gdc 9
* diagnoses in pdc not found in gdc 18
* diagnoses in gdc not found in pdc 179
* note: no diagnoses in tcia

![image](figures/primary_diagnoses.png)
