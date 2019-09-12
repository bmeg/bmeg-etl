
import bmeg.enrichers.drug_enricher as drug_enricher
import logging
from pprint import pprint

# maintain fresh cache for all calls
from uncacher import uncache
uncache(drug_enricher.requests)


NAMES = """
asprin
5-Fluorouracil
(5Z)-7-Oxozeaenol
A-443654
A-770041
Afatinib (1)
Afatinib (2)
AICA Ribonucleotide
AKT inhibitor VIII (1)
Ecotrin
Tomaxifen
Tamoxiten
Abagovomag
Zometa
""".strip().split("\n")

EXPECTED = [
    {'id': 'CID2244', 'id_source': 'PUBCHEM', 'pubchem_id': 'CID2244', 'chebi_id': 'CHEBI:15365', 'chembl_id': 'CHEMBL25', 'drugbank_id': 'DB00945', 'synonym': 'ASPIRIN', 'inchi': 'InChI=1S/C9H8O4/c1-6(10)13-8-5-3-2-4-7(8)9(11)12/h2-5H,1H3,(H,11,12)', 'inchi_key': 'BSYNRYMUTXBXSQ-UHFFFAOYSA-N', 'taxonomy': {'class': 'Benzene and substituted derivatives', 'description': 'This compound belongs to the class of organic compounds known as acylsalicylic acids. These are o-acylated derivatives of salicylic acid.', 'direct-parent': 'Acylsalicylic acids', 'kingdom': 'Organic compounds', 'subclass': 'Benzoic acids and derivatives', 'superclass': 'Benzenoids'}, 'approved_countries': ['Canada', 'EU', 'US'], 'usan_stem_definition': None, 'source_url': 'http://mychem.info/v1/query?q=pubchem.cid:"2244"&fields=chebi.id,chebi.inchi,chebi.inchi_key,chebi.name,chembl.molecule_chembl_id,chembl.pref_name,chembl.inchi,chembl.inchi_key,chembl.molecule_synonyms,chembl.usan_stem_definition,pubchem.cid,pubchem.inchi,pubchem.inchi_key,drugbank.id,drugbank.inchi,drugbank.inchi_key,drugbank.products.approved,drugbank.products.country,drugbank.taxonomy.class,drugbank.taxonomy.direct-parent,drugbank.taxonomy.kingdom,drugbank.taxonomy.subclass,drugbank.taxonomy.superclass,drugbank.taxonomy.description&size=1'},
    {'id': 'CID3385', 'id_source': 'PUBCHEM', 'pubchem_id': 'CID3385', 'chebi_id': 'CHEBI:46345', 'chembl_id': 'CHEMBL185', 'drugbank_id': 'DB00544', 'synonym': 'FLUOROURACIL', 'inchi': 'InChI=1S/C4H3FN2O2/c5-2-1-6-4(9)7-3(2)8/h1H,(H2,6,7,8,9)', 'inchi_key': 'GHASVSINZRGABV-UHFFFAOYSA-N', 'taxonomy': {'class': 'Diazines', 'description': 'This compound belongs to the class of organic compounds known as halopyrimidines. These are aromatic compounds containing a halogen atom linked to a pyrimidine ring. Pyrimidine is a 6-membered ring consisting of four carbon atoms and two nitrogen centers at the 1- and 3- ring positions.', 'direct-parent': 'Halopyrimidines', 'kingdom': 'Organic compounds', 'subclass': 'Pyrimidines and pyrimidine derivatives', 'superclass': 'Organoheterocyclic compounds'}, 'approved_countries': ['Canada', 'US'], 'usan_stem_definition': 'uracil derivatives used as thyroid antagonists and as antineoplastics', 'source_url': 'http://mychem.info/v1/query?q=chembl.pref_name:Fluorouracil&fields=chebi.id,chebi.inchi,chebi.inchi_key,chebi.name,chembl.molecule_chembl_id,chembl.pref_name,chembl.inchi,chembl.inchi_key,chembl.molecule_synonyms,chembl.usan_stem_definition,pubchem.cid,pubchem.inchi,pubchem.inchi_key,drugbank.id,drugbank.inchi,drugbank.inchi_key,drugbank.products.approved,drugbank.products.country,drugbank.taxonomy.class,drugbank.taxonomy.direct-parent,drugbank.taxonomy.kingdom,drugbank.taxonomy.subclass,drugbank.taxonomy.superclass,drugbank.taxonomy.description&size=1'},
    {'id': 'CID9863776', 'id_source': 'PUBCHEM', 'pubchem_id': 'CID9863776', 'chebi_id': 'CHEBI:67559', 'chembl_id': 'CHEMBL1077979', 'drugbank_id': None, 'synonym': '(5Z)-7-Oxozeaenol', 'inchi': 'InChI=1S/C19H22O7/c1-11-5-3-7-14(20)18(23)15(21)8-4-6-12-9-13(25-2)10-16(22)17(12)19(24)26-11/h3-4,6-7,9-11,15,18,21-23H,5,8H2,1-2H3/b6-4+,7-3-/t11-,15-,18+/m0/s1', 'inchi_key': 'NEQZWEXWOFPKOT-BYRRXHGESA-N', 'taxonomy': None, 'approved_countries': [], 'usan_stem_definition': None, 'source_url': 'http://mychem.info/v1/query?q=chebi.name:"(5Z)-7-Oxozeaenol"&fields=chebi.id,chebi.inchi,chebi.inchi_key,chebi.name,chembl.molecule_chembl_id,chembl.pref_name,chembl.inchi,chembl.inchi_key,chembl.molecule_synonyms,chembl.usan_stem_definition,pubchem.cid,pubchem.inchi,pubchem.inchi_key,drugbank.id,drugbank.inchi,drugbank.inchi_key,drugbank.products.approved,drugbank.products.country,drugbank.taxonomy.class,drugbank.taxonomy.direct-parent,drugbank.taxonomy.kingdom,drugbank.taxonomy.subclass,drugbank.taxonomy.superclass,drugbank.taxonomy.description&size=1'},
    {'id': 'CID10172943', 'id_source': 'PUBCHEM', 'pubchem_id': 'CID10172943', 'chebi_id': 'CHEBI:91351', 'chembl_id': 'CHEMBL379300', 'drugbank_id': 'DB08073', 'synonym': 'A-443654', 'inchi': 'InChI=1S/C24H23N5O/c1-15-22-10-16(6-7-24(22)29-28-15)17-9-20(13-26-11-17)30-14-19(25)8-18-12-27-23-5-3-2-4-21(18)23/h2-7,9-13,19,27H,8,14,25H2,1H3,(H,28,29)/t19-/m0/s1', 'inchi_key': 'YWTBGJGMTBHQTM-IBGZPJMESA-N', 'taxonomy': {'class': 'Indoles and derivatives', 'description': 'This compound belongs to the class of organic compounds known as 3-alkylindoles. These are compounds containing an indole moiety that carries an alkyl chain at the 3-position.', 'direct-parent': '3-alkylindoles', 'kingdom': 'Organic compounds', 'subclass': 'Indoles', 'superclass': 'Organoheterocyclic compounds'}, 'approved_countries': [], 'usan_stem_definition': None, 'source_url': 'http://mychem.info/v1/query?q=chembl.pref_name:"A-443654"&fields=chebi.id,chebi.inchi,chebi.inchi_key,chebi.name,chembl.molecule_chembl_id,chembl.pref_name,chembl.inchi,chembl.inchi_key,chembl.molecule_synonyms,chembl.usan_stem_definition,pubchem.cid,pubchem.inchi,pubchem.inchi_key,drugbank.id,drugbank.inchi,drugbank.inchi_key,drugbank.products.approved,drugbank.products.country,drugbank.taxonomy.class,drugbank.taxonomy.direct-parent,drugbank.taxonomy.kingdom,drugbank.taxonomy.subclass,drugbank.taxonomy.superclass,drugbank.taxonomy.description&size=1'},
    {'id': 'CID9549184', 'id_source': 'PUBCHEM', 'pubchem_id': 'CID9549184', 'chebi_id': 'CHEBI:91457', 'chembl_id': 'CHEMBL1970879', 'drugbank_id': None, 'synonym': '9549184', 'inchi': 'InChI=1S/C34H39N9O3/c1-21(44)41-14-16-42(17-15-41)24-9-11-25(12-10-24)43-33-30(32(35)36-20-37-33)31(39-43)23-8-13-26(29(19-23)46-3)38-34(45)28-18-22-6-4-5-7-27(22)40(28)2/h4-8,13,18-20,24-25H,9-12,14-17H2,1-3H3,(H,38,45)(H2,35,36,37)', 'inchi_key': 'ZMNWFTYYYCSSTF-UHFFFAOYSA-N', 'taxonomy': None, 'approved_countries': [], 'usan_stem_definition': None, 'source_url': 'http://mychem.info/v1/query?q=pubchem.cid:"9549184"&fields=chebi.id,chebi.inchi,chebi.inchi_key,chebi.name,chembl.molecule_chembl_id,chembl.pref_name,chembl.inchi,chembl.inchi_key,chembl.molecule_synonyms,chembl.usan_stem_definition,pubchem.cid,pubchem.inchi,pubchem.inchi_key,drugbank.id,drugbank.inchi,drugbank.inchi_key,drugbank.products.approved,drugbank.products.country,drugbank.taxonomy.class,drugbank.taxonomy.direct-parent,drugbank.taxonomy.kingdom,drugbank.taxonomy.subclass,drugbank.taxonomy.superclass,drugbank.taxonomy.description&size=1'},
    {'id': 'CID10184653', 'id_source': 'PUBCHEM', 'pubchem_id': 'CID10184653', 'chebi_id': 'CHEBI:61390', 'chembl_id': 'CHEMBL1173655', 'drugbank_id': 'DB08916', 'synonym': 'AFATINIB', 'inchi': 'InChI=1S/C24H25ClFN5O3/c1-31(2)8-3-4-23(32)30-21-11-17-20(12-22(21)34-16-7-9-33-13-16)27-14-28-24(17)29-15-5-6-19(26)18(25)10-15/h3-6,10-12,14,16H,7-9,13H2,1-2H3,(H,30,32)(H,27,28,29)/b4-3+/t16-/m0/s1', 'inchi_key': 'ULXXDDBFHOBEHA-CWDCEQMOSA-N', 'taxonomy': {'class': 'Diazanaphthalenes', 'description': 'This compound belongs to the class of organic compounds known as quinazolinamines. These are heterocyclic aromatic compounds containing a quianazoline moiety substituted by one or more amine groups.', 'direct-parent': 'Quinazolinamines', 'kingdom': 'Organic compounds', 'subclass': 'Benzodiazines', 'superclass': 'Organoheterocyclic compounds'}, 'approved_countries': ['Canada', 'US'], 'usan_stem_definition': 'tyrosine kinase inhibitors', 'source_url': 'http://mychem.info/v1/query?q=chembl.pref_name:Afatinib&fields=chebi.id,chebi.inchi,chebi.inchi_key,chebi.name,chembl.molecule_chembl_id,chembl.pref_name,chembl.inchi,chembl.inchi_key,chembl.molecule_synonyms,chembl.usan_stem_definition,pubchem.cid,pubchem.inchi,pubchem.inchi_key,drugbank.id,drugbank.inchi,drugbank.inchi_key,drugbank.products.approved,drugbank.products.country,drugbank.taxonomy.class,drugbank.taxonomy.direct-parent,drugbank.taxonomy.kingdom,drugbank.taxonomy.subclass,drugbank.taxonomy.superclass,drugbank.taxonomy.description&size=1'},
    {'id': 'CID10184653', 'id_source': 'PUBCHEM', 'pubchem_id': 'CID10184653', 'chebi_id': 'CHEBI:61390', 'chembl_id': 'CHEMBL1173655', 'drugbank_id': 'DB08916', 'synonym': 'AFATINIB', 'inchi': 'InChI=1S/C24H25ClFN5O3/c1-31(2)8-3-4-23(32)30-21-11-17-20(12-22(21)34-16-7-9-33-13-16)27-14-28-24(17)29-15-5-6-19(26)18(25)10-15/h3-6,10-12,14,16H,7-9,13H2,1-2H3,(H,30,32)(H,27,28,29)/b4-3+/t16-/m0/s1', 'inchi_key': 'ULXXDDBFHOBEHA-CWDCEQMOSA-N', 'taxonomy': {'class': 'Diazanaphthalenes', 'description': 'This compound belongs to the class of organic compounds known as quinazolinamines. These are heterocyclic aromatic compounds containing a quianazoline moiety substituted by one or more amine groups.', 'direct-parent': 'Quinazolinamines', 'kingdom': 'Organic compounds', 'subclass': 'Benzodiazines', 'superclass': 'Organoheterocyclic compounds'}, 'approved_countries': ['Canada', 'US'], 'usan_stem_definition': 'tyrosine kinase inhibitors', 'source_url': 'http://mychem.info/v1/query?q=chembl.pref_name:Afatinib&fields=chebi.id,chebi.inchi,chebi.inchi_key,chebi.name,chembl.molecule_chembl_id,chembl.pref_name,chembl.inchi,chembl.inchi_key,chembl.molecule_synonyms,chembl.usan_stem_definition,pubchem.cid,pubchem.inchi,pubchem.inchi_key,drugbank.id,drugbank.inchi,drugbank.inchi_key,drugbank.products.approved,drugbank.products.country,drugbank.taxonomy.class,drugbank.taxonomy.direct-parent,drugbank.taxonomy.kingdom,drugbank.taxonomy.subclass,drugbank.taxonomy.superclass,drugbank.taxonomy.description&size=1'},
    {'id': 'CID65110', 'id_source': 'PUBCHEM', 'pubchem_id': 'CID65110', 'chebi_id': 'CHEBI:18406', 'chembl_id': 'CHEMBL483849', 'drugbank_id': 'DB01700', 'synonym': 'AICA Ribonucleotide', 'inchi': 'InChI=1S/C9H15N4O8P/c10-7-4(8(11)16)12-2-13(7)9-6(15)5(14)3(21-9)1-20-22(17,18)19/h2-3,5-6,9,14-15H,1,10H2,(H2,11,16)(H2,17,18,19)/t3-,5-,6-,9-/m1/s1', 'inchi_key': 'NOTGFIUVDGNKRI-UUOKFMHZSA-N', 'taxonomy': {'class': 'Imidazole ribonucleosides and ribonucleotides', 'description': 'This compound belongs to the class of organic compounds known as 1-ribosyl-imidazolecarboxamides. These are organic compounds containing the imidazole ring linked to a ribose ring through a 1-2 bond.', 'direct-parent': '1-ribosyl-imidazolecarboxamides', 'kingdom': 'Organic compounds', 'subclass': '1-ribosyl-imidazolecarboxamides', 'superclass': 'Nucleosides, nucleotides, and analogues'}, 'approved_countries': [], 'usan_stem_definition': None, 'source_url': 'http://mychem.info/v1/query?q=chebi.name:"AICA Ribonucleotide"&fields=chebi.id,chebi.inchi,chebi.inchi_key,chebi.name,chembl.molecule_chembl_id,chembl.pref_name,chembl.inchi,chembl.inchi_key,chembl.molecule_synonyms,chembl.usan_stem_definition,pubchem.cid,pubchem.inchi,pubchem.inchi_key,drugbank.id,drugbank.inchi,drugbank.inchi_key,drugbank.products.approved,drugbank.products.country,drugbank.taxonomy.class,drugbank.taxonomy.direct-parent,drugbank.taxonomy.kingdom,drugbank.taxonomy.subclass,drugbank.taxonomy.superclass,drugbank.taxonomy.description&size=1'},
    {'id': 'CID135398501', 'id_source': 'PUBCHEM', 'pubchem_id': 'CID135398501', 'chebi_id': None, 'chembl_id': None, 'drugbank_id': None, 'synonym': '135398501', 'inchi': 'InChI=1S/C34H29N7O/c42-34-39-26-8-4-5-9-31(26)41(34)25-14-16-40(17-15-25)20-22-10-12-24(13-11-22)33-32(23-6-2-1-3-7-23)37-29-18-27-28(36-21-35-27)19-30(29)38-33/h1-13,18-19,21,25,38H,14-17,20H2,(H,39,42)', 'inchi_key': 'IWCQHVUQEFDRIW-UHFFFAOYSA-N', 'taxonomy': None, 'approved_countries': [], 'usan_stem_definition': None, 'source_url': 'http://mychem.info/v1/query?q=pubchem.cid:"135398501"&fields=chebi.id,chebi.inchi,chebi.inchi_key,chebi.name,chembl.molecule_chembl_id,chembl.pref_name,chembl.inchi,chembl.inchi_key,chembl.molecule_synonyms,chembl.usan_stem_definition,pubchem.cid,pubchem.inchi,pubchem.inchi_key,drugbank.id,drugbank.inchi,drugbank.inchi_key,drugbank.products.approved,drugbank.products.country,drugbank.taxonomy.class,drugbank.taxonomy.direct-parent,drugbank.taxonomy.kingdom,drugbank.taxonomy.subclass,drugbank.taxonomy.superclass,drugbank.taxonomy.description&size=1'},
    {'id': 'CID2244', 'id_source': 'PUBCHEM', 'pubchem_id': 'CID2244', 'chebi_id': 'CHEBI:15365', 'chembl_id': 'CHEMBL25', 'drugbank_id': 'DB00945', 'synonym': 'ASPIRIN', 'inchi': 'InChI=1S/C9H8O4/c1-6(10)13-8-5-3-2-4-7(8)9(11)12/h2-5H,1H3,(H,11,12)', 'inchi_key': 'BSYNRYMUTXBXSQ-UHFFFAOYSA-N', 'taxonomy': {'class': 'Benzene and substituted derivatives', 'description': 'This compound belongs to the class of organic compounds known as acylsalicylic acids. These are o-acylated derivatives of salicylic acid.', 'direct-parent': 'Acylsalicylic acids', 'kingdom': 'Organic compounds', 'subclass': 'Benzoic acids and derivatives', 'superclass': 'Benzenoids'}, 'approved_countries': ['Canada', 'EU', 'US'], 'usan_stem_definition': None, 'source_url': 'http://mychem.info/v1/query?q=pubchem.cid:"2244"&fields=chebi.id,chebi.inchi,chebi.inchi_key,chebi.name,chembl.molecule_chembl_id,chembl.pref_name,chembl.inchi,chembl.inchi_key,chembl.molecule_synonyms,chembl.usan_stem_definition,pubchem.cid,pubchem.inchi,pubchem.inchi_key,drugbank.id,drugbank.inchi,drugbank.inchi_key,drugbank.products.approved,drugbank.products.country,drugbank.taxonomy.class,drugbank.taxonomy.direct-parent,drugbank.taxonomy.kingdom,drugbank.taxonomy.subclass,drugbank.taxonomy.superclass,drugbank.taxonomy.description&size=1'},
    {'id': 'CID2733526', 'id_source': 'PUBCHEM', 'pubchem_id': 'CID2733526', 'chebi_id': 'CHEBI:41774', 'chembl_id': 'CHEMBL83', 'drugbank_id': 'DB00675', 'synonym': 'Tamoxifen', 'inchi': 'InChI=1S/C26H29NO/c1-4-25(21-11-7-5-8-12-21)26(22-13-9-6-10-14-22)23-15-17-24(18-16-23)28-20-19-27(2)3/h5-18H,4,19-20H2,1-3H3/b26-25-', 'inchi_key': 'NKANXQFJJICGDU-QPLCGJKRSA-N', 'taxonomy': {'class': 'Stilbenes', 'description': 'This compound belongs to the class of organic compounds known as stilbenes. These are organic compounds containing a 1,2-diphenylethylene moiety. Stilbenes (C6-C2-C6 ) are derived from the common phenylpropene (C6-C3) skeleton building block. The introduction of one or more hydroxyl groups  to a phenyl ring lead to stilbenoids.', 'direct-parent': 'Stilbenes', 'kingdom': 'Organic compounds', 'superclass': 'Phenylpropanoids and polyketides'}, 'approved_countries': ['Canada', 'US'], 'usan_stem_definition': None, 'source_url': 'http://mychem.info/v1/query?q=chembl.pref_name:tamoxifen&fields=chebi.id,chebi.inchi,chebi.inchi_key,chebi.name,chembl.molecule_chembl_id,chembl.pref_name,chembl.inchi,chembl.inchi_key,chembl.molecule_synonyms,chembl.usan_stem_definition,pubchem.cid,pubchem.inchi,pubchem.inchi_key,drugbank.id,drugbank.inchi,drugbank.inchi_key,drugbank.products.approved,drugbank.products.country,drugbank.taxonomy.class,drugbank.taxonomy.direct-parent,drugbank.taxonomy.kingdom,drugbank.taxonomy.subclass,drugbank.taxonomy.superclass,drugbank.taxonomy.description&size=1'},
    {'id': 'CID2733526', 'id_source': 'PUBCHEM', 'pubchem_id': 'CID2733526', 'chebi_id': 'CHEBI:41774', 'chembl_id': 'CHEMBL83', 'drugbank_id': 'DB00675', 'synonym': 'Tamoxifen', 'inchi': 'InChI=1S/C26H29NO/c1-4-25(21-11-7-5-8-12-21)26(22-13-9-6-10-14-22)23-15-17-24(18-16-23)28-20-19-27(2)3/h5-18H,4,19-20H2,1-3H3/b26-25-', 'inchi_key': 'NKANXQFJJICGDU-QPLCGJKRSA-N', 'taxonomy': {'class': 'Stilbenes', 'description': 'This compound belongs to the class of organic compounds known as stilbenes. These are organic compounds containing a 1,2-diphenylethylene moiety. Stilbenes (C6-C2-C6 ) are derived from the common phenylpropene (C6-C3) skeleton building block. The introduction of one or more hydroxyl groups  to a phenyl ring lead to stilbenoids.', 'direct-parent': 'Stilbenes', 'kingdom': 'Organic compounds', 'superclass': 'Phenylpropanoids and polyketides'}, 'approved_countries': ['Canada', 'US'], 'usan_stem_definition': None, 'source_url': 'http://mychem.info/v1/query?q=chembl.pref_name:tamoxifen&fields=chebi.id,chebi.inchi,chebi.inchi_key,chebi.name,chembl.molecule_chembl_id,chembl.pref_name,chembl.inchi,chembl.inchi_key,chembl.molecule_synonyms,chembl.usan_stem_definition,pubchem.cid,pubchem.inchi,pubchem.inchi_key,drugbank.id,drugbank.inchi,drugbank.inchi_key,drugbank.products.approved,drugbank.products.country,drugbank.taxonomy.class,drugbank.taxonomy.direct-parent,drugbank.taxonomy.kingdom,drugbank.taxonomy.subclass,drugbank.taxonomy.superclass,drugbank.taxonomy.description&size=1'},
    {'id': 'CHEMBL1742981', 'id_source': 'CHEMBL', 'pubchem_id': None, 'chebi_id': None, 'chembl_id': 'CHEMBL1742981', 'drugbank_id': None, 'synonym': 'Abagovomab', 'inchi': None, 'inchi_key': None, 'taxonomy': None, 'approved_countries': [], 'usan_stem_definition': 'monoclonal antibodies', 'source_url': 'http://mychem.info/v1/query?q=chembl.pref_name:abagovomab&fields=chebi.id,chebi.inchi,chebi.inchi_key,chebi.name,chembl.molecule_chembl_id,chembl.pref_name,chembl.inchi,chembl.inchi_key,chembl.molecule_synonyms,chembl.usan_stem_definition,pubchem.cid,pubchem.inchi,pubchem.inchi_key,drugbank.id,drugbank.inchi,drugbank.inchi_key,drugbank.products.approved,drugbank.products.country,drugbank.taxonomy.class,drugbank.taxonomy.direct-parent,drugbank.taxonomy.kingdom,drugbank.taxonomy.subclass,drugbank.taxonomy.superclass,drugbank.taxonomy.description&size=1'},
    {'id': 'CID68740', 'id_source': 'PUBCHEM', 'pubchem_id': 'CID68740', 'chebi_id': 'CHEBI:46557', 'chembl_id': 'CHEMBL924', 'drugbank_id': 'DB00399', 'synonym': 'ZOLEDRONIC ACID', 'inchi': 'InChI=1S/C5H10N2O7P2/c8-5(15(9,10)11,16(12,13)14)3-7-2-1-6-4-7/h1-2,4,8H,3H2,(H2,9,10,11)(H2,12,13,14)', 'inchi_key': 'XRASPMIURGNCCH-UHFFFAOYSA-N', 'taxonomy': {'class': 'Organic phosphonic acids and derivatives', 'description': 'This compound belongs to the class of organic compounds known as bisphosphonates. These are organic compounds containing two phosphonate groups linked together through a carbon atoms.', 'direct-parent': 'Bisphosphonates', 'kingdom': 'Organic compounds', 'subclass': 'Bisphosphonates', 'superclass': 'Organic acids and derivatives'}, 'approved_countries': ['Canada', 'EU', 'US'], 'usan_stem_definition': 'calcium metabolism regulators', 'source_url': 'http://mychem.info/v1/query?q=chembl.pref_name:"Zoledronic acid"&fields=chebi.id,chebi.inchi,chebi.inchi_key,chebi.name,chembl.molecule_chembl_id,chembl.pref_name,chembl.inchi,chembl.inchi_key,chembl.molecule_synonyms,chembl.usan_stem_definition,pubchem.cid,pubchem.inchi,pubchem.inchi_key,drugbank.id,drugbank.inchi,drugbank.inchi_key,drugbank.products.approved,drugbank.products.country,drugbank.taxonomy.class,drugbank.taxonomy.direct-parent,drugbank.taxonomy.kingdom,drugbank.taxonomy.subclass,drugbank.taxonomy.superclass,drugbank.taxonomy.description&size=1'}
]


def test_simple(caplog):
    """ straightforward """
    caplog.set_level(logging.DEBUG)

    for name, expected in zip(NAMES, EXPECTED):
        actual = drug_enricher.normalize(name)
        try:
            assert actual, 'Should return value for {}'.format(name)
            assert len(actual) == len(expected), 'Should return same number of fields for {}'.format(name)
            assert actual == expected
        except Exception as e:
            print(str(e))
            pprint(actual)
            pprint(expected)
            assert False, 'Failed'


def test_alias():
    assert 'Tomaxifen' in drug_enricher.ALIASES, 'we should have an alias for Tomaxifen'
    assert drug_enricher.ALIASES['Tomaxifen'].lower() == 'tamoxifen', 'the alias for Tomaxifen should be tamoxifen'
    assert 'Tamoxiten' in drug_enricher.ALIASES, 'we should have an alias for Tamoxiten'
    assert drug_enricher.ALIASES['Tamoxiten'].lower() == 'tamoxifen', 'the alias for Tamoxiten should be tamoxifen'


def test_spell_check():
    """ """
    assert drug_enricher.spell_check('Tamoxiten')[0] == 'tamoxifen'
