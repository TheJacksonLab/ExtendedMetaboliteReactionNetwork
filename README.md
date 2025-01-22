# ExtendedMetaboliteReactionNetwork

This repository contains the python scripts to generate and analyze the Extended Metabolite Reaction Network (EMRN). All the input and output data files relevant to the above codes are present in [Zenodo](https://doi.org/10.5281/zenodo.14719996).

Packages installed:
* rdkit 2023.9.5
* rdkit-pypi 2022.9.4
* pandas 1.5.3
* scikit-learn 1.2.0
* scipy  1.8.1
* opentsne 1.0.1
* openbabel 3.1.1
* tqdm 4.66.4

## Data description
**1. ecmdb.json**
_E. coli_ metabolite database. Downloaded from [ECMDB](https://ecmdb.ca/).

**2. yeast_detected_and_quantified.csv**
_S. cerevisiae_ metabolite database. Downloaded from [YMDB](https://www.ymdb.ca/). Only included metabolites that have been Detected and Quantified.

**PaMet.xlsx**
_P. aeruginosa_ metabolite database. Downloaded from [PAMDB](http://pseudomonas.umaryland.edu/).

**met_all_noMW_common.csv**
Output database after combining and filtering the metabolite databases (Metabolite_db_gen/metabolite_df_generaton.py)

**USPTO_FULL.csv**
The complete USPTO reaction database.

**identified_reactants.csv**
USPTO database after carrying out reaction role mapping to separate the reactants and reagents in the reaction SMILES. Code implemented from [https://github.com/rdkit/rdkit/tree/master/Contrib/RxnRoleAssignment](https://github.com/rdkit/rdkit/tree/master/Contrib/RxnRoleAssignment).

**id_reac_final.csv**
USPTO database after filtering (EMRN_gen/USPTO_filtering_for_EMRN.py). Used for final EMRN generation.

**metr_sus_df.csv**
Database of all molecules in EMRN tagged with which round they are formed in first.

**sus_rxn5.csv**
Database of all the reactions in EMRN. Columns are Met_reactant (metabolite-based reactant in the reaction), Product (product formed in the reaction), Rxn_idx (reaction index based on the index of the given reaction in id_reac), Round (reaction round).

**OMG_monomers.csv, OMG_polymers.csv**
Open Macromolecular Genome databases. Downloaded from [https://zenodo.org/records/7556992](https://zenodo.org/records/7556992).

## 1. Metabolite_db_gen

Scripts for combining and filtering the metabolite database, which will form the core of the EMRN.

## 2. EMRN_gen

Scripts for filtering the reaction-role mapped USPTO database, and generating the EMRN database.

## 3. scscore_mw

Scripts for calculating the Synthetic Complexity Score and Molecular Weight, which have been the selected complexity metrics for the molecules in this study.

## 4. appl_codes

Scripts for some of the analysis done in this work. Functional group distribution of EMRN molecules across reaction rounds and with respect to OMG functionalities, calculation of Morgan fingerpints and t-SNE and the analysis of EMRN molecules with respect to the selected petrochemical building blocks, bio-based platform chemicals and popular commodity monomers.

