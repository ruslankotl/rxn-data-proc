# rxn-data-proc
## Scripts for the 1st year report

Reaxys data was downloaded as an RDFile, which was filtered to only contain formal borylation reaction using `borylation_cleanup.py`.

Where the script failed due to presence of Markush structures, those were removed using `markush_remover.py`.

SMILES strings can be processed in the same manner using functions from `borylation.py`.

To convert RDFiles into CSV files `rdf_to_csv.py` was used.

Automated resolution of reagents, catalysts, and solvents was performed using `csv_and_discard.py`. Any entry with at least one unresolved chemical entity was discarded.

## T5Chem model training

Setup, training, and prediction instructions are provided at https://yzhang.hpc.nyu.edu/T5Chem/tutorial.html or called using `t5chem {train, predict} -h` upon successful installation.
