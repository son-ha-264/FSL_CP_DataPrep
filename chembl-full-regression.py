import argparse
import os
import sqlite3
from itertools import compress

import numpy as np
import progressbar
from scipy.io import mmwrite
from scipy.sparse import lil_matrix

from utils.queries import query_all_compounds, query_compounds_pchembl_value, query_assays_pchembl


def save_sparse_matrices(out_path, filename, values, counts):
    # mat_dict = {}
    # counts_coo = counts.tocoo()
    values_coo = values.tocoo()
    # coo_to_dict(counts_coo, mat_dict, "counts")
    # values_coo.data = values_coo.data / counts_coo.data
    # coo_to_dict(values_coo, mat_dict, "values")
    # --
    # np.savez(os.path.join(out_path, filename + ".npz"), **mat_dict)
    mmwrite(os.path.join(out_path, filename + "-values"), values_coo)
    # mmwrite(os.path.join(out_path, filename + "-counts"), counts_coo)


def coo_to_dict(coo, dct, name):
    dct[name + "_row"] = coo.row
    dct[name + "_col"] = coo.col
    dct[name + "_data"] = coo.data
    dct[name + "_shape"] = coo.shape


def save_index(filename, inchi, smiles):
    with open(filename, "w") as f:
        f.write("INDEX,INCHIKEY,SMILES\n")
        for i, values in enumerate(zip(inchi, smiles)):
            f.write("{},{},{}\n".format(i, values[0], values[1]))


parser = argparse.ArgumentParser()
parser.add_argument("--db", default=None, type=str, help="path to ChEMBL sqlite database ")
parser.add_argument("--out", default=".", type=str, help="path to output directory")
parser.add_argument("--version", default="24", type=str, help="chembl version to use as prefix for output files")
args = parser.parse_args()

# File Paths
chembl_db_path = "/home/son.ha/FSL_CP_DataPrep/sql/chembl_29_sqlite/chembl_29.db"
out_path = "/home/son.ha/FSL_CP_DataPrep/temp"
chembl_version_prefix = "chembl29"
"""
chembl_db_path = args.db
out_path = args.out
chembl_version_prefix = "chembl{}".format(args.version)
"""

# Open ChemBL
db = None
try:
    db = sqlite3.connect(chembl_db_path)
    cursor = db.cursor()
    # Query Compounds
    print("Preparing Compound Index")
    compound_result = cursor.execute(query_all_compounds).fetchall()
    compounds = [c[0] for c in compound_result]
    compounds_smiles = [c[1] for c in compound_result]
    compound_dict = {c: i for i, c in enumerate(compounds)}
    # Query assays
    print("Preparing Assay Index")
    assays = [a[0] for a in cursor.execute(query_assays_pchembl).fetchall()]
    assay_dict = {a: i for i, a in enumerate(assays)}
    # Query relevant pchembl values
    print("Query all values")
    pchembl_values = cursor.execute(query_compounds_pchembl_value).fetchall()

    # Store number of compounds and assays
    n_assays = len(assays)
    n_compounds = len(compounds)

    print("Found {} assays for {} compounds".format(n_assays, n_compounds))

    # Set up progress bar
    _pbw = ['Sanitizing ChEMBL:', progressbar.ETA()]
    progress = progressbar.ProgressBar(widgets=_pbw, maxval=len(pchembl_values)).start()

    # Label matrix
    regression_values = lil_matrix((n_compounds, n_assays), dtype=float)
    counts = lil_matrix((n_compounds, n_assays), dtype=int)

    # Query assays for compounds
    for i, measurement in enumerate(pchembl_values):
        inchikey = measurement[0]
        std_type = measurement[1]
        pchembl = measurement[2]
        assay = measurement[3]

        # find cell in label matrix
        row = compound_dict[inchikey]
        col = assay_dict[assay]

        # set value
        regression_values[row, col] += pchembl
        counts[row, col] += 1

        # update progress
        progress.update(i)

    # remove zero rows
    print("Pruning...")
    nnz_rows = regression_values.getnnz(1) > 0
    compounds = list(compress(compounds, nnz_rows))
    compounds_smiles = list(compress(compounds_smiles, nnz_rows))
    classification_labels = regression_values[nnz_rows]
    counts = counts[nnz_rows]

    print("Saving...")
    # regression_labels.data = regression_labels.data / counts.data
    save_sparse_matrices(out_path, "{}-regression".format(chembl_version_prefix), regression_values, counts)
    save_index(os.path.join(out_path, "{}-regression-compound-index.csv".format(chembl_version_prefix)), compounds, compounds_smiles)
    np.savetxt(os.path.join(out_path, "{}-regression-assay-index.csv".format(chembl_version_prefix)), assays, fmt="%d", header="ASSAY_ID", comments="")

finally:
    progress.finish()
    if db:
        db.close()
    print("Fin")
