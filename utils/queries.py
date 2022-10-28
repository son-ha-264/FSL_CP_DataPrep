assays_for_compound_smiles = "select s.standard_inchi_key, s.canonical_smiles, " \
                             "ay.assay_id, ay.assay_organism, ay.description, act.standard_type, act.standard_units, act.standard_relation, " \
                             "act.published_value, act.published_units, act.pchembl_value, act.activity_comment " \
                             "from compound_structures s, activities act, assays ay where canonical_smiles=?" \
                             " and act.molregno = s.molregno and act.assay_id = ay.assay_id"

assays_for_compound_inchikey = "select s.standard_inchi_key, s.canonical_smiles, " \
                               "ay.assay_id, ay.assay_organism, ay.description, act.standard_type, act.standard_units, act.standard_relation, " \
                               "act.published_value, act.published_units, act.pchembl_value, act.activity_comment " \
                               "from compound_structures s, activities act, assays ay where standard_inchi_key=?" \
                               " and act.molregno = s.molregno and act.assay_id = ay.assay_id"

query_all_compounds = "select standard_inchi_key, canonical_smiles from compound_structures"

query_compounds_pchembl_value = "select standard_inchi_key, act.standard_type, act.pchembl_value, act.assay_id " \
                                "from compound_structures s, activities act " \
                                "where act.standard_type in (\"pIC50\", \"pEC50\", \"IC50\", \"EC50\", \"-Log EC50\", \"-Log IC50\", \"Potency\", \"GI50\", \"AC50\", \"ED50\") and " \
                                "act.pchembl_value is not NULL and act.pchembl_value BETWEEN 4 and 10 and act.standard_relation=\"=\" and act.molregno = s.molregno"

query_assays_pchembl = "select distinct(act.assay_id) " \
                                     "from compound_structures s, activities act " \
                                     "where act.standard_type in (\"pIC50\", \"pEC50\", \"IC50\", \"EC50\", \"-Log EC50\", \"-Log IC50\", \"Potency\", \"GI50\", \"AC50\", \"ED50\") and " \
                                     "act.pchembl_value is not NULL and act.pchembl_value BETWEEN 4 and 10 and act.standard_relation=\"=\" and act.molregno = s.molregno"

query_compounds_active_inactive = "select standard_inchi_key, act.standard_type, act.activity_comment, act.assay_id " \
                                  "from compound_structures s, activities act " \
                                  "where act.standard_type in (\"pIC50\", \"pEC50\", \"IC50\", \"EC50\", \"-Log EC50\", \"-Log IC50\", \"Potency\", \"GI50\", \"AC50\", \"ED50\") and " \
                                  "act.activity_comment in (\"active\", \"Active\", \"inactive\", \"Inactive\", \"Not active\", \"No activity\", \"Not Active (inhibition < 50% @ 10 uM and thus dose-reponse curve not measured)\") and act.molregno = s.molregno"

query_assays_active_inactive = "select distinct(act.assay_id) " \
                               "from compound_structures s, activities act " \
                               "where act.standard_type in (\"pIC50\", \"pEC50\", \"IC50\", \"EC50\", \"-Log EC50\", \"-Log IC50\", \"Potency\", \"GI50\", \"AC50\", \"ED50\") and " \
                               "act.activity_comment in (\"active\", \"Active\", \"inactive\", \"Inactive\", \"Not active\", \"No activity\", \"Not Active (inhibition < 50% @ 10 uM and thus dose-reponse curve not measured)\") and act.molregno = s.molregno"

