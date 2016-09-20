import os
import sys

import json
from pprint import pprint

infile = sys.argv[1]

with open(infile) as data_file:    
    data = json.load(data_file)

#pprint(data)

demogr_keys = ["updated_datetime", "created_datetime", "gender", "state", "submitter_id", "year_of_birth", "race", "demographic_id", "ethnicity", "year_of_death"]

diagn_keys = ["classification_of_tumor", "last_known_disease_status", 
        "updated_datetime", "primary_diagnosis", "submitter_id", "tumor_stage", 
        "age_at_diagnosis", "vital_status", "morphology", "days_to_death", 
        "days_to_last_known_disease_status", "days_to_last_follow_up", "state", 
        "days_to_recurrence", "diagnosis_id", "tumor_grade", 
        "tissue_or_organ_of_origin", "days_to_birth", "progression_or_recurrence",
        "prior_malignancy", "site_of_resection_or_biopsy", "created_datetime"]

exposr_keys = ["cigarettes_per_day", "weight", "updated_datetime", 
        "alcohol_history", "alcohol_intensity", "bmi", "years_smoked", 
        "height", "created_datetime", "state", "exposure_id", "submitter_id"]

f = open(infile + '-out.tsv','w')
header = "case_id" + '\t' + '\t'.join(demogr_keys) + '\t' + '\t'.join(diagn_keys) + '\t' + '\t'.join(exposr_keys) + '\n'
f.write(header)


for d in data:
    #print d
    case_id = d["case_id"]
    demogr = [str(d["demographic"][k]) for k in demogr_keys if "demographic" in d]
    diagn = [str(d["diagnoses"][0][k]) for k in diagn_keys if "diagnoses" in d]
    exposr = [str(d["exposures"][0][k]) for k in exposr_keys if "exposures" in d]
    
    entry = case_id  + '\t' + '\t'.join(demogr)  + '\t' + '\t'.join(diagn)  + '\t' + '\t'.join(exposr) + '\n'
    f.write(entry)

    
    #print exposr
    #sys.exit(0)



f.close() 

