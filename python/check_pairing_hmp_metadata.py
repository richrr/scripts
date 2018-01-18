import sys
import os
from utils import * 


patterns = read_file(sys.argv[1])

skin = ["Left_Antecubital_fossa", "Left_Retroauricular_crease", "Right_Antecubital_fossa", "Right_Retroauricular_crease"]
allowed = ["Stool", "Posterior_fornix"]
allowed.extend(skin)
print allowed

rsid_site_dict = dict() # patient and actual site hmpbodysubsite
rsid_tissue_dict = dict() # patient and tissue
rsid_psn_dict = dict() # patient and psn
rsid_sampid_dict = dict() # patient and sample id

map = ['\t'.join(["sampid", "patient", "psn", "site", "tissue"])]

for p in patterns:
    p = p.strip()
    conts = p.split('\t')
    
    # keep only females and visit 1
    if len(conts) != 15 or conts[11] != 'female' or conts[14] != '1':
        continue
    
    sampid = conts[0]
    patient = conts[1]
    psn = conts[2]
    site = conts[12]
    
    
    # keep only those sites that we are interested to study
    if site not in allowed:
        continue
    
    print p
    #print conts

    
    tissue = ''
    if site in skin:
        tissue = 'Skin'
    elif site == 'Posterior_fornix':
        tissue = 'Vagina'
    elif site == 'Stool':
        tissue = 'Stool'
    
    if patient in rsid_site_dict:
        rsid_site_dict[patient].append(site)
        rsid_tissue_dict[patient].append(tissue)
        rsid_psn_dict[patient].append(psn)
        rsid_sampid_dict[patient].append(sampid)
    else:
        rsid_site_dict[patient] = [site]
        rsid_tissue_dict[patient] = [tissue]
        rsid_psn_dict[patient] = [psn]
        rsid_sampid_dict[patient] = [sampid]

    map.append('\t'.join([sampid, patient, psn, site, tissue]))    
    

writeLIST_to_file(map, "skin-stool-vagina-map.txt")
writeDICT_to_file(rsid_site_dict, "rsid_site_dict.txt")
writeDICT_to_file(rsid_tissue_dict, "rsid_tissue_dict.txt")
writeDICT_to_file(rsid_psn_dict, "rsid_psn_dict.txt")
writeDICT_to_file(rsid_sampid_dict, "rsid_sampid_dict.txt")

choose_paired_samples = ["RSID"]
choose_st_vag_samples = ["RSID"]
# select the samples that have skin, stool, and vagina samples
for k,v in rsid_tissue_dict.items():
    if not k:
        continue
    if "Skin" in v and "Stool" in v and "Vagina" in v:
        print k
        choose_paired_samples.append(k)
    if "Stool" in v and "Vagina" in v:
        choose_st_vag_samples.append(k)

        
writeLIST_to_file(choose_paired_samples, "samples-w-skin-stool-vagina.txt")
writeLIST_to_file(choose_st_vag_samples, "samples-w-stool-vagina.txt")
