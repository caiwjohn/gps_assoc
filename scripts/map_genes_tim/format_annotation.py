import pandas as pd

with open("Ptrichocarpa_210_v3.0_cut_annotation.txt") as f:
    xx=pd.read_csv(f, sep = '\t')
    # print xx

at_TF_family = {}

with open("Ath_TF_list.txt") as f:
    f.readline()
    for line in f:
        line1 = line.strip().split('\t')
        at_TF_family[line1[0]]=line1[2]


pt_TF_family = {}

with open("Ptr_TF_list.txt") as f:
    f.readline()
    for line in f:
        line2 = line.strip().split('\t')
        pt_TF_family[line2[1]]=line2[2]
# print pt_TF_family

xx['at_tf_family'] = xx.Best_hit_arabi_name.map(at_TF_family)
xx['pt_tf_family'] = xx.locusName.map(pt_TF_family)

xr = xx.drop_duplicates(subset='locusName', keep='first')
# print xr

xr.to_csv('format_annotation_with_tf.txt', sep = '\t', header=True, index = False, na_rep='-')

# print xx
