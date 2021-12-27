import sys
import pandas as pd
import numpy as np
import csv
import gzip

anno_file = sys.argv[1] 
anc_file = sys.argv[2] 
phe_file = sys.argv[3] 
out_file = sys.argv[4] 


# read the annotation file
print("read files...")


anno = pd.read_table(anno_file, names=["geneName","name","chrom","strand","txStart","txEnd","cdsStart","cdsEnd","exonCount", "exonStarts", "exonEnds"])
phe = pd.read_csv(phe_file, sep=" ")


with gzip.open(anc_file, mode='rt') as file:
    tsv_file = csv.reader(file, delimiter="\t")
    header = np.array(next(tsv_file))
    common_samples = np.array(list(set(list(map(lambda x: str(x), list(phe.iid)))) & set(header)))
    risk_count = np.zeros((anno.shape[0], len(common_samples)))
    boo = np.isin(header, common_samples)
    print("construct gene-level aggregated count tables...")
    for line in tsv_file:
        in_window = (anno.chrom == ("chr"+line[0])) & (anno.txStart < int(line[1])) & (anno.txEnd > int(line[1]))
        risk_count[np.array(in_window),] = risk_count[np.array(in_window),] + np.array(line)[boo].astype("int")


selected_genes = np.sum(risk_count, axis = 1) != 0

tbl_anno = anno.loc[selected_genes, ].reset_index(drop=True)
tbl_count = pd.DataFrame(risk_count[selected_genes, ], columns = list(common_samples), dtype="int")

print("write table")
output = pd.concat([tbl_anno, tbl_count], axis=1)
output.to_csv(out_file,sep='\t',index=False)

