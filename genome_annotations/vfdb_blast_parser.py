import pandas as pd
import collections
import os

gene_df_all = pd.DataFrame(columns=['gene'])
for file in os.listdir('./'):
    if file.endswith('blastp_vfdb.txt'):
        sample_id=file[:-16]
        print(sample_id)
        gene=[]
        with open('./%s'%file,'r') as f:
            for line in f.readlines():
                i=line.split('\t')
                if float(i[5])>=80 and float(i[12])>=70:
                    gene.append(i[1])
            gene_counter=collections.Counter(gene)
            gene_df=pd.DataFrame(gene_counter.items(), columns=['gene', '%s'%sample_id])        
        gene_df_all=gene_df_all.merge(gene_df,how='outer',on='gene').fillna(0)
gene_df_all
with open('../90vf_count.txt','w') as f:
    gene_df_all.to_csv(f, header = True, index = False, sep = "\t")
