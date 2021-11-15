import pandas as pd
import os,collections

gene_df_all = pd.DataFrame(columns=['gene'])
for file in os.listdir('./'):
    if file.endswith('VFDB_blastp.txt'):
        sample_id=file[:-16]
    #    print(sample_id)
        with open('./%s_58vf_blastp.txt'%sample_id,'r') as g:
            results={}
            for line in g.readlines():
                n=line.split('\t')
                if float(n[5])>=80 and float(n[12])>=70:
                    results[n[0]]=n[1]                
        gene=[]
        with open('./%s'%file,'r') as f:
            for line in f.readlines():
                i=line.split('\t')
                if float(i[5])>=80 and float(i[12])>=70:
                    if i[0] not in results:
                        gene.append(i[1])
                    else:
                        print(results[i[0]])
            gene_counter=collections.Counter(gene)
            gene_df=pd.DataFrame(gene_counter.items(), columns=['gene', '%s'%sample_id])        
        gene_df_all=gene_df_all.merge(gene_df,how='outer',on='gene').fillna(0)
gene_df_all
with open('../52gennomes_VFDB_count1.txt','w') as f:
    gene_df_all.to_csv(f, header = True, index = False, sep = "\t")
history
