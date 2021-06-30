import os, collections
import pandas as pd

iall_arg_df=pd.DataFrame(columns=['arg'])
print(all_arg_df)
for file in os.listdir('./'):
    if file.endswith('fna_blastn.txt'):
        sample=file[:-25]
        print(sample)
        with open('./%s'%file,'r') as f:
            args=[]
            for line in f.readlines():
                line=line.strip()
                i=line.split('\t')
                if float(i[2])>=80 and float(i[11])>=100:
                    args.append((i[1]).split('_')[0])
            arg_counter=collections.Counter(args)
            arg_df=pd.DataFrame(arg_counter.items(), columns=['arg', '%s'%sample])
            print(arg_df)
            all_arg_df=all_arg_df.merge(arg_df,how='outer',on='arg').fillna(0)
    print(all_arg_df)
with open('../100genomes_resfinder.txt','w') as f:
    all_arg_df.to_csv(f, header = True, index = False, sep = "\t")
