quast_dfs=[]
for file in os.listdir('./'):
    if file.endswith('quast'):
        sample=file[:-6]
        print(sample)
        quast_df=pd.read_table('./%s/report.tsv'%file,sep = '\t',engine = 'python', header = 0)
        print(quast_df.columns)
        quast_df=quast_df.set_index('Assembly')
        quast_dfs.append(quast_df)
    quast_df_all=pd.concat(quast_dfs, axis=1)
quast_df_all
with open('../quast_report_all.txt','w') as f:
    quast_df_all.to_csv(f, header = True, index = True, sep = "\t")
