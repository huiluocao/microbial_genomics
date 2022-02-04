import os
p_contigs=[]
for file in os.listdir('./'):
    if file.endswith('plasmidfinder_blastn.txt'):
        sample=file[:-25]
        print(sample)
        with open(file,'r') as f:
            for line in f.readlines():
                i=line.split('\t')
                if float(i[2])>=90:
                    print(sample,i[0],i[1])
                    contig='%s:%s'%((i[1]).split('_')[0],i[0])
                    p_contigs.append([sample,contig])
p_contigs
for key, val in p_contigs:
    p_contigs_dict.setdefault(key,[]).append(val)
for k,v in p_contigs_dict.items():
    with open('../all_inc_type.txt','a') as f:
        f.write(k+'\t'+';'.join(v)+'\n')
