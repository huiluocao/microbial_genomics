p_contigs=[]
p_contigs_dict={}
for file in os.listdir('./'):
    if os.path.isdir(file):
        print(file)
        try:
            with open('./%s/inctyper.csv'%file) as f:
                for line in f.readlines()[1:]:
                    i=line.split(',')
                    contig='%s:%s'%(((i[5]).strip('"')).rsplit('_', 1)[0],(i[2]).strip('"'))
                    print((i[1]).strip('"'),contig)
                    p_contigs.append([(i[1]).strip('"'),contig])
        except:
            pass
            
for key, val in p_contigs:
    p_contigs_dict.setdefault(key,[]).append(val)
for k,v in p_contigs_dict.items():
    with open('./all_inc_type.txt','a') as f:
        f.write(k+'\t'+';'.join(v)+'\n')
