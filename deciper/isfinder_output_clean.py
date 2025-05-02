import os

for file in os.listdir('./'):
    if file.endswith('_clean.txt'):
        basef = os.path.basename(file)
        sample_id = basef.split('_')[0]
        print(sample_id)
        df = pd.read_table(file,sep = '\t', engine = 'python', header = None)
        df.columns = ['0','1','2','3','4','5','6','7','8','9','10','11']
        grouped = df.groupby('0')
        keep_blast=[]        
        for name, group in grouped:
            print(name)
            group=group.values.tolist()
            if len(group) == 1:
                keep_blast.append(group[0])
            while len(group) > 1:
                print(group[0])
                keep_blast.append(group[0])
                items = group[0]
                s = items[6]
                e = items[7]
                score = items[11]
                index_keep = []
                index_del = [0]
                for i in range(1,len(group)):
                    #print(i)
                    elements=group[i]
                    #del group[0]
                    if int(elements[6]) >= int(s) and int(elements[6]) < int(e) and int(elements[7]) <= int(e):
                        index_del.append(i)
                    elif  int(elements[6]) < int(e) and elements[6] >= s and int(elements[7]) > int(e):
                        if elements[11] < score:
                            index_del.append(i)
                        else:
                            index_keep.append(i)
                            #print(group[i])
                            index_del.append(i)
                    elif elements[6] <= s and s <= elements[7] <= e:
                        if elements[11] < score:
                            index_del.append(i)
                        else:
                            #print(group[i])
                            index_keep.append(i)
                            index_del.append(i)
                    else:
                        pass
                print(index_del)
                print(index_keep)
                if len(index_keep) >= 1:
                    for i in index_keep:
                        keep_blast.append(group[i])
                #for n in index_del:
                group = [i for i in group if group.index(i) not in index_del]
                print(len(group))
    with open('./%s_ISfinder_clean_final.txt'%sample_id,'w') as f:
        for i in keep_blast:
            f.write('\t'.join(str(j) for j in i)+'\n')
