import os

for file in os.listdir('./'):
    if file.endswith('serotypes__genes__EcOH__results.txt'):
        sample=file.split('_serotypes')[0]
        #print(sample)
        with open(file,'r') as f:
            O=[]
            for line in f.readlines()[1:]:
                #print(line)
                i=line.split('\t')
                #print(i[1])
                H=((i[1])[5:]).split('_')[0]
                for n in range(2,len(i)):
                    #print(n)
                    O.append(((i[n])[4:]).split('_')[0])
                O=(list(set(O)))
        with open('../serotype_srst_summary.txt','a') as g:
            g.write(sample+'\t'+H+'\t'+'/'.join(O)+'\n')
