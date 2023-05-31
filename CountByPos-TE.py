import sys
from scipy import stats

expected_m6A_per_bp=0.00038144753045829213
ID={}
desc={}
gene={}
chrom={}
start={}
end={}
m6A={}
ID=0
f=open('ce10_dfam.nrph.hits')
for l in f:
    if l[0]!='#':
        ID+=1
        sl=l.strip().split()
#        print(sl)
        Strand=sl[8] 
        Gene=sl[2]
        Desc=l.strip()
        Chrom=sl[0][-1]
        Ev=float(sl[4])
        if Strand == '+':
            Start=sl[11]
            End=sl[12]
        else:
            Start=sl[12]
            End=sl[11]
#        print(Gene, Chrom, Start, End)
        gene[ID]=Gene
        desc[ID]=Desc
        chrom[ID]=Chrom.lower()
        start[ID]=int(Start)
        end[ID]=int(End)
        m6A[ID]=0
    else:
        header=l.strip()[1:]

IN_name=sys.argv[1]
f=open(IN_name)
for l in f:
    if l[0]=='>':
        sl=l.strip().split('|')
        Chrom=sl[0][1:].lower()
        Pos=int(sl[1])
#        print(Chrom, Pos)
        for ID in start:
            if chrom[ID]==Chrom:
                if Pos >= start[ID]and Pos <= end[ID]:
                    m6A[ID]+=1

print('ID|TE|length|m6A|m6A/bp|is_greater|p|'+'|'.join(header.split('\t')))
for ID in sorted(m6A):
    leng=end[ID]-start[ID]
    m6A_per_bp=m6A[ID]/float(leng)
    diff=m6A_per_bp-expected_m6A_per_bp > 0
    if diff==True:
        p=stats.binom_test((m6A[ID]), leng, expected_m6A_per_bp, alternative='greater')
    else:
        p=stats.binom_test((m6A[ID]), leng, expected_m6A_per_bp, alternative='less')
    Description='|'.join(desc[ID].split('\t'))
    print(str(ID)+'|'+gene[ID]+'|'+str(leng)+'|'+str(m6A[ID])+'|'+str(m6A_per_bp)+'|'+str(diff)+'|'+str(p)+'|'+Description)
