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
IN_name='Caenorhabditis_elegans.WBcel235.cdna.all.fa'
f=open(IN_name)
for l in f:
    if l[0]=='>':
        sl=l.strip().split()
#        print(sl)
        ID=sl[0][1:]
        Gene=sl[6].split(':')[-1]
        Desc=l.strip()[1:]
        Chrom=sl[2].split(':')[2]
        Start=sl[2].split(':')[3]
        End=sl[2].split(':')[4]
#        print(Gene, Chrom, Start, End)
        gene[ID]=Gene
        desc[ID]=Desc
        chrom[ID]=Chrom.lower()
        start[ID]=int(Start)
        end[ID]=int(End)
        m6A[ID]=0

#f=open('../TE-analysis/n2d1_m6A_context.fa-TopH.bl')
f=open(sys.argv[1])
for l in f:
    if l[0]=='>':
        sl=l.strip().split('|')
        Chrom=sl[0][1:].lower()
        Pos=int(sl[1])
#        print(Chrom, Pos)
        for ID in start:
#            print(chrom[ID],Chrom)
            if chrom[ID]==Chrom:
                if Pos >= start[ID]and Pos <= end[ID]:
                    m6A[ID]+=1

print('ID|gene|length|m6A|m6A/bp|is_greater|p|description')
for ID in sorted(m6A):
    leng=end[ID]-start[ID]
    m6A_per_bp=m6A[ID]/float(leng)
    diff=m6A_per_bp-expected_m6A_per_bp > 0
    if diff==True:
        p=stats.binom_test((m6A[ID]), leng, expected_m6A_per_bp, alternative='greater')
    else:
        p=stats.binom_test((m6A[ID]), leng, expected_m6A_per_bp, alternative='less')
    print(ID+'|'+gene[ID]+'|'+str(leng)+'|'+str(m6A[ID])+'|'+str(m6A_per_bp)+'|'+str(diff)+'|'+str(p)+'|'+desc[ID])
