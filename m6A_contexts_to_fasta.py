
uniq_cont={}
f=open('BAM_d5/m54282_200107_024536.subreads.bam-pbalign.bam-ipdSummary-minCoverage1.gff')
for l in f:
    if l[0] != '#':
        sl=l.strip().split()
        chrom=sl[0]
        mod=sl[2]
        pos=int(sl[3])
        context=sl[8].split(';')[1].split('=')[1]
        if mod=='m6A':
            if chrom not in uniq_cont:
                uniq_cont[chrom]={} 
            uniq_cont[chrom][pos]=context

for chrom in sorted(uniq_cont):
    for pos in sorted(uniq_cont[chrom]):
        print('>'+chrom+'|'+str(pos))
        print(uniq_cont[chrom][pos])
