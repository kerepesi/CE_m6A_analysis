import sys
from scipy import stats

expected_m6A_per_bp=0.00038144753045829213

sumLeng={}
summ6A={}
h=1
f=open(sys.argv[1])
for l in f:
    if h==1:
        h=0
    else:
        sl=l.strip().split('|')
        Repeat=sl[1]
        Leng=int(sl[2])
        m6A=int(sl[3])
        Ev=float(sl[11])
#        if Ev < 1e-3:
        if Repeat not in sumLeng:
            sumLeng[Repeat]=0
            summ6A[Repeat]=0
        sumLeng[Repeat]+=Leng
        summ6A[Repeat]+=m6A

print('repeat_group|sum_leng|sum_m6A|sum_m6A_per_bp|is_greater|p')
for Repeat in sorted(sumLeng):
    sum_m6A_per_bp=summ6A[Repeat]/sumLeng[Repeat]
    diff=sum_m6A_per_bp-expected_m6A_per_bp > 0
    if diff==True:
        p=stats.binom_test(summ6A[Repeat], sumLeng[Repeat], expected_m6A_per_bp, alternative='greater')
    else:
        p=stats.binom_test(summ6A[Repeat], sumLeng[Repeat], expected_m6A_per_bp, alternative='less')
    print(Repeat+'|'+str(sumLeng[Repeat])+'|'+str(summ6A[Repeat])+'|'+str(sum_m6A_per_bp)+'|'+str(diff)+'|'+str(p))
