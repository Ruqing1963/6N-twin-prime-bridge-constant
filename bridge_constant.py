#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Arithmetic structure of the bridge constant C(d)  (Part VIII).

Measures, on S_K twin centres, the bridge constant
    C(d) = r(d|omega) / P(N+d twin|omega)   (omega-independent, Part VI)
for d=1..DMAX, computes the Hardy-Littlewood admissibility S(d) of the separation
(pure CRT singular-series factor), and verifies

    C(d) * S(d) = C0   (constant in d)                       [arithmetic structure]
    C0 = 1/rho,  rho = twin-centre line density               [C0 identified]

via the omega-merged two-point relation P_merged(d) = rho * S(d). Emits
Cd_arith_S{K}.csv with columns d,6d,C_meas,S_adm,C_times_S,C0_over_S.

  S(d) = prod_{q>3} [ (#{r: r,r+d both q-safe}/q) / ((q-2)/q)^2 ],
         q-safe meaning r not in dead(q)={+-6^{-1} mod q}.

USAGE: python bridge_constant.py   (default S10, ~12 min); MAXK=9 for S9.
Requires: numpy.
"""
import numpy as np, math, os, csv
from collections import defaultdict
def primes_upto(n):
    s=np.ones(n+1,bool); s[:2]=False
    for i in range(2,int(math.isqrt(n))+1):
        if s[i]: s[i*i::i]=False
    return np.nonzero(s)[0].astype(np.int64)
MAXK=int(os.environ.get("MAXK",10))
LO=10**(MAXK-1)//6+1; HI=10**MAXK//6; SEG=4_000_000
PB=int(math.isqrt(6*HI+250))+1; BP=primes_upto(PB)
DMAX=30
N_list=[]; om_list=[]
n=LO; import time; t0=time.time()
while n<=HI:
    nh=min(n+SEG,HI+1); sz=nh-n
    rem=np.arange(n,nh,dtype=np.int64); ob=np.zeros(sz,np.int16)
    for p in BP:
        if p*p>nh-1: break
        f=((n+p-1)//p)*p
        if f>=nh: continue
        idx=np.arange(f-n,sz,p)
        if idx.size==0: continue
        sub=rem[idx]; m=(sub%p)==0
        while m.any(): sub[m]//=p; m=(sub%p)==0
        rem[idx]=sub
        if p>3: ob[idx]+=1
    ob[rem>1]+=1
    Narr=np.arange(n,nh,dtype=np.int64)
    vlo=6*n-1; vhi=6*(nh-1)+1; span=vhi-vlo+1
    comp=np.zeros(span,bool); sq=int(math.isqrt(vhi))+1
    for p in BP:
        if p>sq: break
        st=max(p*p,((vlo+p-1)//p)*p)
        if st>vhi: continue
        comp[st-vlo:span:p]=True
    tw=(~comp[(6*Narr-1)-vlo])&(~comp[(6*Narr+1)-vlo])
    pos=np.nonzero(tw)[0]
    N_list.append(Narr[pos]); om_list.append(ob[pos])
    n=nh
N_arr=np.concatenate(N_list); om_arr=np.concatenate(om_list).astype(np.int16)
ntw=len(N_arr); Nrange=HI-LO+1; rho=ntw/Nrange
print(f"S{MAXK} twins {ntw:,}; scan {time.time()-t0:.0f}s; twin line density rho={rho:.6f}")
def is_twin(vals):
    idx=np.searchsorted(N_arr,vals); idx=np.clip(idx,0,len(N_arr)-1)
    return N_arr[idx]==vals
def primes_list(n):
    s=[True]*(n+1); s[0]=s[1]=False
    for i in range(2,int(n**.5)+1):
        if s[i]:
            for j in range(i*i,n+1,i): s[j]=False
    return [i for i in range(5,n+1) if s[i]]
def adm(dd,qmax=5000):
    p=1.0
    for q in primes_list(qmax):
        inv=pow(6,q-2,q); dead={inv%q,(-inv)%q}
        safe=sum(1 for r in range(q) if r not in dead and (r+dd)%q not in dead)/q
        p*=safe/((q-2)/q)**2
    return p
GMAX=DMAX+1
g=np.diff(N_arr); omL=om_arr[:-1]
gap_count=defaultdict(lambda: np.zeros(GMAX+1))
for i in range(len(g)):
    gg=g[i]
    if 1<=gg<=GMAX: gap_count[int(omL[i])][gg]+=1
overall=np.zeros(GMAX+1)
for om in gap_count: overall+=gap_count[om]
base_pref=overall/overall[1:GMAX+1].sum()
omegas=[om for om in range(1,7) if (om_arr==om).sum()>=20000 and gap_count[om][1:GMAX+1].sum()>=5000]
print(f"omegas used: {omegas}")
Cvals={}
for d in range(1,DMAX+1):
    ratios=[]
    for om in omegas:
        tot=gap_count[om][1:GMAX+1].sum(); cnt=gap_count[om][d]
        if cnt<200: continue
        r=(cnt/tot)/base_pref[d]
        surv=is_twin(N_arr[om_arr==om]+d).mean()
        if surv>0: ratios.append(r/surv)
    if len(ratios)>=3: Cvals[d]=np.mean(ratios)
A={d:adm(d) for d in Cvals}
prod=np.array([Cvals[d]*A[d] for d in sorted(Cvals)])
C0=prod.mean()
print(f"\nC(d)*S(d): mean C0={C0:.3f}, CV={100*prod.std()/prod.mean():.2f}%")
print(f"C0 = {C0:.2f}  vs  1/rho = {1/rho:.2f}  (agree {100*(1-abs(1/rho-C0)/C0):.1f}%)")
print(f"{'d':>3}{'C(d)':>9}{'S(d)':>9}{'C*S':>9}")
for d in sorted(Cvals):
    print(f"{d:>3}{Cvals[d]:>9.2f}{A[d]:>9.3f}{Cvals[d]*A[d]:>9.2f}")
with open(f'Cd_arith_S{MAXK}.csv','w',newline='') as f:
    w=csv.writer(f); w.writerow(['d','6d','C_meas','S_adm','C_times_S','C0_over_S'])
    for d in sorted(Cvals):
        w.writerow([d,6*d,f'{Cvals[d]:.3f}',f'{A[d]:.4f}',f'{Cvals[d]*A[d]:.3f}',f'{C0/A[d]:.3f}'])
print(f"\n[ok] wrote Cd_arith_S{MAXK}.csv ; C0=1/rho identity: C0={C0:.2f}, 1/rho={1/rho:.2f}")
