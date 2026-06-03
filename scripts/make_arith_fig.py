#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Build the 3-panel bridge-constant figure from ../data/Cd_arith_S10.csv (produced
by bridge_constant.py with default MAXK=10). Panels: C(d) vs C0/S(d); C(d) linear
in 1/S(d) through origin (slope C0); product C(d)*S(d) flat at C0.
"""
import csv, numpy as np
import matplotlib; matplotlib.use('Agg')
import matplotlib.pyplot as plt
rows=list(csv.DictReader(open('../data/Cd_arith_S10.csv')))
d=np.array([int(r['d']) for r in rows]); C=np.array([float(r['C_meas']) for r in rows])
S=np.array([float(r['S_adm']) for r in rows]); prod=np.array([float(r['C_times_S']) for r in rows])
pred=np.array([float(r['C0_over_S']) for r in rows]); C0=prod.mean()
fig,axes=plt.subplots(1,3,figsize=(17,4.8))
axes[0].plot(d,C,'o',color='#185FA5',ms=8,label='measured $C(d)$',zorder=4)
axes[0].plot(d,pred,'x',color='#c0392b',ms=9,mew=2,label=r'$C_0/\mathfrak{S}(d)$',zorder=3)
axes[0].set_xlabel('centre-step $d$',fontsize=11); axes[0].set_ylabel('$C(d)$',fontsize=11)
axes[0].set_title('measured $C(d)$ vs closed form $C_0/\\mathfrak{S}(d)$',fontsize=11)
axes[0].legend(fontsize=10); axes[0].grid(alpha=.25)
invS=1/S
axes[1].plot(invS,C,'o',color='#185FA5',ms=8,zorder=3)
xs=np.linspace(0,invS.max()*1.05,50); axes[1].plot(xs,C0*xs,'--',color='#c0392b',lw=1.8,label=f'slope $C_0={C0:.2f}$')
axes[1].set_xlabel(r'$1/\mathfrak{S}(d)$',fontsize=11); axes[1].set_ylabel('$C(d)$',fontsize=11)
axes[1].set_title('$C(d)$ linear in $1/\\mathfrak{S}(d)$ through origin',fontsize=11)
axes[1].legend(fontsize=10); axes[1].grid(alpha=.25)
axes[2].plot(d,prod,'D-',color='#2ca25f',lw=1.5,ms=6)
axes[2].axhline(C0,color='gray',ls=':',lw=1.3,label=f'$C_0={C0:.2f}=1/\\rho$')
axes[2].set_ylim(C0*0.95,C0*1.05)
axes[2].set_xlabel('centre-step $d$',fontsize=11); axes[2].set_ylabel(r'$C(d)\cdot\mathfrak{S}(d)$',fontsize=11)
axes[2].set_title(f'product constant (CV {100*prod.std()/prod.mean():.2f}%)',fontsize=11)
axes[2].legend(fontsize=10); axes[2].grid(alpha=.25)
plt.suptitle('Bridge constant in $S_{10}$:  $C(d)=C_0/\\mathfrak{S}(d)$, $C_0=1/\\rho\\approx62.75$ (inverse twin density)',fontsize=12.5,y=1.02)
plt.tight_layout()
plt.savefig('fig_paper8_arith.pdf',bbox_inches='tight')
plt.savefig('fig_paper8_arith.png',dpi=160,bbox_inches='tight')
print("figure saved")
