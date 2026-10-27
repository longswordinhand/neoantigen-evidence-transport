#!/usr/bin/env python3
from pathlib import Path
import numpy as np, pandas as pd
from sklearn.metrics import roc_auc_score, average_precision_score
ROOT=Path(__file__).resolve().parents[1]
D=pd.read_csv(ROOT/'data/benchmark_features_only.csv')
METHODS=['MHCflurry2.2-Affinity','MHCflurry2.2-Presentation','MixMHCpred3.0','PRIME2.1','BigMHC-EL']
def macro(g,m):
    a=[]
    for _,z in g.groupby(['trial','patient_group']):
        if z.response.nunique()==2: a.append(roc_auc_score(z.response,z[m]))
    return float(np.mean(a)) if a else np.nan, len(a)
rows=[]
for scope,g in list(D.groupby('trial'))+[('POOLED',D)]:
  for m in METHODS:
    ma,n=macro(g,m)
    rows.append(dict(scope=scope,method=m,auc=roc_auc_score(g.response,g[m]),ap=average_precision_score(g.response,g[m]),macro_patient_auc=ma,eligible_patients=n,n=len(g),positive=int(g.response.sum()),patients=g[['trial','patient_group']].drop_duplicates().shape[0]))
out=pd.DataFrame(rows)
out.to_csv(ROOT/'results/reproduced_metrics.csv',index=False)
ref=pd.read_csv(ROOT/'results/reference_metrics.csv')
ref=ref[ref.method.isin(METHODS)].sort_values(['scope','method']).reset_index(drop=True)
chk=out.sort_values(['scope','method']).reset_index(drop=True)
cols=['auc','ap','macro_patient_auc']
mx=max(float(np.nanmax(np.abs(chk[c]-ref[c]))) for c in cols)
print(out.to_string(index=False))
print(f'\nmax_abs_metric_difference_vs_reference={mx:.12g}')
if mx>1e-10: raise SystemExit('Metric reproduction mismatch')
print('PASS: primary point metrics reproduce reference outputs.')
