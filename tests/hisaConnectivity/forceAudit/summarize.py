#!/usr/bin/env python3
"""Project independently integrated surface forces and compare to DAFoam."""
import argparse
import json
from pathlib import Path
import numpy as np
p=argparse.ArgumentParser();p.add_argument('--case',type=Path,required=True);p.add_argument('--reference',type=Path,required=True);a=p.parse_args()
f=np.genfromtxt(a.case/'force_audit_faces.csv',delimiter=',',names=True)
q=13402.323345186078; area=45.5; alpha=np.deg2rad(2.)
d=np.array([np.cos(alpha),0,np.sin(alpha)]);l=np.array([-np.sin(alpha),0,np.cos(alpha)])
sf=np.column_stack([f['s'+x] for x in 'xyz']);pn=np.column_stack([f['fp'+x] for x in 'xyz']);tv=np.column_stack([f['fv'+x] for x in 'xyz'])
r={'q':q,'reference_area':area,'alpha_deg':2.,'faces':len(f),'area_vector':sf.sum(axis=0).tolist()}
for name,v in [('pressure',pn.sum(axis=0)),('viscous',tv.sum(axis=0)),('total',(pn+tv).sum(axis=0)),('gauge_pressure',(pn-26499.89786492552*sf).sum(axis=0))]:
 r[name]={'force':v.tolist(),'cd':float(v@d/(q*area)),'cl':float(v@l/(q*area))}
ref=json.loads(a.reference.read_text())['functions'];r['reference']=ref
r['difference']={k:r['total'][k]-ref[k] for k in ['cd','cl']}
cp=(f['p']-26499.89786492552)/q
r['cp_range']=[float(cp.min()),float(cp.max())]
r['pressure_drag_fraction']=r['pressure']['cd']/r['total']['cd']
r['largest_cp_faces']=[{'face_index':int(i),'xyz':[float(f[x][i]) for x in 'xyz'],'cp':float(cp[i])} for i in np.argsort(cp)[-5:]]
(a.case/'summary.json').write_text(json.dumps(r,indent=2)+'\n');print(json.dumps(r,indent=2))
