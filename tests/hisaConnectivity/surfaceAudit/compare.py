#!/usr/bin/env python3
"""Compare saved CGNS cell-centred wall Cp to Foam face data."""
import ctypes as c
import json
from pathlib import Path
import numpy as np
from scipy.spatial import cKDTree
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
root=Path('/mdouser/test_aerostructure')
out=Path(__file__).resolve().parent
lib=c.CDLL('/home/dafoamuser/dafoam/packages/CGNS-4.5.0/opt-gfortran/lib/libcgns.so')
fn=c.c_int()
def check(code):
    if code: raise RuntimeError('CGNS read failed: '+str(code))
check(lib.cg_open(str(root/'reports/rigid_comparison/adflow_tight/rigid2_000_surf.cgns').encode(),0,c.byref(fn)))
nz=c.c_int();check(lib.cg_nzones(fn,1,c.byref(nz)))
centres=[];pressures=[];zones=[]
for z in range(1,nz.value+1):
    name=c.create_string_buffer(33);size=(c.c_int*9)();check(lib.cg_zone_read(fn,1,z,name,size))
    if not name.value.startswith(b'NSWall'):continue
    ni,nj=size[0],size[1];lo=(c.c_int*2)(1,1);hi=(c.c_int*2)(ni,nj)
    coords=[]
    for ax in 'XYZ':
        arr=np.zeros(ni*nj);check(lib.cg_coord_read(fn,1,z,('Coordinate'+ax).encode(),4,lo,hi,arr.ctypes.data_as(c.c_void_p)))
        coords.append(arr.reshape((ni,nj),order='F'))
    xyz=np.stack(coords,axis=-1);ctr=(xyz[:-1,:-1]+xyz[1:,:-1]+xyz[:-1,1:]+xyz[1:,1:])/4
    a=np.zeros((ni-1)*(nj-1));hi=(c.c_int*2)(ni-1,nj-1)
    check(lib.cg_field_read(fn,1,z,1,b'CoefPressure',4,lo,hi,a.ctypes.data_as(c.c_void_p)))
    verts=[xyz[:-1,:-1],xyz[1:,:-1],xyz[1:,1:],xyz[:-1,1:]]
    accum=np.zeros_like(ctr);weights=np.zeros(ctr.shape[:-1])
    for k in range(4):
        v0,v1=verts[k],verts[(k+1)%4]
        weight=np.linalg.norm(np.cross(v0-ctr,v1-ctr),axis=-1)
        accum+=weight[...,None]*(v0+v1+ctr)/3;weights+=weight
    ctr=accum/weights[...,None]
    centres.extend(ctr.reshape((-1,3),order='F'));pressures.extend(a);zones.extend([name.value.decode()]*len(a))
lib.cg_close(fn)
xyz=np.array(centres);cpad=np.array(pressures)
f=np.genfromtxt(root/'reports/rigid_comparison/hisa_adiabatic_lowre/force_audit_faces.csv',names=True,delimiter=',')
hxyz=np.column_stack([f[x] for x in 'xyz']);dist,idx=cKDTree(xyz).query(hxyz)
print('Maximum face-centre mismatch:',max(dist))
assert len(set(idx))==len(idx), 'Face matching is not one-to-one'
assert max(dist)<0.01, 'Face centres differ too much'
q=13402.323345186078;S=45.5;alpha=np.deg2rad(2);d=np.array([np.cos(alpha),0,np.sin(alpha)])
sf=np.column_stack([f['s'+x] for x in 'xyz']);area=np.linalg.norm(sf,axis=1)
cp=(f['p']-26499.89786492552)/q;delta=cp-cpad[idx]
pressure_delta=delta*(sf@d)/S
fv=np.column_stack([f['fv'+x] for x in 'xyz']);normal=sf/area[:,None]
shear=fv/area[:,None];tangent=shear-(shear*normal).sum(axis=1)[:,None]*normal
cf=np.linalg.norm(tangent,axis=1)/q
rows=np.column_stack([hxyz,dist,cpad[idx],cp,delta,pressure_delta,cf])
np.savetxt(out/'matched_faces.csv',rows,delimiter=',',header='x,y,z,match_distance,Cp_ADflow,Cp_HiSA,delta_Cp,delta_CD_pressure,Cf_HiSA',comments='')
# Relative chord coordinate from actual per-strip face-centre bounds, for locating regions only.
y=np.round(hxyz[:,1],6);xc=np.zeros(len(f))
for yy in np.unique(y):
 mask=y==yy;x=hxyz[mask,0];xc[mask]=(x-x.min())/max(np.ptp(x),1e-12)
fig,axes=plt.subplots(2,3,figsize=(14,8),layout='constrained')
upper=sf[:,2]<0
for row,mask in enumerate([upper,~upper]):
 for col,(val,title) in enumerate([(cpad[idx],'ADflow Cp'),(cp,'HiSA Cp (not converged)'),(delta,'HiSA - ADflow Cp')]):
  ax=axes[row,col];sc=ax.scatter(hxyz[mask,1],hxyz[mask,0],c=val[mask],s=40,cmap='coolwarm',vmin=-1.7 if col<2 else -1.5,vmax=2.7 if col<2 else 1.5)
  ax.set(xlabel='span y [m]',ylabel='x [m]',title=title+(' / Sf_z < 0' if row==0 else ' / Sf_z >= 0'));fig.colorbar(sc,ax=ax)
fig.savefig(out/'surface_maps.png',dpi=180);plt.close(fig)
fig,axes=plt.subplots(1,2,figsize=(11,4),layout='constrained')
for ax,mask in zip(axes,[upper,~upper]):
 sc=ax.scatter(hxyz[mask,1],hxyz[mask,0],c=cf[mask],s=40,cmap='viridis',vmin=0,vmax=max(cf));fig.colorbar(sc,ax=ax);ax.set(xlabel='span y [m]',ylabel='x [m]',title='HiSA tangential Cf / '+('Sf_z < 0' if mask is upper else 'Sf_z >= 0'))
fig.savefig(out/'hisa_friction_map.png',dpi=180);plt.close(fig)
fig,axes=plt.subplots(1,2,figsize=(12,4),layout='constrained')
upper=sf[:,2]<0 # outward fluid-domain normals point into solid
for ax,mask,title in [(axes[0],upper,'Upper-facing surface (Sf_z < 0)'),(axes[1],~upper,'Other faces (Sf_z >= 0)')]:
 ax.scatter(cpad[idx][mask],cp[mask],s=22);ax.plot([-2,3],[-2,3],'k--');ax.set(xlabel='ADflow Cp',ylabel='HiSA Cp',title=title)
fig.savefig(out/'cp_comparison.png',dpi=180);plt.close(fig)
order=np.argsort(abs(pressure_delta))[::-1]
r={'faces':len(f),'max_matching_distance':float(max(dist)),'unique_matches':len(set(idx)),'adflow_cp_range':[float(cpad.min()),float(cpad.max())],'hisa_cp_range':[float(cp.min()),float(cp.max())],'mapped_adflow_pressure_cd':float(sum(cpad[idx]*(sf@d)/S)),'hisa_pressure_cd':float(sum(cp*(sf@d)/S)),'pressure_cd_difference':float(sum(pressure_delta)),'largest_pressure_difference_faces':[{'face':int(i),'xyz':hxyz[i].tolist(),'ADflow_Cp':float(cpad[idx[i]]),'HiSA_Cp':float(cp[i]),'delta_CD':float(pressure_delta[i]),'ADflow_zone':zones[idx[i]]} for i in order[:15]],'top10_absolute_contribution_fraction':float(sum(abs(pressure_delta[order[:10]]))/sum(abs(pressure_delta)))}
r['pressure_difference_by_adflow_zone']={z:float(sum(pressure_delta[np.array(zones)[idx]==z])) for z in sorted(set(zones))}
r['highest_hisa_cf_faces']=[{'face':int(i),'xyz':hxyz[i].tolist(),'Cf':float(cf[i])} for i in np.argsort(cf)[-5:][::-1]]
(out/'summary.json').write_text(json.dumps(r,indent=2)+'\n');print(json.dumps(r,indent=2))
