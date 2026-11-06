#!/usr/bin/env python3
"""Prepare an isolated continuation without altering the source case."""
import argparse,re,shutil
from pathlib import Path
p=argparse.ArgumentParser();p.add_argument('--source',type=Path,required=True);p.add_argument('--time',required=True);p.add_argument('--output',type=Path,required=True);p.add_argument('--cfl',type=float,required=True);p.add_argument('--steps',type=int,required=True);p.add_argument('--gauss-gradient',action='store_true');a=p.parse_args()
a.output.mkdir(parents=True,exist_ok=False)
for name in ['constant','system']:shutil.copytree(a.source/name,a.output/name)
shutil.copytree(a.source/a.time,a.output/'0',ignore=shutil.ignore_patterns('uniform','pseudoCo*','rhoR','rhoUR','rhoER','C','Cx','Cy','Cz'))
f=a.output/'system/controlDict';s=f.read_text()
for key,value in [('startFrom','startTime'),('startTime','0'),('endTime',str(a.steps)),('writeInterval','500'),('purgeWrite','0')]:
 s,n=re.subn(r'(?m)^'+key+r'\s+[^;]+;',key+' '+value+';',s);assert n==1
f.write_text(s)
f=a.output/'system/fvSolution';s=f.read_text()
for key,value in [('pseudoCoNum',str(a.cfl)),('pseudoCoNumMax',str(a.cfl))]:
 s,n=re.subn(r'(?m)^(\s*)'+key+r'\s+[^;]+;',lambda m:m[1]+key+' '+value+';',s);assert n==1
for key,value in [('resetPseudo','true'),('pseudoCoNumMaxIncreaseFactor','1'),('pseudoCoNumMinDecreaseFactor','1')]:
 if re.search(r'\b'+key+r'\s+',s):s=re.sub(r'\b'+key+r'\s+[^;]+;',key+' '+value+';',s)
 else:s=s.replace('pseudoTime\n{','pseudoTime\n{\n    '+key+' '+value+';')
f.write_text(s)
if a.gauss_gradient:
 f=a.output/'system/fvSchemes';s=f.read_text();assert 'faceLeastSquares linear' in s;s=s.replace('faceLeastSquares linear','Gauss linear');f.write_text(s)
