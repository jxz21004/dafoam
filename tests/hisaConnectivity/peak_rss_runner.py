#!/usr/bin/env python3
"""Measure the solver Python process itself, including its loaded native libraries."""
from pathlib import Path
import resource
import runpy
import sys
runner=Path(sys.argv[1]).resolve()
args=sys.argv[2:]
case=Path(args[args.index('--case')+1]).resolve()
sys.argv=[str(runner)]+args
try:
    runpy.run_path(str(runner),run_name='__main__')
finally:
    (case/'peak_rss_kb.txt').write_text(str(resource.getrusage(resource.RUSAGE_SELF).ru_maxrss)+'\n')
