#!/usr/bin/env python

import os, sys, pathlib

base_path = pathlib.Path(__file__).parent.absolute()
base_path = os.path.dirname(base_path)

data_path = os.path.join(base_path, 'data')
src_path = os.path.join(base_path, 'src')

N = int(sys.argv[1])

submit_prefix = 'sbatch --cpus-per-task=1 --time=4:00:00 --partition=ccr,quick --mem=8g'

for i in range(N):    
    err_log = os.path.join(data_path, 'log', 'error.%d' % i)
    out_log = os.path.join(data_path, 'log', 'output.%d' % i)
    
    if os.path.exists(err_log): os.remove(err_log)
    if os.path.exists(out_log): os.remove(out_log)
    
    os.system(' '.join([submit_prefix, '-o', out_log, '-e', err_log , os.path.join(src_path, 'correlation.py'), str(i), str(N)]))
