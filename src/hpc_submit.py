#!/usr/bin/env python

import os, sys, pathlib

base_path = pathlib.Path(__file__).parent.absolute()
base_path = os.path.dirname(base_path)

src_path = os.path.join(base_path, 'src')
data_path = os.path.join(base_path, 'data')
log_path = os.path.join(data_path, 'log')
out_path = os.path.join(data_path, 'output')

N = int(sys.argv[1])

submit_prefix = 'sbatch --cpus-per-task=1 --time=4:00:00 --partition=ccr,quick --mem=8g'

if not os.path.exists(log_path): os.mkdir(log_path)
if not os.path.exists(out_path): os.mkdir(out_path)

for i in range(N):    
    err_log = os.path.join(log_path, 'error.%d' % i)
    out_log = os.path.join(log_path, 'output.%d' % i)
    
    if os.path.exists(err_log): os.remove(err_log)
    if os.path.exists(out_log): os.remove(out_log)
    
    os.system(' '.join([submit_prefix, '-o', out_log, '-e', err_log , os.path.join(src_path, 'run_wrapper.sh'), str(i), str(N)]))
