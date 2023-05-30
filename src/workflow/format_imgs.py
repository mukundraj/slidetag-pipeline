#!/usr/bin/env python
# print('python', snakemake.input[0])

import os
import glob,os

print('formatting imgs, ip:', snakemake.input[0], 'op:', snakemake.output[0])
print('formatting imgs, ip:', snakemake.input[1], 'op:', snakemake.output[0])




# os.system(f'cp -r {snakemake.input[0]} {snakemake.output[0]}')
# os.system(f'cp -r {snakemake.input[1]} {snakemake.output[0]}')

inp_stags_dir = snakemake.input[0]
fmted_stags_dir = f'{snakemake.output[0]}/stags'

os.system(f'mkdir -p {fmted_stags_dir}')
os.chdir(inp_stags_dir)
for file in glob.glob("*.png"):
    ip_file_path = os.path.join(inp_stags_dir, file)
    op_file_name = file.split('.')[0] + '.tif'
    op_file_path = os.path.join(fmted_stags_dir, op_file_name)
    os.system(f'convert {ip_file_path} {op_file_path}')

inp_nis_dir = snakemake.input[1]
fmted_nis_dir = f'{snakemake.output[0]}/nissls'

os.system(f'mkdir -p {fmted_nis_dir}')
os.chdir(inp_nis_dir)
for file in glob.glob("*.png"):
    ip_file_path = os.path.join(inp_nis_dir, file)
    op_file_name = file.split('.')[0] + '.tif'
    op_file_path = os.path.join(fmted_nis_dir, op_file_name)
    os.system(f'convert {ip_file_path} {op_file_path}')




