#!/usr/bin/env python
# print('python', snakemake.input[0])

import os

print('formatting imgs, ip:', snakemake.input[0], 'op:', snakemake.output[0])

os.system(f'mkdir -p {snakemake.output[0]}')

os.system(f'cp -r {snakemake.input[0]}/nissls {snakemake.output[0]}')
os.system(f'cp -r {snakemake.input[0]}/stags {snakemake.output[0]}')
