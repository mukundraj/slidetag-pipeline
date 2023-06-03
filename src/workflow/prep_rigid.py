from pathlib import Path
import sys
path_root = Path(__file__).parents[2]
sys.path.append(str(path_root))

import os
from src.python.slicer_mrml_gen import get_sub_text

print ('inpy', {snakemake.input.rigid_imgs_dir})

stag_imgs_dir = f'{snakemake.input.rigid_imgs_dir}/stags'

# read filenames in stag_imgs_dir
stag_imgs = os.listdir(stag_imgs_dir)
stag_imgs = [f for f in stag_imgs if f.endswith('.tif')]

# remove extension from stag_imgs
stag_imgs = [f.split('.')[0] for f in stag_imgs]




sub_text = get_sub_text(stag_imgs)



# generate mrml from template
with open("./templates/rigid/rigid.mrml", "rt") as fin:
    with open(f'{snakemake.output.mrml}', "wt") as fout:
        for line in fin:
            # line = line.replace('fname_nis', f'{snakemake.wildcards.fname}')
            # line = line.replace('fname_stag', f'stag_{snakemake.wildcards.fname}')
            line = line.replace('fname_nis', f'{snakemake.config["dataname"]}_nissl')
            line = line.replace('fname_stag', f'bead_plot')
            line = line.replace('<Hidden></Hidden>', f'{sub_text}')
            fout.write(line)

# os.system(f'touch {snakemake.output[0]}')
# os.system(f'touch {snakemake.output[1]}')
# os.system(f'touch {snakemake.output[2]}')
# os.system(f'touch {snakemake.output[3]}')

os.system(f'cp {snakemake.input.from_fids} {snakemake.output.from_fids}')
os.system(f'cp {snakemake.input.to_fids} {snakemake.output.to_fids}')
os.system(f'cp {snakemake.input.tfm1} {snakemake.output.tfm1}')
