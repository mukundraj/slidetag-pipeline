
import os
print ('inpy', {snakemake.output.mrml})


with open("./templates/rigid.mrml", "rt") as fin:

# generate mrml from template
with open("./templates/rigid/rigid.mrml", "rt") as fin:
    with open(f'{snakemake.output.mrml}', "wt") as fout:
        for line in fin:
            # line = line.replace('fname_nis', f'{snakemake.wildcards.fname}')
            # line = line.replace('fname_stag', f'stag_{snakemake.wildcards.fname}')
            line = line.replace('fname_nis', f'{snakemake.config["dataname"]}_nissl')
            line = line.replace('fname_stag', f'bead_plot')
            fout.write(line)

# os.system(f'touch {snakemake.output[0]}')
# os.system(f'touch {snakemake.output[1]}')
# os.system(f'touch {snakemake.output[2]}')
# os.system(f'touch {snakemake.output[3]}')

os.system(f'cp {snakemake.input.from_fids} {snakemake.output.from_fids}')
os.system(f'cp {snakemake.input.to_fids} {snakemake.output.to_fids}')
os.system(f'cp {snakemake.input.tfm1} {snakemake.output.tfm1}')
