
import os
print ('inpy', {snakemake.output.mrml})


with open("./templates/rigid.mrml", "rt") as fin:
    with open(f'{snakemake.output.mrml}', "wt") as fout:
        for line in fin:
            line = line.replace('fname_nis', f'{snakemake.wildcards.fname}')
            line = line.replace('fname_stag', f'stag_{snakemake.wildcards.fname}')
            fout.write(line)

# os.system(f'touch {snakemake.output[0]}')
# os.system(f'touch {snakemake.output[1]}')
# os.system(f'touch {snakemake.output[2]}')
# os.system(f'touch {snakemake.output[3]}')

os.system(f'cp templates/F.mrk.json {snakemake.output.from_fids}')
os.system(f'cp templates/T.mrk.json {snakemake.output.to_fids}')
os.system(f'cp templates/Tfm1.txt {snakemake.output.tfm1}')
