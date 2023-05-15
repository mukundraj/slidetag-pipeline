
import os



# copy stag image after applying tfms and lift bbox
# os.system(f'touch {snakemake.output.nissl}')
os.system(f'cp {snakemake.input.nissl} {snakemake.output.nissl}')

# copy nissl image after cropping bbox after transform stag coords
# os.system(f'touch {snakemake.output.stag}')
os.system(f'cp {snakemake.input.stag} {snakemake.output.stag}')

# os.system(f'touch {snakemake.output.mrml}')
# os.system(f'touch {snakemake.output.from_fids}')
# os.system(f'touch {snakemake.output.to_fids}')

# copy template files
with open("./templates/warp/warp.mrml", "rt") as fin:
    with open(f'{snakemake.output.mrml}', "wt") as fout:
        for line in fin:
            line = line.replace('fname_nis', f'{snakemake.wildcards.fname}')
            line = line.replace('fname_stag', f'stag_{snakemake.wildcards.fname}')
            fout.write(line)

os.system(f'cp templates/warp/F.mrk.json {snakemake.output.from_fids}')
os.system(f'cp templates/warp/T.mrk.json {snakemake.output.to_fids}')
# os.system(f'touch {snakemake.output.tfm2}') # don't add this to output, empty file takes slicer to undesired state
