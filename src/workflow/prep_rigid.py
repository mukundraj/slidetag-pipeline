
# with open("Stud.txt", "rt") as fin:
#     with open("out.txt", "wt") as fout:
#         for line in fin:
#             fout.write(line.replace('A', 'Orange'))
import os
print ('inpy', {snakemake.output[0]})


os.system(f'touch {snakemake.output[0]}')
os.system(f'touch {snakemake.output[1]}')
os.system(f'touch {snakemake.output[2]}')
os.system(f'touch {snakemake.output[3]}')
