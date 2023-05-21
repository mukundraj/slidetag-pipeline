#!/bin/bash

# bash prep.sh <dataset_name> <optional:config_template>
# bash prep.sh a5 ./config/config.yaml

if [ "$#" -le 1 ]; then
    echo "Too few number of parameters. At least dataset name and path to config file template are required." >& 2
    exit 2
fi

# prep the config file - populate username and dataset name
python ./src/workflow/prep_init.py $1 $2 $(whoami)


# cp snakefile
# cp ./src/workflow/Snakefile ./build/Snakefile

# the relative path of the source file must be specified relative to the link file directory, not to the current directory
ln -s ../src/workflow/Snakefile build/Snakefile
# ln -s ../src/workflow/format_imgs.py build/format_imgs.py

# run the snakemake initial rule for unzipping images
# snakemake --cores 1 injest_imgs

# build warping code in c++
cd build
if [[ $OSTYPE == "linux-gnu" ]]; then
    cmake -DITK_DIR=/usr/src/InsightToolkit-5.3.0/build ..
else
    cmake -DITK_DIR=/Users/mraj/Desktop/work/pkgsources/ITK/build ..
fi
make 
cd ..

echo 'inbash' $1 $2 

snakemake --cores 1 --snakefile build/Snakefile format_imgs
