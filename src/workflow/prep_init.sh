#!/bin/bash

# bash prep.sh <dataset_name> <optional:config_template>
# bash prep.sh a5 ./config/config.yaml

if [ "$#" -le 0 ]; then
    echo "Too few number of parameters. At least dataset name and path to config file template are required." >& 2
    exit 2
fi

# # prep the config file - populate username and dataset name
# python ./src/workflow/prep_init.py $1 $2 $(whoami)


# cp snakefile
# cp ./src/workflow/Snakefile ./build/Snakefile
pipeline_root=$(dirname "$0")/../..
echo $pipeline_root

data_root=$PWD

mkdir -p build # make build dir if directory doesn't exist
# the relative path of the source file must be specified relative to the link file directory, not to the current directory
ln -nsf $pipeline_root/src/workflow/Snakefile build/Snakefile
# ln -s ../src/workflow/format_imgs.py build/format_imgs.py

# run the snakemake initial rule for unzipping images
# snakemake --cores 1 injest_imgs

# build warping code in c++
cd $pipeline_root/build
if [[ $OSTYPE == "linux-gnu" ]]; then
    python $pipeline_root/src/workflow/prep_init.py $1 $pipeline_root/templates/config.yaml $(whoami) $data_root $pipeline_root #prep config file
    cmake -DITK_DIR=/usr/src/InsightToolkit-5.3.0/build ..
else
    python $pipeline_root/src/workflow/prep_init.py $1 $pipeline_root/templates/config_osx.yaml $(whoami) $data_root $pipeline_root #prep config file
    cmake -DITK_DIR=/Users/mraj/Desktop/work/pkgsources/ITK/build ..
fi
make 

# return to data directory
cd -

ln -nsf $pipeline_root/build/cmapper2 build/cmapper2

echo 'inbash' $1 $2 

snakemake --use-conda --cores 1 --snakefile ./build/Snakefile prep_rigid
