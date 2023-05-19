## Steps for slide tag image alignment

### Sample data to test pipeline

### Part 0: Setup (only first time)

- Start instance

```
  git clone https://github.com/mukundraj/slidetag-pipeline.git
  cd slidetag-pipeline && mkdir build
```

### Part 1: Initialization (run once for each dataset)

- copy nissls to instance

```
  gcloud compute scp nissls.zip <INSTANCE_NAME>:~/
```

- copy slide tab images to instance

```
  gcloud compute scp stags.zip <INSTANCE_NAME>:~/
```

```
  conda activate snakemake
  bash src/workflow/prep_init.sh <DATASET_NAME> templates/config.yaml
```

### Part 2: Rigid alignment

- perform alignment in slicer
- run command 2:

```
bash src/workflow/prep_rigid.sh
```

### Part 3: Nonlinear alignment

- perform alignment in slicer
- run command 3

```
bash src/workflow/render_plots.sh
```
