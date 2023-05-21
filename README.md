## Sample data to try pipeline

## Steps for slide tag image alignment

### Part 0: Setup chrome remote desktop

### Part 1: Date initialization (run once for each dataset)

- copy nissls to instance

```
  gcloud compute scp nissls.zip st-alignment:~/
```

- copy slide tab images to instance

```
  gcloud compute scp stags.zip st-alignment:~/
```

- activate conda environment needed for alignment workflow scripts

```
  conda activate snakemake
```

- run following command in terminal to format images

```
  bash src/workflow/prep_init.sh DATASET_NAME templates/config.yaml
```

### Part 2: Rigid alignment

- run following command to prepare slicer files for rigid alignment:

```
bash src/workflow/prep_rigid.sh
```

- identify 3 pairs of fiducial points in slicer

### Part 3: Nonlinear alignment

- run following command to prepare slicer files for nonlinear alignment:

```
bash src/workflow/prep_warp.sh
```

- identify as many pairs of fiducials as needed in slicer

### Part 4: Create plots

```
bash src/workflow/render_plots.sh
```

## Appendix: Instance setup steps (only needed for instance setup, ignore otherwise)

- Create instance
- Run setup script
