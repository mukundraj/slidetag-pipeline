## Sample data to try pipeline

Download [sample data](https://drive.google.com/drive/folders/1Fzp5OZB1giX962EspsRdKo92qgRRl0o8?usp=drive_link)

## Steps for slide tag image alignment

### Part 0: Setup chrome remote desktop

- Follow steps [here](https://support.google.com/chrome/answer/1649523?hl=en&co=GENIE.Platform%3DDesktop) to connect to instance via remote desktop. This step assumes GCP instance is already set up. If not, see appendix to set up GCP instance.

### Part 1: Data initialization (run once for each dataset)

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

## Appendix

### A1: Instance setup steps (only needed for instance setup, ignore otherwise)

- Run following script on local machine to create instance machine

```
bash src/setup_instance.sh
```

- Run setup script on instance machine

```
bash src/startup.sh
```

### A2: Use following command arguments for local testing of pipeline on osx

- note that different config template file is passed as argument than in case of linux instance

**Part 1: Data initialization command\_**

```
  bash src/workflow/prep_init.sh DATASET_NAME templates/config_osx.yaml
```
