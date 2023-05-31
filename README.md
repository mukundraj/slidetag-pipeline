## Sample data to try pipeline

Download [sample data](https://drive.google.com/drive/folders/1Fzp5OZB1giX962EspsRdKo92qgRRl0o8?usp=drive_link)

## Steps for slide tag image alignment

### Part 0: Setup

- **_Connect to remote desktop_** Follow steps [here](https://support.google.com/chrome/answer/1649523?hl=en&co=GENIE.Platform%3DDesktop) to connect to instance via remote desktop. This step assumes GCP instance is already set up. If not, see appendix to set up GCP instance.

- **_Activate snakemake conda environment_** needed for alignment workflow scripts. This command should be run before running any of the scripts in the following parts of pipeline.

```
  conda activate snakemake
```

### Part 1: Data initialization (run once for each dataset)

- copy dataset to GCP instance as zip file containing 3 items - nissl image, seurat object in qs format, bead coordinates in csv file.

```
  gcloud compute scp DATASET_NAME.zip st-alignment:~/
```

- run following command in terminal to extract and prepare data to start alignment process.

```
  bash src/workflow/prep_init.sh DATASET_NAME
```

### Part 2: Rigid alignment

- identify 3 pairs of fiducial points in slicer
- run following script to perform rigid alignment and prepare for nonlinear alignment

```
bash src/workflow/prep_warp.sh
```

### Part 3: Nonlinear alignment

- identify as many pairs of fiducials as needed in slicer
- run following script to perform nonlinear alignment and prepare overlay plots

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
