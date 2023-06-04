#!/bin/bash

snakemake --use-conda --cores 1 --snakefile build/Snakefile render_plots
