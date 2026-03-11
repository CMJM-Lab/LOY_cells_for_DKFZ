#!/bin/bash

pyGenomeTracks \
    --tracks ./gene_tracks.ini \
    --region chr5:122061195-122080259 \ #LOX: chr5:122061195-122080259; THY1: chr11:119410000-119425000
    --outFileName ./gene_plot.pdf \
    --fontSize 10 \
    --dpi 300 \
    --width 40
