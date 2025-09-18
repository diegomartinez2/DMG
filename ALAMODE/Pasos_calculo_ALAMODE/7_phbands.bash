#!/bin/bash
echo "first build the 'anphon_phband.in' with 'Construir_alm_anphon_file.py'"
anphon anphon_phband.in
echo "change system.bands for the *.bands output from anphon"
plotband.py system.bands
