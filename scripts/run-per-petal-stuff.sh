#!/bin/bash

source /global/cfs/cdirs/desi/software/desi_environment.sh master

datapath='/global/cfs/cdirs/desi/users/awan/'
night=20191206

for expid in ${datapath}/${night}/*
do
    python /global/homes/a/awan/desi/desi-bao-sys/scripts/create-per-petal-spectra.py \
                        --data-path=${datapath} \
                        --night=${night} \
                        --expid=$(basename ${expid})

    python /global/homes/a/awan/desi/desi-bao-sys/scripts/group-and-coadd-spectra.py \
                        --data-path=${datapath} \
                        --night=${night} \
                        --expid=$(basename ${expid})
done

