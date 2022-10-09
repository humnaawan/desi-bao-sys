#!/bin/bash

create_per_petal_spectra=1
night=20191206

get_qso_list=1
group_and_coadd_spectra=1
nside=16

get_zfile=1
# --------------------------------------------------------------------------
source /global/cfs/cdirs/desi/software/desi_environment.sh master

basepath='/global/cfs/cdirs/desi/users/awan/'
datapath=${basepath}'exposures/'

# ---------------------------------------------------------
if [ $create_per_petal_spectra == 1 ];
then
    for expid in ${datapath}/${night}/*
    do
         printf '\n## running create-per-petal-spectra.py ...\n'
         python /global/homes/a/awan/desi/desi-bao-sys/scripts/create-per-petal-spectra.py \
                        --data-path=${datapath} \
                        --night=${night} \
                        --expid=$(basename ${expid})
    done
fi
# ---------------------------------------------------------
if [ $get_qso_list == 1 ];
then
     printf '\n## running get-qso-exposures-list.py ...\n'
     python /global/homes/a/awan/desi/desi-bao-sys/scripts/get-qso-exposures-list.py \
                        --data-path=${datapath} \
                        --nside=${nside}
fi
# ---------------------------------------------------------
if [ $group_and_coadd_spectra == 1 ];
then
     printf '\n## running group-and-coadd-spectra.py ...\n'
     python /global/homes/a/awan/desi/desi-bao-sys/scripts/group-and-coadd-spectra.py \
                        --exposures-list-path=${datapath}'qso-exposures-list_nside'${nside}'.csv' \
                        --exposures-path=${basepath} \
                        --nside=${nside} \
                        --outdir=${datapath}
fi
# ---------------------------------------------------------
if [ $get_zfile == 1 ];
then
     printf '\n## running save-zfile-per-pixel.py ...\n'
     python /global/homes/a/awan/desi/desi-bao-sys/scripts/get-zfile-per-pixel.py \
                        --simspec-path=${datapath} \
                        --coadds-path=${datapath} \
                        --nside=${nside}
fi
