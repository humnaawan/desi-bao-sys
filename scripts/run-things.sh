#!/bin/bash

create_per_petal_spectra=1
night=20191206

get_qso_list=1
group_and_coadd_spectra=1
nside=16

get_zfile=1
get_zcatalog=1
get_drqcatalog=1
run_picca=1
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
# ---------------------------------------------------------
if [ $get_zcatalog == 1 ];
then
     printf '\n## running desi_zcatalog ...\n'
     # https://github.com/desihub/desispec/blob/main/bin/desi_zcatalog
     desi_zcatalog --indir=${datapath} \
                    --outfile=${datapath}'zcat.fits' \
                    --minimal \
                    --prefix='ztrue'\
                    > ${datapath}'out_desi_zcatalog.log'
fi
# ---------------------------------------------------------
if [ $get_drqcatalog == 1 ];
then
     printf '\n## running get_drqcatalog.py ...\n'
     module load python
     conda activate picca_pip
     python /global/homes/a/awan/desi/desi-bao-sys/scripts/get-drqcatalog.py \
                    --zcat-path=${datapath}'zcat.fits' \
                    --outdir=${datapath} \
                    > ${datapath}'out_drqcatalog.log'
fi
# ---------------------------------------------------------
if [ $run_picca == 1 ];
then
     printf '\n## running run_picca.py ...\n'
     module load python
     conda activate picca_pip
     python /global/homes/a/awan/desi/desi-bao-sys/scripts/run-picca.py \
                         --outdir=${basepath}'picca-output/' \
                         --zcat-drq-path=${datapath}'zcat_drq.fits' \
                         --spectra-path=${datapath}'spectra-16/' \
                         --get-deltas \
                         --calc-corrs \
                         --get-dmas \
                         --fit-bao \
                         --inis-path=${basepath}/'picca-inis/' \
                         > ${basepath}'picca-output/out_runpicca.log'
fi
