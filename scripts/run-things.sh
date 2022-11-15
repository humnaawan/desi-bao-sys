#!/bin/bash

create_per_petal_spectra=1

get_qso_list=1
group_and_coadd_spectra=1
nside=16

get_zfile=1
get_zcatalog=1
get_drqcatalog=1
run_picca=1
diagnostic_plots=1
# --------------------------------------------------------------------------
source /global/cfs/cdirs/desi/software/desi_environment.sh master

path_base='/global/cfs/cdirs/desi/users/awan/'
path_sims=${path_base}'2019-sims/'
path_exposures=${path_sims}'exposures/'

# ---------------------------------------------------------
if [ $create_per_petal_spectra == 1 ];
then
     for night in ${path_exposures}/*
     do
          night=$(basename ${night})
          printf '\n## looking for all exposures in night=%s ..' ${night}
          for expid in ${path_exposures}/${night}/*
          do
               printf '\n## looking at night=%s and exposure=%s ..\n' ${night} $(basename ${expid})
               printf '## running create-per-petal-spectra.py ...\n'
               python /global/homes/a/awan/desi/desi-bao-sys/scripts/create-per-petal-spectra.py \
                             --data-path=${path_exposures} \
                             --night=${night} \
                             --expid=$(basename ${expid})
          done
     done
fi
# ---------------------------------------------------------
if [ $get_qso_list == 1 ];
then
     printf '\n## running get-qso-exposures-list.py ...\n'
     python /global/homes/a/awan/desi/desi-bao-sys/scripts/get-qso-exposures-list.py \
                        --data-path=${path_exposures} \
                        --outdir=${path_sims} \
                        --nside=${nside}
fi
# ---------------------------------------------------------
if [ $group_and_coadd_spectra == 1 ];
then
     printf '\n## running group-and-coadd-spectra.py ...\n'
     python /global/homes/a/awan/desi/desi-bao-sys/scripts/group-and-coadd-spectra.py \
                        --exposures-list-path=${path_sims}'qso-exposures-list_nside'${nside}'.csv' \
                        --exposures-path=${path_sims} \
                        --nside=${nside} \
                        --outdir=${path_sims}
fi
# ---------------------------------------------------------
if [ $get_zfile == 1 ];
then
     printf '\n## running save-zfile-per-pixel.py ...\n'
     python /global/homes/a/awan/desi/desi-bao-sys/scripts/get-zfile-per-pixel.py \
                        --simspec-path=${path_exposures} \
                        --coadds-path=${path_sims} \
                        --nside=${nside}
fi
# ---------------------------------------------------------
if [ $get_zcatalog == 1 ];
then
     printf '\n## running desi_zcatalog ...\n'
     # https://github.com/desihub/desispec/blob/main/bin/desi_zcatalog
     desi_zcatalog --indir=${path_sims} \
                    --outfile=${path_sims}'zcat.fits' \
                    --minimal \
                    --prefix='ztrue'\
                    > ${path_sims}'out_desi_zcatalog.log'
fi
# ---------------------------------------------------------
if [ $get_drqcatalog == 1 ];
then
     printf '\n## running get_drqcatalog.py ...\n'
     module load python
     conda activate picca_pip
     python /global/homes/a/awan/desi/desi-bao-sys/scripts/get-drqcatalog.py \
                    --zcat-path=${path_sims}'zcat.fits' \
                    --outdir=${path_sims} \
                    > ${path_sims}'out_drqcatalog.log'
fi
# ---------------------------------------------------------
if [ $run_picca == 1 ];
then
     printf '\n## running run_picca.py ...\n'
     module load python
     conda activate picca_pip
     python /global/homes/a/awan/desi/desi-bao-sys/scripts/run-picca.py \
                         --outdir=${path_sims} \
                         --zcat-drq-path=${path_sims}'zcat_drq.fits' \
                         --spectra-path=${path_sims}'spectra-16/' \
                         --get-deltas \
                         --calc-corrs \
                         --get-dmas \
                         --fit-bao \
                         --inis-path=${path_base}/'picca-inis/' \
                         > ${path_sims}'out_runpicca.log'
fi
# ---------------------------------------------------------
if [ $diagnostic_plots == 1 ];
then
     printf '\n## running create-diagnostic-plots.py ...\n'
     module load python
     conda activate picca_pip
     python /global/homes/a/awan/desi/desi-bao-sys/scripts/create-diagnostic-plots.py \
                         --outdir=${path_sims} \
                         --exposures-path=${path_exposures} \
                         --qso-cat-path=${path_sims}'qso-exposures-list_nside'${nside}'.csv' \
                         --zcat-path=${path_sims}'zcat.fits' \
                         --spectra-path=${path_sims}'spectra-16/' \
                         --nside=16
fi
