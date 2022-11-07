# ------------------------------------------------------------------------
# this script runs picca; different intermediate steps are available
# as options.
# following https://desi.lbl.gov/trac/wiki/LymanAlphaWG/how_to_run_picca
# ------------------------------------------------------------------------
import sys
import os
import picca
import time
import datetime
import logging
import subprocess

from optparse import OptionParser

parser = OptionParser()
parser.add_option('--outdir',
                  dest='outdir', type='str',
                  help='path to directory where things should be saved.')
parser.add_option('--zcat-drq-path',
                  dest='path_drq', type='str',
                  help='path to zcat-drq file.')
parser.add_option('--spectra-path',
                  dest='path_spec16', type='str',
                  help='path to spectra-<nside> folder.')
parser.add_option('--get-deltas',
                  action='store_true', dest='get_deltas', default=False,
                  help='run picca_deltas and plot delta-attributes.')
parser.add_option('--calc-corrs',
                  action='store_true', dest='calc_corrs', default=False,
                  help='run picca_cf.')
parser.add_option('--plot-corrs',
                  action='store_true', dest='plot_corrs', default=False,
                  help='plot correlations; if calc-corrs, this would automatically be set to True.')
parser.add_option('--get-dmas',
                  action='store_true', dest='get_dmas', default=False,
                  help='run picca_dmas (distortion matrices) and picca_export.')
parser.add_option('--fit-bao',
                  action='store_true', dest='fit_bao', default=False,
                  help='picca_fitter2; this will require inis-path set.')
parser.add_option('--plot-baofit',
                  action='store_true', dest='plot_baofit', default=False,
                  help='plot bao-fit; if fit-bao, this would automatically be set to True.')
parser.add_option('--inis-path',
                  dest='inis_path', type='str', default=None,
                  help='run to the ini files.')

(options, args) = parser.parse_args()
outdir = options.outdir
path_drq = options.path_drq
path_spec16 = options.path_spec16
get_deltas = options.get_deltas
calc_corrs = options.calc_corrs
plot_corrs = options.plot_corrs
get_dmas = options.get_dmas
fit_bao = options.fit_bao
plot_baofit = options.plot_baofit
inis_path = options.inis_path
# ------------------------------------------------------------------------
start0 = time.time()
print('\n## inputs: %s' % options)
# set up the logger
temp = f'save-run-picca_{datetime.datetime.now()}'.replace(' ', '_').split('.')[0]
logging.basicConfig(filename=f'{outdir}/log_{temp}.log',
                    level=logging.DEBUG, filemode='w',
                   )
logging.info(f'running {sys.argv[0]}')
logging.info(f'\n## inputs: {options}\n')

logging.getLogger('matplotlib.font_manager').disabled = True
logging.getLogger("imported_module").setLevel(logging.WARNING)
#logging.info(f'\n ## inputs: {options}\n')

if calc_corrs and not plot_corrs:
    logging.info(f'## setting plot_corrs = True since calc_corrs = True')
    plot_corrs = True
if fit_bao and not plot_baofit:
    logging.info(f'## setting plot_baofit = True since fit_bao = True')
    plot_baofit = True
if False:
    outdir = '/global/cfs/cdirs/desi/users/awan/picca-output/'
    path_drq = '/global/cfs/cdirs/desi/users/awan/exposures/zcat_drq.fits'
    path_spec16 = '/global/cfs/cdirs/desi/users/awan/exposures/spectra-16/'
    get_deltas = False
    calc_corrs = False
    plot_corrs = True
    get_dmas = False
    fit_bao = False
    plot_baofit = False
    inis_path = f'{outdir}/../picca-inis/'  # assume one level up from outdir

# create the plots directory
plotsdir = f'{outdir}/plots'
if not os.path.exists(plotsdir):
    logging.info(f'## creating {plotsdir}')
    os.makedirs(plotsdir)

# lets set up the delta folder
path_deltas = f'{outdir}/delta_lya'
if not os.path.exists(path_deltas):
    logging.info(f'## creating {path_deltas}')
    os.makedirs(path_deltas)
path_deltas = None   # i dont think this path will be needed after this
# now the subfolders
path_deltas_delta = f'{outdir}/delta_lya/delta'
if not os.path.exists(path_deltas_delta):
    logging.info(f'## creating {path_deltas_delta}')
    os.makedirs(path_deltas_delta)
path_deltas_log = f'{outdir}/delta_lya/log'
if not os.path.exists(path_deltas_log):
    logging.info(f'## creating {path_deltas_log}')
    os.makedirs(path_deltas_log)

data_tag = 'lyalya_lyalya'

if get_deltas:
    start1 = time.time()
    logging.info(f'\n-- running picca continuum fitting ...')
    print(f'\n-- running picca continuum fitting ...')

    # picca continuum fitting
    bash_command = 'picca_deltas.py ' + \
                    f'--drq {path_drq} ' + \
                    f'--in-dir {path_spec16} ' + \
                    f'--out-dir {path_deltas_delta} ' + \
                    f'--mode desi ' + \
                    f'--log {path_deltas_log}/input.log ' + \
                    f'--iter-out-prefix {path_deltas_log}/delta_attributes'
                    #f'--nproc 32 '
    logging.info(f'-- running bash_command =\n{bash_command}\n')
    process = subprocess.Popen(bash_command.split(), stdout=subprocess.PIPE)
    output, error = process.communicate()
    if error is not None:
        print(error)
    logging.info(f'-- command output:\n{output}\n')
    logging.info(f'-- command error:\n{error}\n')

    # lets plot some of the data saved
    from helpers_plots import plot_delta_attributes
    fname = plot_delta_attributes(data_path=path_deltas_log,
                                  outdir=plotsdir)
    logging.info(f'-- saved {fname}')
    logging.info(f'## time taken for get_deltas: {(time.time() - start1)/60: .2f} min')

if calc_corrs:
    start1 = time.time()
    logging.info(f'\n-- calculating correlations using picca ...')
    print(f'\n-- calculating correlations using picca ...')
    path_cfs = f'{outdir}/correlations'
    if not os.path.exists(path_cfs):
        logging.info(f'## creating {path_cfs}')
        os.makedirs(path_cfs)

    # run picca_cf
    bash_command = 'picca_cf.py ' + \
                    f'--in-dir {path_deltas_delta}/ ' + \
                    f'--out {path_cfs}/cf_{data_tag}.fits.gz '
                    #f'--nproc 32 '
    logging.info(f'running bash_command =\n{bash_command}\n')
    process = subprocess.Popen(bash_command.split(), stdout=subprocess.PIPE)
    output, error = process.communicate()
    if error is not None:
        print(error)
    logging.info(f'-- command output:\n{output}\n')
    logging.info(f'-- command error:\n{error}\n')
    logging.info(f'## time taken for calc_corrs: {(time.time() - start1)/60: .2f} min')

if plot_corrs:
    path_cfs = f'{outdir}/correlations'
    if not os.path.exists(path_cfs):
        message = f'{path_cfs} not found - this means that cfs dont exist already!?!'
        logging.critical(f'## {message}')
        raise ValueError(message)

    from helpers_plots import plots_cfs
    # lets plot some of the data saved
    fname = plots_cfs(data_path=path_cfs,
                      outdir=plotsdir, data_tag=data_tag)
    logging.info(f'-- saved {fname}')

if get_dmas:
    start1 = time.time()
    logging.info(f'\n-- dealing with distortion matrices ...')
    print(f'\n-- dealing with distortion matrices ...')
    path_cfs = f'{outdir}/correlations'
    if not os.path.exists(path_cfs):
        message = f'{path_cfs} not found - this means that cfs dont exist already!?!'
        logging.critical(f'## {message}')
        raise ValueError(message)

    # run picca_dmat
    bash_command = 'picca_dmat.py ' + \
                    f'--in-dir {path_deltas_delta}/ ' + \
                    f'--out {path_cfs}/dmat_{data_tag}.fits.gz ' + \
                    '--rej 0.99 '
                    #f'--nproc 32 '
    logging.info(f'-- running bash_command =\n{bash_command}\n')
    process = subprocess.Popen(bash_command.split(), stdout=subprocess.PIPE)
    output, error = process.communicate()
    if error is not None:
        print(error)
    logging.info(f'-- command output:\n{output}\n')
    logging.info(f'-- command error:\n{error}\n')

    # run picca_export
    bash_command = 'picca_export.py ' + \
                    f'--data {path_cfs}/cf_{data_tag}.fits.gz ' + \
                    f'--dmat {path_cfs}/dmat_{data_tag}.fits.gz ' + \
                    f'--out {path_cfs}/cf_exp_{data_tag}.fits.gz '
                    #f'--nproc 32 '
    logging.info(f'-- running bash_command =\n{bash_command}\n')
    process = subprocess.Popen(bash_command.split(), stdout=subprocess.PIPE)
    output, error = process.communicate()
    if error is not None:
        print(error)
    logging.info(f'-- command output:\n{output}\n')
    logging.info(f'-- command error:\n{error}\n')

    if False:
        # not sure if need to run this for these mocks
        # run picca_dmat
        bash_command = 'picca_metal_dmat.py ' + \
                        f'--in-dir {path_deltas_delta}/ ' + \
                        f'--out {path_cfs}/metal_dmat_{data_tag}.fits.gz ' + \
                        '--rej 0.99 ' + \
                        '--abs-igm SiII\(1260\) SiIII\(1207\) SiII\(1193\) SiII\(1190\)'
                        #f'--nproc 32 '
        logging.info(f'-- running bash_command =\n{bash_command}\n')
        process = subprocess.Popen(bash_command.split(), stdout=subprocess.PIPE)
        output, error = process.communicate()
        if error is not None:
            print(error)
        logging.info(f'-- command output:\n{output}\n')
        logging.info(f'-- command error:\n{error}\n')

    logging.info(f'## time taken for distortion matrices: {(time.time() - start1)/60: .2f} min')

if fit_bao:
    import shutil
    start1 = time.time()
    logging.info(f'\n-- bao fitting using picca ...')
    print(f'\n-- bao fitting using picca ...')
    path_fits = f'{outdir}/baofits'
    if not os.path.exists(path_fits):
        logging.info(f'## creating {path_fits}')
        os.makedirs(path_fits)

    for file in [f'config_{data_tag}.ini', f'chi2_{data_tag}.ini']:
        filepath = f'{inis_path}/{file}'
        if not os.path.exists(filepath):
            message = f'{filepath} not found'
            logging.critical(f'## {message}')
            raise ValueError(message)
        else:
            shutil.copy(filepath, f'{path_fits}/{file}')
            logging.info(f' -- copied {file} from {inis_path} to {path_fits}')

    # run picca_fitter
    bash_command = f'picca_fitter2.py {path_fits}/chi2_{data_tag}.ini'
    logging.info(f'running bash_command =\n{bash_command}\n')
    process = subprocess.Popen(bash_command.split(), stdout=subprocess.PIPE)
    output, error = process.communicate()
    if error is not None:
        print(error)
    logging.info(f'-- command output:\n{output}\n')
    logging.info(f'-- command error:\n{error}\n')
    logging.info(f'## time taken for bao fitting: {(time.time() - start1)/60: .2f} min')

if plot_baofit:
    # plot
    path_cfs = f'{outdir}/correlations'
    if not os.path.exists(path_cfs):
        message = f'{path_cfs} not found - this means that cfs dont exist already!?!'
        logging.critical(f'## {message}')
        raise ValueError(message)

    path_fits = f'{outdir}/baofits'
    if not os.path.exists(path_fits):
        message = f'{path_fits} not found!?!'
        logging.critical(f'## {message}')
        raise ValueError(message)

    from helpers_plots import plot_baofit
    fname = plot_baofit(path_cf=f'{path_cfs}/cf_exp_{data_tag}.fits.gz',
                        outdir=plotsdir, data_tag=data_tag,
                        path_baofit=f'{path_fits}/result_{data_tag}.h5',
                        power=2,
                        mus=[1., 0.95, 0.8, 0.5, 0],
                        fit_rmax=180.,
                        absMus=False)

    logging.info(f'-- saved {fname}')


logging.info('\n## all done')
logging.info(f'## time taken: {(time.time() - start0)/60: .2f} min')
