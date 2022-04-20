# ------------------------------------------------------------------------
# this script takes in the fibermap and simspec simulation files and
# produces the spectra. if using the debug option, will produce the
# per-petal fibermap and simspec files.
# most of the code is based on quickspectra/quickquasars
# ------------------------------------------------------------------------
from astropy.io import fits
from astropy.table import Table
import astropy.units as u
import numpy as np
import matplotlib.pyplot as plt

import desisim.simexp
import desimodel.io
import desispec.io
from desispec.spectra import Spectra
from desispec.resolution import Resolution
from desispec.frame import Frame
from desispec.specscore import compute_and_append_frame_scores

import logging # would use desi logger but dont know where it saves stuff
import sys
import time
import datetime

from optparse import OptionParser
parser = OptionParser()
parser.add_option('--data-path',
                  dest='data_path', type='str',
                  help='path to directory with all the data \
                        (i.e., where folder name = night)')
parser.add_option('--night',
                  dest='night', type='str',
                  help='night to consider')
parser.add_option('--expid',
                  dest='expid', type='str',
                  help='exposure ID to consider')
parser.add_option('--debug',
                  action='store_true', dest='debug', default=False,
                  help='Run things in debug mode \
                        => some intermediate plots are created and \
                        per-petal fibermap and simsec files are saved.')
parser.add_option('--exptime', default=None,
                  dest='exptime', type='float',
                  help='exposure time in seconds; \
                        will override the value in simspec file; \
                        valid only in debug mode.')
parser.add_option('--skyerr', default=0,
                  dest='skyerr', type='float',
                  help='skyerr value; default = 0')
(options, args) = parser.parse_args()
print('\n## inputs: %s' % options)

data_path = options.data_path
night = options.night
expid = options.expid
debug = options.debug
exptime = options.exptime
skyerr = options.skyerr
# ------------------------------------------------------------------------
start0 = time.time()

# one night's data
#night = '20191206'
#expid = '00000344'
#path = f'/global/cscratch1/sd/awan/desi/{night}/{expid}/'
path = f'{data_path}/{night}/{expid}/'

# set up the logger
temp = f'{datetime.datetime.now()}'.replace(' ', '_').split('.')[0]
if debug: temp = f'debug_{temp}'
logging.basicConfig(filename=f'{path}/log_{temp}.log',
                    level=logging.DEBUG, filemode='w',
                    #format='%(asctime)s %(message)s', datefmt='%m/%d/%Y %I:%M:%S %p'
                   )
logging.info(f'running {sys.argv[0]}')
logging.info(f'\n## inputs: {options}\n')

# ------------------------------------------------------------------------
# fibermap data
fibermap_fname = f'{path}/fibermap-{expid}.fits'
fibermap = Table.read(fibermap_fname, format='fits')
hdul_fibermap = fits.open(fibermap_fname)

for key in ['FIBER', 'TARGETID', 'DEVICE_LOC']:
    uniq, cts = np.unique(fibermap[key], return_counts=True)
    logging.info(f'{len(uniq)} unique {key} with {np.unique(cts)} unique counts across')

# ------------------------------------------------------------------------
# simspec data
simspec_fname = f'{path}/simspec-{expid}.fits'
simspec = Table.read(simspec_fname, format='fits', hdu='TRUTH')
simspec_hdul = fits.open(simspec_fname)
obsconditions = Table.read(simspec_fname, format='fits', hdu='OBSCONDITIONS')
obsconditions_dict = [dict(zip(obsconditions.colnames, row)) for row in obsconditions][0]
if debug and exptime is not None:
    logging.info(f'# updating exptime to be {exptime}s.')
    obsconditions_dict['EXPTIME'] = exptime

# update the exptime in the fibermap
fibermap['EXPTIME'][:] = obsconditions_dict['EXPTIME']
#print(f'obsconditions_dict: {obsconditions_dict}')
# ------------------------------------------------------------------------------------
def save_petal_spectra(simulator, flux_this_petal, spectra_fibermap, spectra_filename,
                       seed=10, skyerr=0, nfiber_to_plot=400):
    """
    spectra_fibermap is essentially fibermap with added EXPID and NIGHT columns
    (ref: https://github.com/desihub/desisim/blob/8308fc44cdc86aea14155b0db5c6f529eeea8423/py/desisim/scripts/quickspectra.py#L103)
    """
    logging.info(f'## running save_spectra_for_petal ...')
    scale = 1e17          # not exactly sure where this is coming from
    specdata = None       # initialize variable for the Spectra object
    nspec = flux_this_petal.shape[0] # number of spectra
    # initialize random state
    random_state = np.random.RandomState(seed)
    # skyscale
    skyscale = skyerr * random_state.normal(size=simulator.num_fibers)

    # set up the Resolution for each camera
    resolution = {}
    for camera in simulator.instrument.cameras:
        R = Resolution(camera.get_output_resolution_matrix())
        resolution[camera.name] = np.tile(R.to_fits_array(), [nspec, 1, 1]) # need to check if this is right since quickgen seems to have a different way to handle this

    # create a Spectra object for each camera
    for table in simulator.camera_output :
        wave = table['wavelength'].astype(float)
        # base flux: observed-flux + noise-electrons * flux-calibration
        flux = (table['observed_flux'] + \
                table['random_noise_electrons'] * table['flux_calibration']
               ).T.astype(float)
        # add sky if available
        if np.any(skyscale):
            flux += ((table['num_sky_electrons'] * skyscale) * \
                     table['flux_calibration']
                    ).T.astype(float)
        # flux finverse variance for Spectra
        ivar = table['flux_inverse_variance'].T.astype(float)
        band  = table.meta['name'].strip()[0]

        # scale things - again dont know why
        flux = flux * scale
        ivar = ivar / scale**2
        # mask - not sure wyh its zero
        mask  = np.zeros(flux.shape).astype(int)
        # set up the Spectra object for this camera
        spec = Spectra(bands=[band],
                       wave={band: wave},
                       flux={band: flux}, ivar={band: ivar},
                       resolution_data={band: resolution[band]},
                       mask={band: mask},
                       fibermap=spectra_fibermap,
                       #meta=meta,
                       single=True
                      )
        if specdata is None :
            specdata = spec
        else :
            specdata.update(spec)

        if debug:
            plt.clf()
            plt.plot(wave_simspec, flux_this_petal[nfiber_to_plot, :], 'k-', label='truth')
            plt.plot(wave, flux[nfiber_to_plot, :], 'x', label=f'input {band}')
            #for band_ in specdata.bands:
            #    plt.plot(specdata.wave[band_], specdata.flux[band_][nfiber_to_plot, :],
            #             '.', label=f'Spectra {band_}')
            plt.plot(specdata.wave[band], specdata.flux[band][nfiber_to_plot, :],
                         '.', label=f'Spectra {band}')

            #plt.ylim(-0.5, 2.5)
            plt.legend(bbox_to_anchor=(1,1))
            plt.xlabel('wavelength')
            plt.ylabel('flux')
            plt.title(f'nfiber in this petal: {nfiber_to_plot}')
            plt.show()
    # write out the file
    if skyerr != 0:
        spectra_filename = f'{spectra_filename}_skyerr{skyerr}.fits'
    else:
        spectra_filename = f'{spectra_filename}.fits'
    desispec.io.write_spectra(spectra_filename, specdata)

    logging.info(f'## ** wrote {spectra_filename}')
# ------------------------------------------------------------------------
fibermap_by_petal = fibermap.group_by('PETAL_LOC')

wave_simspec = simspec_hdul[0].data
flux_all_fibers = simspec_hdul[1].data

# some extra stuff for specsim
cols_to_add = ['NIGHT', 'EXPID', 'TILEID']
data_to_add = []
for col in cols_to_add:
    data_to_add.append(np.int32(hdul_fibermap[1].header[col]))

# work with each petal
petal_inds = range(10)
if debug: petal_inds = [0]

logging.info(f'looping over {len(petal_inds)} petals\n')
for petal_ind in petal_inds:
    start1 = time.time()
    fibermap_this_petal = fibermap_by_petal.groups[petal_ind]
    petal = np.unique(fibermap_this_petal['PETAL_LOC'].quantity)[0]
    logging.info(f'## working with petal {petal}')

    fibers_this_petal = fibermap_this_petal['FIBER'].quantity
    flux_this_petal = flux_all_fibers[fibers_this_petal, :] # assuming row index = fiber number NEED TO UDPATE

    nfiber = 400
    if debug:
        plt.plot(wave_simspec, flux_this_petal[nfiber, :])
        plt.xlabel('wavelength')
        plt.ylabel('flux')
        #plt.ylim(0, 2.5)
        plt.title(f'nfiber {fibers_this_petal[nfiber+1]}; petal {petal_ind}')
        plt.show()

    # add missing columns
    spectra_fibermap_this_petal = desispec.io.util.add_columns(fibermap_this_petal,
                                                               cols_to_add,
                                                               data_to_add
                                                              )
    # set up the simulator object
    simulator = desisim.simexp.simulate_spectra(wave=wave_simspec,
                                                flux=flux_this_petal,
                                                fibermap=fibermap_this_petal,
                                                obsconditions=obsconditions_dict
                                               )
    # set up the output filename
    spectra_filename = f'{path}spectra-{petal_ind}-{expid}'

    if debug: spectra_filename += f'_exptime-{simulator.observation.exposure_time.value}'

    # save the spectra for this petal
    save_petal_spectra(simulator=simulator,
                           flux_this_petal=flux_this_petal,
                           spectra_fibermap=spectra_fibermap_this_petal,
                           spectra_filename=spectra_filename,
                           seed=10,
                           skyerr=skyerr,
                           nfiber_to_plot=nfiber
                          )
    logging.info(f'## time taken for petal {petal_ind}: {(time.time() - start1)/60: .2f} min')

    # cframe
    nspec_per_petal = 500
    waves = {}
    resolution = {}
    for camera in simulator.instrument.cameras:
        waves[camera.name] = (camera.output_wavelength.to(u.Angstrom).value.astype(np.float32))
        maxbin = len(waves[camera.name])
        cframe_observedflux = np.zeros((nspec_per_petal, 3, maxbin))  # calibrated object flux
        cframe_ivar = np.zeros((nspec_per_petal, 3, maxbin)) # inverse variance of calibrated object flux
        cframe_rand_noise = np.zeros((nspec_per_petal, 3 ,maxbin)) # random Gaussian noise to calibrated flux

        R = Resolution(camera.get_output_resolution_matrix())
        resolution[camera.name] = np.tile(R.to_fits_array(), [nspec_per_petal, 1, 1])
        #print(resolution[camera.name].shape)

    fluxunits = 1e-17 * u.erg / (u.s * u.cm ** 2 * u.Angstrom)
    for j in range(nspec_per_petal):
        for i, output in enumerate(simulator.camera_output):
            num_pixels = len(output)
            # Get results for our flux-calibrated output file.
            cframe_observedflux[j, i, :num_pixels] = 1e17 * output['observed_flux'][:,0]
            cframe_ivar[j, i, :num_pixels] = 1e-34 * output['flux_inverse_variance'][:,0]

            # Fill brick arrays from the results.
            camera = output.meta['name']

            # Use the same noise realization in the cframe and frame, without any
            # additional noise from sky subtraction for now.
            cframe_rand_noise[j, i, :num_pixels] = 1e17 * (
                output['flux_calibration'][:,0] * output['random_noise_electrons'][:,0])

    #
    armName = {"b":0, "r":1, "z":2}

    for channel in 'brz':
        start = min(fibermap_this_petal['FIBER'])
        end = max(fibermap_this_petal['FIBER'])+1
        #print(start, end)

        num_pixels = len(waves[channel])
        camera = "{}{}".format(channel, petal_ind)
        # Write cframe file
        cframeFileName = desispec.io.findfile("cframe", night, int(expid), camera)
        cframeFlux = cframe_observedflux[:, armName[channel],:num_pixels] + \
                                    cframe_rand_noise[:, armName[channel],:num_pixels]
        cframeIvar = cframe_ivar[:, armName[channel], :num_pixels]

        name = cframeFileName.split('/')[-1]
        cframeFileName = f'{path}/{name}'

        meta = obsconditions_dict
        meta['CAMERA'] = camera
        cframe = Frame(waves[channel], cframeFlux, cframeIvar, \
                       resolution_data=resolution[channel],
                       spectrograph=petal_ind,
                       fibermap=fibermap_this_petal,
                       meta=meta
                      )
        compute_and_append_frame_scores(cframe)
        desispec.io.frame.write_frame(cframeFileName, cframe)

    desisim.specsim._simulators.clear()
    logging.info('done with this petal.\n\n')

logging.info('## all done')
logging.info(f'## time taken: {(time.time() - start0)/60: .2f} min')