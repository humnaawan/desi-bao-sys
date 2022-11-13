# ------------------------------------------------------------------------
# this script takes in the fibermap and simspec simulation files and
# produces the spectra and cframe files. if using the debug option, will
# just run things for one petal; all outputs will then have _debug tag.
#
# most of the code is based on quickspectra/quickquasars
# ------------------------------------------------------------------------
from astropy.io import fits
from astropy.table import Table
import astropy.units as u
import numpy as np
from settings import *

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
temp = f'create-per-pixel-spectra_{datetime.datetime.now()}'.replace(' ', '_').split('.')[0]
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
hdul_fibermap = fits.open(fibermap_fname)
fibermap = Table(hdul_fibermap[hdul_fibermap.index_of('FIBERMAP')].data)

for key in ['FIBER', 'TARGETID', 'DEVICE_LOC']:
    uniq, cts = np.unique(fibermap[key], return_counts=True)
    logging.info(f'{len(uniq)} unique {key} with {np.unique(cts)} unique counts across')

# ------------------------------------------------------------------------
# simspec data
simspec_fname = f'{path}/simspec-{expid}.fits'
hdul_simspec = fits.open(simspec_fname)
simspec = Table(hdul_simspec[hdul_simspec.index_of('TRUTH')].data)
obsconditions = Table(hdul_simspec[hdul_simspec.index_of('OBSCONDITIONS')].data)
obsconditions_dict = [dict(zip(obsconditions.colnames, row)) for row in obsconditions][0]
if debug and exptime is not None:
    logging.info(f'# updating exptime to be {exptime}s.')
    obsconditions_dict['EXPTIME'] = exptime

# update the exptime in the fibermap
fibermap['EXPTIME'][:] = obsconditions_dict['EXPTIME']
#print(f'obsconditions_dict: {obsconditions_dict}')
# ------------------------------------------------------------------------------------
def save_data_for_this_petal(simulator, nspec, fibermap, spectra_fibermap, spectra_filename,
                             seed=10, skyerr=0, debug=False):
    #
    logging.info(f'## running save_spectra_for_petal ...')
    scale = 1e17          # not exactly sure where this is coming from but this goes into the flux units
    # initialize random state
    random_state = np.random.RandomState(seed)
    # skyscale
    skyscale = skyerr * random_state.normal(size=simulator.num_fibers)

    # ---------------------------------------
    # first deal with the cframe file
    # set up the Resolution for each camera, alongside the waves for each channel (3 of them); need for cframe files
    waves, resolution = {}, {}
    for camera in simulator.instrument.cameras:
        waves[camera.name] = (camera.output_wavelength.to(u.Angstrom).value.astype(np.float32))
        maxbin = len(waves[camera.name])
        cframe_observedflux = np.zeros((nspec, 3, maxbin))      # calibrated object flux
        cframe_ivar = np.zeros((nspec, 3, maxbin))              # inverse variance of calibrated object flux
        cframe_noise = np.zeros((nspec, 3, maxbin))             # noise to calibrated flux

        R = Resolution(camera.get_output_resolution_matrix())
        resolution[camera.name] = np.tile(R.to_fits_array(), [nspec, 1, 1]) # need to check if this is right since quickgen seems to have a different way to handle this

    # store the flux etc for cframe file
    for j in range(nspec): #:
        for i, output in enumerate(simulator.camera_output):
            num_pixels = len(output)
            # get the flux
            cframe_observedflux[j, i, :num_pixels] = scale * output['observed_flux'][:, j]
            # get the variance
            cframe_ivar[j, i, :num_pixels] = (1/scale**2) * output['flux_inverse_variance'][:, j]
            # noise realization; with additional noise from sky subtraction if applicable
            cframe_noise[j, i, :num_pixels] = scale * (output['flux_calibration'][:, j] * output['random_noise_electrons'][:, j])
            if np.any(skyscale):
                cframe_noise[j, i, :num_pixels] += scale * ((output['num_sky_electrons'][:, j] * skyscale) * output['flux_calibration'][:, j])

    # now lets create the cframe file
    arm_name = {'b': 0, 'r': 1, 'z': 2}
    for channel in 'brz':
        start = min(fibermap['FIBER'])
        end = max(fibermap['FIBER'])+1

        num_pixels = len(waves[channel])
        camera = "{}{}".format(channel, petal_ind)
        # write the cframe file
        cframe_fname = desispec.io.findfile('cframe', night, int(expid), camera)
        cframe_flux_thischannel = cframe_observedflux[:, arm_name[channel], :num_pixels] + \
                                    cframe_noise[:, arm_name[channel], :num_pixels]
        cframe_ivar_thischannel = cframe_ivar[:, arm_name[channel], :num_pixels]

        name = cframe_fname.split('/')[-1]
        if skyerr != 0:
            name = name.split('.fits')[0]
            name = f'{name}_skyerr{skyerr}.fits'
        if debug:
            name = name.split('.fits')[0]
            name = f'{name}_debug.fits'
        cframe_fname = f'{path}/{name}'

        meta = obsconditions_dict
        meta['CAMERA'] = camera
        cframe = Frame(waves[channel], cframe_flux_thischannel, cframe_ivar_thischannel, \
                       resolution_data=resolution[channel],
                       spectrograph=petal_ind,
                       fibermap=fibermap,
                       meta=meta
                      )
        compute_and_append_frame_scores(cframe)
        desispec.io.frame.write_frame(cframe_fname, cframe)
        logging.info(f'## wrote {cframe_fname}')

    # -------------------------------------
    # now save the Spectra object
    specdata = None
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

        # scale things
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

    # write out the file
    if skyerr != 0:
        spectra_filename = f'{spectra_filename}_skyerr{skyerr}.fits'
    else:
        spectra_filename = f'{spectra_filename}.fits'
    desispec.io.write_spectra(spectra_filename, specdata)

    logging.info(f'## ** wrote {spectra_filename}')
# ------------------------------------------------------------------------
fibermap_by_petal = fibermap.group_by('PETAL_LOC')

wave_simspec = hdul_simspec[hdul_simspec.index_of('WAVE')].data
# lets make sure that the fluxes are sorted by fiber number
inds = np.argsort(hdul_simspec[hdul_simspec.index_of('FIBERMAP')].data['FIBER'])
flux_all_fibers = hdul_simspec[hdul_simspec.index_of('FLUX')].data[inds, :]

# some extra stuff for specsim
cols_to_add = ['NIGHT', 'EXPID', 'TILEID']
data_to_add = []
for col in cols_to_add:
    data_to_add.append(np.int32(hdul_fibermap[hdul_fibermap.index_of('FIBERMAP')].header[col]))

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
    flux_this_petal = flux_all_fibers[fibers_this_petal, :] # row index = fiber number give the sort above

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
    save_data_for_this_petal(simulator=simulator,
                             nspec=flux_this_petal.shape[0],
                             spectra_fibermap=spectra_fibermap_this_petal,
                             spectra_filename=spectra_filename,
                             fibermap=fibermap_this_petal,
                             seed=10, skyerr=skyerr, debug=debug)
    logging.info(f'## time taken for petal {petal_ind}: {(time.time() - start1)/60: .2f} min')


    desisim.specsim._simulators.clear()
    logging.info('done with this petal.\n\n')

logging.info('## all done')
logging.info(f'## time taken: {(time.time() - start0)/60: .2f} min')