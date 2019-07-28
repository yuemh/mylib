
import os, sys
import numpy as np
from pyraf import iraf

# read directories
dir_data = os.getcwd()

# load iraf packages

iraf.load('noao')
iraf.load('imred')
iraf.load('ccdred')
iraf.load('onedspec')
iraf.load('twodspec')
iraf.load('longslit')
iraf.load('apextract')

def filename_prep(direct, suffix='.fits'):
    bias_dir = direct + '/bias'
    flat_dir = direct + '/flat'
    arc_dir = direct + '/arc'
    sci_dir = direct + '/sci'
    std_dir = direct + '/std'

    # biases
    biases = []
    for root, namedir, names in os.walk(bias_dir):
        for name in names:
            if not name[-len(suffix):]==suffix:
                continue
            else:
                biases.append(root + '/' + name)

    # flats
    flats = []
    for root, namedir, names in os.walk(flat_dir):
        for name in names:
            if not name[-len(suffix):]==suffix:
                continue
            else:
                flats.append(root + '/' + name)

    # arcs
    arcs = []
    for root, namedir, names in os.walk(arc_dir):
        for name in names:
            if not name[-len(suffix):]==suffix:
                continue
            else:
                arcs.append(root + '/' + name)

    # science frames
    scis = []
    for root, namedir, names in os.walk(sci_dir):
        for name in names:
            if not name[-len(suffix):]==suffix:
                continue
            else:
                scis.append(root + '/' + name)

    # standard stars
    stds = []
    for root, namedir, names in os.walk(std_dir):
        for name in names:
            if not name[-len(suffix):]==suffix:
                continue
            else:
                stds.append(root + '/' + name)

    files_dict = {'bias': biases,\
                  'flat': flats,\
                  'arc': arcs,\
                  'sci': scis,\
                  'std': stds}

    # create database directory

    database_dir = './database'
    database_dir_arc = database_dir + '/id./arc'
    database_dir_sci = database_dir + '/id./sci'
    database_dir_std = database_dir + '/id./std'

    database_dir_list = [database_dir_arc, database_dir_sci, database_dir_std]

    for direct in database_dir_list:
        if not os.path.exists(direct):
            os.system('mkdir -p ' + direct)

    return files_dict

def images_with_suffix(namelist):
    trimmed = [name[:-5] + '_t.fits' for name in namelist]
    trimmed_bias = [name[:-5] + '_tb.fits' for name in namelist]
    trimmed_bias_flat = [name[:-5] + '_tbf.fits' for name in namelist]
    extracted = [name[:-5] + '_tbfe.fits' for name in namelist]
    extracted_dispcor = [name[:-5] + '_tbfed.fits' for name in namelist]

    return [trimmed, trimmed_bias, trimmed_bias_flat, extracted, extracted_dispcor]

def image_list_to_str(namelist):
    str_to_return = ''
    for name in namelist:
        str_to_return = str_to_return + name + ','

    return str_to_return[:-1]

def trim_image_iraf(imagelist, outputlist, cut_range):
    if not len(imagelist)==len(outputlist):
        raise ValueError('The number of input and output images is different')

    for index in range(len(imagelist)):
        input_image = imagelist[index] + cut_range
        output_image = outputlist[index]

        if os.path.exists(output_image):
            os.system('rm ' + output_image)

        iraf.imcopy(input_image, output_image)

    return outputlist

def combine_bias(bias_images, output_bias_image, cut_range=''):
    bias_images_trimmed = [name[:-5] + '_t.fits' for name in bias_images]
    output = trim_image_iraf(bias_images, bias_images_trimmed, cut_range)

    image_to_combine_str = image_list_to_str(bias_images_trimmed)

    iraf.imcombine(input=image_to_combine_str, output=output_bias_image,\
                   combine="average", project="no")

    for index in range(len(bias_images_trimmed)):
        os.system('rm ' + bias_images_trimmed[index])

    return output_bias_image

def combine_flat(flat_images, output_flat_image, output_nflat_image, bias_image,\
                 cut_range=''):
    flat_images_trimmed = [name[:-5] + '_t.fits' for name in flat_images]
    flat_images_trimmed_bias = [name[:-5] + '_tb.fits' for name in flat_images]

    output = trim_image_iraf(flat_images, flat_images_trimmed, cut_range)

    for index in range(len(flat_images_trimmed)):
        tname = flat_images_trimmed[index]
        tbname = flat_images_trimmed_bias[index]
        if os.path.exists(tbname):
            os.system('rm ' + tbname)
        iraf.imarith(tname, '-', bias_image, tbname)

    image_to_combine_str = image_list_to_str(flat_images_trimmed_bias)

    if os.path.exists(output_flat_image):
        os.system('rm ' + output_flat_image)

    iraf.imcombine(image_to_combine_str, output_flat_image, combine='median')
    iraf.response(calibrat=output_flat_image,\
                  normaliz=output_flat_image, response=output_nflat_image,\
                  interac='no', functio="spline3", order=16, low_rej=3, high_re=3)

    for index in range(len(flat_images_trimmed)):
        os.system('rm ' + flat_images_trimmed[index])
        os.system('rm ' + flat_images_trimmed_bias[index])

    return output_nflat_image

def do_bias_flat(imagelist, bias_image, flat_image, b_imagelist, bf_imagelist):
    for index in range(len(imagelist)):
        name = imagelist[index]
        b_name = b_imagelist[index]
        bf_name = bf_imagelist[index]

        if os.path.exists(b_name):
            os.system('rm ' + b_name)

        if os.path.exists(bf_name):
            os.system('rm ' + bf_name)

        iraf.imarith(name, '-', bias_image, b_name)
        iraf.imarith(b_name, '/', flat_image, bf_name)
        os.system('rm ' + b_name)

def wavelength_calib(arc_images, bias_image, flat_image, ref_spec,\
                     cut_range='', filename_only=False):
    arc_images_trimmed, arc_images_trimmed_bias,\
        arc_images_trimmed_bias_flat, arc_extracted,\
            arc_extracted_dispcor = images_with_suffix(arc_images)

    if filename_only:
        return arc_extracted

    else:
        # trim
        trim_image_iraf(arc_images, arc_images_trimmed, cut_range=cut_range)

        # bias & flat
    #    do_bias_flat(arc_images_trimmed, bias_image, flat_image,\
    #                 arc_images_trimmed_bias, arc_images_trimmed_bias_flat)

        # extract the 1-D spec
        arc_images_str = image_list_to_str(arc_images_trimmed)
        arc_extracted_str = image_list_to_str(arc_extracted)

        iraf.apall(arc_images_str,\
                   ref=ref_spec, output=arc_extracted_str,\
                   format="multispec", recenter="no", trace="no", back="no", intera="no")

        # identify
        iraf.identify(arc_extracted[0], coordlist="linelists$henear.dat")

        # re-identify
        # iraf.reidentify(output_spec[0], output_spec_str)

        return arc_extracted

def standard_star_step1(std_images, bias_image, flat_image,\
                        cut_range='', filename_only=False,\
                       line=600):
    std_images_trimmed, std_images_trimmed_bias,\
        std_images_trimmed_bias_flat, std_extracted,\
            std_extracted_dispcor = images_with_suffix(std_images)

    if filename_only:
        return std_extracted
    else:
        std_images_str = image_list_to_str(std_images_trimmed_bias_flat)
        std_extracted_str = image_list_to_str(std_extracted)
        std_extracted_dispcor_str = image_list_to_str(std_extracted_dispcor)

        # trim
        trim_image_iraf(std_images, std_images_trimmed, cut_range=cut_range)

        # bias & flat
        do_bias_flat(std_images_trimmed, bias_image, flat_image,\
                     std_images_trimmed_bias, std_images_trimmed_bias_flat)

        # Get the aperture
        for index in range(len(std_images_trimmed_bias_flat)):
            std_image = std_images_trimmed_bias_flat[index]
            iraf.apall(std_image, output=std_extracted[index], intera="yes",\
                      line=line, t_function="spline3")

        return std_extracted

def standard_star_step2(std_images, arc_extracted, std_name,\
                        caldir="onedstds$irscal/", sens_output='sens.fits',\
                        filename_only=False):
    std_images_trimmed, std_images_trimmed_bias,\
        std_images_trimmed_bias_flat, std_extracted,\
            std_extracted_dispcor = images_with_suffix(std_images)

    if filename_only:
        return sens_output

    else:
        # Edit wavelength info
        for index in range(len(std_extracted)):
            iraf.hedit(std_extracted[index], "REFSPEC1", add="yes", value=arc_extracted)

            if os.path.exists(std_extracted_dispcor[index]):
                os.system('rm '+std_extracted_dispcor[index])

            iraf.dispcor(std_extracted[index], std_extracted_dispcor[index])
            iraf.setairmass(std_extracted_dispcor[index])

        # measure the sensfunc
        print('standrad')
        iraf.standard(std_extracted_dispcor[0], output=std_name+'.fits',\
                      caldir="onedstds$irscal/", star_name=std_name, answer="no")

        if os.path.exists(sens_output):
            os.system('rm ' + sens_output)
        iraf.sensfunc(std_name+'.fits', sens_output, ignorea="yes")

        return sens_output

def science_frames_step1(sci_images, bias_image, flat_image,\
                        cut_range='', lower=-5, upper=5, width=5, line=600,\
                        filename_only=False):
    sci_images_trimmed, sci_images_trimmed_bias,\
        sci_images_trimmed_bias_flat, sci_extracted,\
            sci_extracted_dispcor = images_with_suffix(sci_images)

    if filename_only:
        return sci_extracted
    else:
        # trim
        trim_image_iraf(sci_images, sci_images_trimmed, cut_range=cut_range)

        # bias & flat
        do_bias_flat(sci_images_trimmed, bias_image, flat_image,\
                     sci_images_trimmed_bias, sci_images_trimmed_bias_flat)

        sci_images_str = image_list_to_str(sci_images_trimmed_bias_flat)
        sci_extracted_str = image_list_to_str(sci_extracted)

        # Get the aperture
        iraf.apall(sci_images_str, output=sci_extracted_str, intera="yes",\
                  lower=lower, upper=upper, width=width, t_function='spline3',\
                  line=line)
        return sci_extracted

def science_frames_step2(sci_images, arc_extracted):
    sci_images_trimmed, sci_images_trimmed_bias,\
        sci_images_trimmed_bias_flat, sci_extracted,\
            sci_extracted_dispcor = images_with_suffix(sci_images)

    sci_final = [name[:-5] + '_final.fits' for name in sci_images]

    sci_images_str = image_list_to_str(sci_images_trimmed_bias_flat)
    sci_extracted_str = image_list_to_str(sci_extracted)
    sci_extracted_dispcor_str = image_list_to_str(sci_extracted_dispcor)

    for index in range(len(sci_extracted)):
        iraf.hedit(sci_extracted[index], "REFSPEC1", add="yes", value=arc_extracted)

        if os.path.exists(sci_extracted_dispcor[index]):
            os.system('rm '+sci_extracted_dispcor[index])

        iraf.dispcor(sci_extracted[index], sci_extracted_dispcor[index])
        iraf.setairmass(sci_extracted_dispcor[index])

        if os.path.exists(sci_final[index]):
            os.system('rm '+sci_final[index])

        ## work on final filename ##
        final_idx = 0
        while True:
            if os.path.exists(sci_final[index][:-5]+'_%d.fits'%(final_idx)):
                final_idx += 1
            else:
                break

        final_name = sci_final[index][:-5]+'_%d.fits'%(final_idx)
        print(final_name)

        iraf.calibrate(sci_extracted_dispcor[index], final_name)

        os.system('rm ' + sci_extracted_dispcor[index])
        os.system('rm ' + sci_extracted[index])
        os.system('rm ' + sci_images_trimmed_bias_flat[index])
        os.system('rm ' + sci_images_trimmed_bias[index])
        os.system('rm ' + sci_images_trimmed[index])

def full_reduction(filenames=[], arcnames=[], std_name='', cut_range='', step=0,\
                  extract_params={}):
    """
    Do the full spectra reduction.
    Step 1: Combine flat and bias
    Step 2: Wavelength calibration with arc frame
    Step 3: Flux calibration with standard star
    Step 4: Generate the calibrated images.
    """
    files = filename_prep('.', suffix='.fits')

    if len(filenames)==0:
        filenames = files['sci']
    if len(arcnames)==0:
        arcnames = files['arc']

    if not os.path.exists('./bias.fits'):
        bias = combine_bias(files['bias'], './bias.fits',\
                           cut_range=cut_range)

    if not os.path.exists('./flat.fits'):
        flat = combine_flat(files['flat'], './rflat.fits',\
                        './flat.fits', './bias.fits',\
                           cut_range=cut_range)

    # step 0: extract the standard star

    std_extracted = standard_star_step1(files['std'],\
                            './bias.fits', './flat.fits',\
                            cut_range=cut_range, filename_only=(step>0),\
                            **extract_params)

    # step 1: wavelength calibration

    arc_extracted = wavelength_calib(arcnames,\
                     bias_image='./bias.fits', flat_image='./flat.fits',\
                     ref_spec=std_extracted, cut_range=cut_range,\
                    filename_only=(step>1))

    # step 2: get the response curve

    sens = standard_star_step2(files['std'],\
                              arc_extracted=arc_extracted[0],
                              std_name=std_name, filename_only=(step>2))

    # step 3: extract the science frame

    science_frames_step1(filenames, './bias.fits', './flat.fits', cut_range=cut_range,\
                        **extract_params)

    # step 4: wavelength / flux calibration

    science_frames_step2(filenames, arc_extracted[0])

def default_params_instru(instrument, seeing=1, line=600):
    known_instru = ['LDSS3', 'RED']

    if not instrument in known_instru:
        raise ValueError('Unknown instrument name %s. Accepted instruments includes %s'\
                         %(instrument, known_instru))

    if instrument=='LDSS3':
        pixscale = 0.189
        axis = 1

    elif instrument=='RED':
        pixscale = 0.25
        axis = 1

    upper = int(seeing / pixscale)
    lower = -int(seeing / pixscale)
    width = int(seeing / pixscale)

    return {'upper': upper, 'lower': lower, 'line': line}

def pre_quick_reduction():
    files = filename_prep('.')

    if len(filenames)==0:
        filenames = files['sci']
    if len(arcnames)==0:
        filenames = files['arc']

    if not os.path.exists('./bias.fits'):
        bias = combine_bias(files['bias'], './bias.fits',\
                           cut_range=cut_range)

    if not os.path.excists('./flat.fits'):
        flat = combine_flat(files['flat'], './rflat.fits',\
                        './flat.fits', './bias.fits',\
                           cut_range=cut_range)

    if not os.path.exists('./arc.fits'):
        std_extracted = standard_star_step1(files['std'],\
                                './bias.fits', './flat.fits',\
                                cut_range=cut_range)
        arc_extracted = wavelength_calib(arcnames,\
                         bias_image='./bias.fits', flat_image='./flat.fits',\
                         ref_spec=std_extracted, cut_range=cut_range)

    if not os.path.exists('./sens.fits'):
        sens = standard_star_step2(files['std'],\
                                  arc_extracted=arc_extracted,
                                  std_name=std_name)

def quick_reduction(filename):
    """
    Do a quick reduction.
    Assume that you have the first three steps done.
    """
    science_frames_step1(files['sci'], './bias.fits', './flat.fits')
    science_frames_step2(files['sci'], './arc.fits')

def main():
    full_reduction(filenames=['./sci/ccd0173c1.fits'],\
                   std_name='hiltner600', cut_range='[50:1000, 1000:4000]')

if __name__=='__main__':
    main()
