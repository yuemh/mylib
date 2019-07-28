from __future__ import print_function
import os, sys
import numpy as np
from astropy.table import Table
from astropy.coordinates import SkyCoord
import requests
from PIL import Image
from io import BytesIO
import pylab

class image_cutout(object):

    def __init__(self, image, north=0):
        self.image = image

        ymax, xmax = image.shape
        self.ymax = ymax
        self.xmax = xmax
        self.size = np.sqrt(ymax**2 + xmax**2)

    def plot_pa(self, ax, theta):
        #fig, ax = plt.subplots()
        #ax.imshow(self.image, interpolation='None', origin='lower', cmap='gray')

        theta = 0
        line = np.linspace(-self.size/2, self.size/2, 100)
        x_to_plot = line * np.cos(theta * np.pi / 180)
        y_to_plot = line * np.sin(theta * np.pi / 180)

        ax.plot(x_to_plot, y_to_plot, 'w-')

    def measure_pa(self):
        do_sth = 1

class ps1_cutout_downloader(object):
    def __init__(self):
        self.pixscale = 0.25

    def geturl(self, ra, dec, size=240, output_size=None, filters="grizy",\
               format="fits", color=False):

        """Get URL for images in the table

        ra, dec = position in degrees
        size = extracted image size in pixels (0.25 arcsec/pixel)
        output_size = output (display) image size in pixels (default = size).
                      output_size has no effect for fits format images.
        filters = string with filters to include
        format = data format (options are "jpg", "png" or "fits")
        color = if True, creates a color image (only for jpg or png format).
                Default is return a list of URLs for single-filter grayscale images.
        Returns a string with the URL
        """

        # disable the function of downloading multiple filters at once
        if len(filters)>1:
            color = True

        if color and format == "fits":
            raise ValueError("color images are available only for jpg or png formats")
        if format not in ("jpg", "png", "fits"):
            raise ValueError("format must be one of jpg, png, fits")
        table = self.getimages(ra, dec, size=size, filters=filters)
        url = ("https://ps1images.stsci.edu/cgi-bin/fitscut.cgi?"
               "ra={ra}&dec={dec}&size={size}&format={format}").format(**locals())
        if output_size:
            url = url + "&output_size={}".format(output_size)
        # sort filters from red to blue
        flist = ["yzirg".find(x) for x in table['filter']]
        table = table[np.argsort(flist)]

        if color:
            if len(table) > 3:
                # pick 3 filters
                table = table[[0, len(table) // 2, len(table) - 1]]
            for i, param in enumerate(["red", "green", "blue"]):
                url = url + "&{}={}".format(param, table['filename'][i])
        else:
            urlbase = url + "&red="
            url = []
            for filename in table['filename']:
                url.append(urlbase + filename)
        return url

    def getimages(self, ra, dec, size=240, filters="grizy"):
        """Query ps1filenames.py service to get a list of images

        ra, dec = position in degrees
        size = image size in pixels (0.25 arcsec/pixel)
        filters = string with filters to include
        Returns a table with the results
        """

        service = "https://ps1images.stsci.edu/cgi-bin/ps1filenames.py"
        url = ("{service}?ra={ra}&dec={dec}&size={size}&format=fits"
               "&filters={filters}").format(**locals())
        table = Table.read(url, format='ascii')
        return table

    def download_images(self, ra, dec, size=240, output_size=None, \
                        filters="i", format="fits", color=False,\
                        outname_root='J', outname_precision=2):

        size = int(size)
        url = self.geturl(ra, dec, size=size, output_size=output_size, \
                               filters=filters, format=format, color=color)

        coord = SkyCoord(ra=ra, dec=dec, unit='deg')
        coordstr = coord.to_string('hmsdms', sep='', precision=outname_precision)

        outname = outname_root + coordstr.replace(' ', '') + '.%s.%s'%(filters, format)
        os.system('wget -O %s "%s"'%(outname, url[0]))


class legacy_cutout_downloader(object):
    def __init__(self):
        dummy = 1
        self.pixscale = 0.262

    def geturl(self, ra, dec, size, filters='grz', pixscale=0.262,
               layer='mzls+bass-dr6', format='fits'):
        """

        :param ra: float, ra of the object, in deg
        :param dec: float, dec of the object, in deg
        :param size: int, size of the image, in pixels
        :param filters: str, the filters of the images (can be g, r, z)
        :param pixscale: float, the size of the pixel of the output image
        :param layer: which layer to download. Can be 'decals-dr7', 'mzls+bass-dr6', or the two plus '-resid' for
                    residual image
        :param format: can be 'jpg', 'fits'
        :return: the url for downloading
        """
        service = "http://legacysurvey.org/viewer/"
        url = ["{service}{format}-cutout?ra={ra}&dec={dec}"
               "&size={size}&layer={layer}&pixscale={pixscale}"
               "&bands={filters}".format(**locals())]

        #table = Table.read(url, format='ascii')
        return url

    def download_images(self, ra, dec, size=200, \
                        filters='z', layer='mzls+bass-dr6', format='fits', \
                        outname_root='J', outname_precision=2):
        """

        :param ra: float, ra of the object, in deg
        :param dec: float, dec of the object, in deg
        :param size: int, size of the image, in pixels
        :param output_filename: str, the filename(s) to save
        :param filters:
        :param layer:
        :param fmt:
        :return:
        """
        size = int(size)
        url = self.geturl(ra, dec, size, filters=filters, pixscale=self.pixscale,\
                          layer=layer, format=format)

        coord = SkyCoord(ra=ra, dec=dec, unit='deg')
        coordstr = coord.to_string('hmsdms', sep='', precision=outname_precision)
        outname = outname_root + coordstr.replace(' ', '') + '.%s.%s' % (filters, format)

        os.system('wget -O %s "%s"' % (outname, url[0]))

### make global variables
PS1_cutout_downloader = ps1_cutout_downloader()
Legacy_cutout_downloader = legacy_cutout_downloader()

#PS1_cutout_downloader.download_images(0, 0, 20)
class des_cutout_downloader(object):
    def __init__(self):
        dummy = 1