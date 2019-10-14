from __future__ import print_function
import os, sys
import numpy as np
from astropy.table import Table
from astropy.coordinates import SkyCoord
import requests
from PIL import Image
from io import BytesIO
import pylab

from pyvo.dal import sia
#DEF_ACCESS_URL = "https://datalab.noao.edu/sia/des_dr1"
#svc = sia.SIAService(DEF_ACCESS_URL)

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

    def geturl(self, ra, dec, size, filters="grizy", output_size=None,\
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
        sizepix = int(size / self.pixscale)

        # disable the function of downloading multiple filters at once
        if len(filters)>1:
            color = True

        if color and format == "fits":
            raise ValueError("color images are available only for jpg or png formats")
        if format not in ("jpg", "png", "fits"):
            raise ValueError("format must be one of jpg, png, fits")
        table = self.getimages(ra, dec, sizepix=sizepix, filters=filters)
        url = ("https://ps1images.stsci.edu/cgi-bin/fitscut.cgi?"
               "ra={ra}&dec={dec}&size={sizepix}&format={format}").format(**locals())
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

    def getimages(self, ra, dec, sizepix, filters):
        """Query ps1filenames.py service to get a list of images

        ra, dec = position in degrees
        size = image size in pixels (0.25 arcsec/pixel)
        filters = string with filters to include
        Returns a table with the results
        """

        service = "https://ps1images.stsci.edu/cgi-bin/ps1filenames.py"
        url = ("{service}?ra={ra}&dec={dec}&size={sizepix}&format=fits"
               "&filters={filters}").format(**locals())
        table = Table.read(url, format='ascii')
        return table

    def download_images(self, ra, dec, size=6, output_size=None, \
                        filters="i", format="fits", color=False,\
                        outdir='.', outname_root='J', outname_precision=2):

#        size = int(size)
        url = self.geturl(ra, dec, size=size, output_size=output_size, \
                               filters=filters, format=format, color=color)

        coord = SkyCoord(ra=ra, dec=dec, unit='deg')
        coordstr = coord.to_string('hmsdms', sep='', precision=outname_precision)

        outname = outdir + '/' + outname_root\
                + coordstr.replace(' ', '') + '.%s.%s'%(filters, format)
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
        sizepix = int(size / self.pixscale)
        service = "http://legacysurvey.org/viewer/"
        url = ["{service}{format}-cutout?ra={ra}&dec={dec}"
               "&size={sizepix}&layer={layer}&pixscale={pixscale}"
               "&bands={filters}".format(**locals())]

        #table = Table.read(url, format='ascii')
        return url

    def download_images(self, ra, dec, size=6, \
                        filters='z', layer='dr8', format='fits', \
                        outdir='.', outname_root='J', outname_precision=2):
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
#        size = int(size)
        url = self.geturl(ra, dec, size, filters=filters, pixscale=self.pixscale,\
                          layer=layer, format=format)

        coord = SkyCoord(ra=ra, dec=dec, unit='deg')
        coordstr = coord.to_string('hmsdms', sep='', precision=outname_precision)
        outname = outdir + '/' + outname_root\
                + coordstr.replace(' ', '') + '.%s.%s' % (filters, format)

        os.system('wget -O %s "%s"' % (outname, url[0]))


class des_cutout_downloader(object):
    def __init__(self):
        dummy = 1
        self.url_dr1 = "https://datalab.noao.edu/sia/des_dr1"
        self.svc_dr1 = sia.SIAService(self.url_dr1)

    def geturl(self, ra, dec, size, filters):

        fov = size / 3600
        svc = self.svc_dr1
        imgTable = svc.search((ra,dec), (fov/np.cos(dec*np.pi/180), fov), verbosity=2).table
        sel0 = imgTable['obs_bandpass'].astype(str)==filters

        sel = sel0 & ((imgTable['proctype'].astype(str)=='Stack')\
                      & (imgTable['prodtype'].astype(str)=='image')) # basic selection
        Table = imgTable[sel]
        if (len(Table)>0):
            row = Table[np.argmax(Table['exptime'].data.data.astype('float'))] # pick image with longest exposure time
            url = row['access_url'].decode() # get the download URL

        else:
            print(imgTable)
            print('No image available.')
            url = None

        return url

    def download_images(self, ra, dec, size=6, filters='z',\
                        outdir='.', outname_root='J', outname_precision=2):
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
        url = self.geturl(ra, dec, size, filters)

#        print(url)
        coord = SkyCoord(ra=ra, dec=dec, unit='deg')
        coordstr = coord.to_string('hmsdms', sep='', precision=outname_precision)
        outname = outdir + '/' + outname_root\
                + coordstr.replace(' ', '') + '.%s.fits' % (filters)

#        print(outname)
#        print('wget -O %s "%s"' % (outname, url))

        if url:
            os.system('wget -O %s "%s"' % (outname, url))
        else:
            print('No image available.')

### make global variables
PS1_cutout_downloader = ps1_cutout_downloader()
Legacy_cutout_downloader = legacy_cutout_downloader()
DES_cutout_downloader = des_cutout_downloader()

#DES_cutout_downloader.download_images(10, 0)

#PS1_cutout_downloader.download_images(0, 0, 20)

