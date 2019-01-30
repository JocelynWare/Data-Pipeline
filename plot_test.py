'''
Created on 29 Jan 2019

@author: Jacob Parnell
'''

import matplotlib.pyplot as plt
import glob
import astropy.io.fits as pyfits
import numpy as np
from scipy import ndimage

from flat_fielding import onedim_pixtopix_variations_single_order

## DO NOT RUN THIS CODE FOR ALL EXTRACTED SOLAR SPECTRA AND PRINT ALL PLOTS - WILL RUN OUT OF MEMORY

## loads the extracted 2-d images and places in a list - only for the quick-extracted solar images
plot_path = '/Users/Jacob/Desktop/data_for_mq/'
plot_list = sorted(glob.glob(plot_path + '*solar*' + '*quick*.fits'))

stellar_object = np.arange(0,len(plot_list))
flux_load = np.load(plot_path + 'flux.npy') .item() # added this line so that we can run this script without having to re-run the entire basic_reduction_script

for j in stellar_object:

    solar_spectrum = pyfits.getdata(plot_list[j])
    nx,ny = solar_spectrum.shape # set the size of the array to nx,ny, note that nx = 42, ny = 4096 but the image is the opposite

    nxx = np.arange(0,nx) # change the y-value to adjust the orders (default is nx)

    for i in nxx:
        one_dim_solar = solar_spectrum[i, :]
        nx = np.arange(0, ny)

        ## plot - just the 1-dim spectrum with polynomial fit (x-axis: pixels, y-axis: intensity)
        ## calculate polynomial and fit 5th order - works for second order, later orders fit seems to drop off
        z = np.polyfit(nx, one_dim_solar, 5)
        f = np.poly1d(z)
        ## calculate new x's and new y's
        x_new = np.linspace(nx[0], nx[-1], len(nx))
        y_new = f(x_new)
        ## plot
        plot_name = plot_list[j].split('/')
        plt.xlim(0, ny)
        plt.xlabel('Pixels')
        plt.ylabel('Intensity')
        plt.title(plot_name[-1].split('.')[0] + ' : ' + 'Order: ' + str(i+2))
        plt.plot(nx, one_dim_solar, x_new, y_new, 'r', linewidth=2)
        plt.show()
        plt.close()

        ## save file as the name of the graph + row, also include a new path into a new folder
        # save_plots = plot_path + 'Plots/'
        # plt.savefig(save_plots + plot_name[-1].split('.')[0] + '.png')

        ## find the pixel sensitivity
        f_flat = flux_load.values()[i]
        smoothed_flat, pix_sens = onedim_pixtopix_variations_single_order(f_flat)
        ## plot - commented out unless you want to see the pixel sensitivity
        plt.xlim(0, ny)
        plt.xlabel('Pixels')
        plt.ylabel('Relative Sensitivity')
        plt.title('Pixel Sensitivity')
        plt.plot(nx, pix_sens)
        plt.show()

        ## divide this new pixel sensitivity by the science 1-dim spectra to produce a smoother spectrum
        ## is this really smoother?
        smooth_1dim = one_dim_solar / pix_sens
        plt.xlim(0, ny)
        plt.xlabel('Pixels')
        plt.ylabel('Intensity')
        plt.title('Smoothed: ' + plot_name[-1].split('.')[0] + ' : ' + 'Order: ' + str(i+2))
        plt.plot(nx, smooth_1dim, x_new, y_new, 'r', linewidth=2)
        plt.show()
        break # added to skip all orders and all elements - just prints the 2nd order for the first element
    break # added to skip all orders and all elements - just prints the 2nd order for the first element