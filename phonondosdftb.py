#!/usr/bin/env python

import sys
import numpy as np
import matplotlib as mpl 
import matplotlib.pyplot as plt 
from matplotlib import ticker, gridspec

#import util as u

#eV2wn = u.energy_conversion['wn'] / u.energy_conversion['eV']


#
# Lorentzian gamma : FWHM
# Gaussian sigma : Std Dev
#

def lorentzian(SpecRange, freq, Gamma):
    
    lsfunc = 1 / np.pi * Gamma / (Gamma**2 + (SpecRange - freq)**2)

    return lsfunc


def damped_lorentzian(SpecRange, freq, Gamma):
    
    lsfunc = 1 / np.pi * SpecRange**2 * Gamma / ((Gamma * SpecRange)**2 + (SpecRange**2 - freq**2)**2)

    return lsfunc


def gaussian(SpecRange, freq, Sigma):

    expfac = -0.5 * ((SpecRange - freq) / Sigma)**2
    lsfunc = 1 / Sigma * np.sqrt(2 * np.pi) * np.exp(expfac)

    return lsfunc


def calcspecden(freqs, ls="lor", widths=None, step=0.5):

    # Width in wavenumbers, default value is 5.0 cm^-1 for each transitions
    # widths is an array, so it possible to have different broadenings for each
    # transition
    if not widths:
        widths = np.array(len(freqs) * [5.0])

    if ls == "gau":
        ls = gaussian

    elif ls == "lor":
        ls = lorentzian

    elif ls == "dlor":
        ls = damped_lorentzian

    SpecRange = np.arange(0, 3500, step)
    SpecDen = np.zeros(len(SpecRange))

    for i in range(len(freqs)):
        SpecDen += ls(SpecRange, freqs[i], widths[i])

    result = np.c_[SpecRange, SpecDen]

    return result


def plot(data, unit="wn", figname="specden.svg"):

    minloc = ticker.AutoMinorLocator(4)

    if unit == "eV":
        x = data[:,1]
        xlabel = 'E / eV'
        xmaj = ticker.MultipleLocator(0.05)

    elif unit == "wn":
        x = data[:,0]
        xlabel = r'$\tilde{\nu}$ / cm$^{-1}$'
        xmaj = ticker.MultipleLocator(500)

    specden = data[:,-1]

    # Set up a plot
    fig = plt.figure() 
    gs = gridspec.GridSpec(1, 1)
    gs.update(wspace=0.1, hspace=0.1)

    # Plot options
    ax0 = plt.subplot(gs[0])
    ax0.set_xlabel(xlabel)
    ax0.xaxis.set_major_locator(xmaj)
    ax0.xaxis.set_minor_locator(minloc)
    ax0.yaxis.set_minor_locator(minloc)
    ax0.set_ylabel(r'$J(\omega)$ / cm$^{-1}$')
    ax0.yaxis.set_label_coords(-0.1, 0.5)
    ax0.tick_params(axis='both', which='major', length=7, pad=10)
    ax0.tick_params(axis='both', which='minor', length=3)
    ax0.plot(x, specden)
    ax0.relim()
    ax0.autoscale_view(True,True,True)
    ax0.set_aspect(0.65/ax0.get_data_ratio())
    plt.minorticks_on()

    plt.minorticks_on()
    mpl.rcParams.update({'font.size': 18})
    plt.gcf().subplots_adjust(bottom=0.15)
    plt.savefig(figname, dpi=600)
    plt.close()
    # plt.show()

    return


if __name__ == '__main__':

    import textwrap
    import argparse as arg

    parser = arg.ArgumentParser(description='Convolutes Stick Spectra',
                                formatter_class=arg.ArgumentDefaultsHelpFormatter)
    
    #
    # Input files
    #
    inp = parser.add_argument_group("Input Data")
    
    inp.add_argument('-i', '--input',
                     default="results.txt", type=str, dest="InputFile",
                     help='''Input file''')
    
    #
    # Spectra Options
    #
    spec = parser.add_argument_group("Spectra Convolution Options")
    
    spec.add_argument('--ls',
                     default="lor", type=str, choices=["gau", "lor", "dlor"],
                     dest="LineShape", help='''Spectral LineShape.''')
    
    spec.add_argument('--lw',
                     default=[5.0], type=float, nargs='+', dest="LineWidth",
                     help='''Spectral LineWidth in wavenumbers (gamma for Lorentzian,
                     sigma for Gaussian LineShape.''')
    
    spec.add_argument('--unit',
                     default="wn", type=str, choices=["eV", "wn"], dest="SpecUnit",
                     help='''X axis unit for plotting Spectra.''')

    #
    # Output Options
    #
    out = parser.add_argument_group("Output Options")

    out.add_argument('-o', '--out',
                     default="Spec", type=str, dest="OutPref",
                     help='''Output Prefix for files''')

    out.add_argument('--figext',
                     default=None, type=str, choices=["svg", "png", "eps"],
                     dest="FigExt", help='''Format for image output''')
    
    #
    # Parse and create the Options Dictionary
    #
    args = parser.parse_args()
    Opts = vars(args)

    if Opts['OutPref']:
        Opts['OutPref'] = Opts['OutPref'] + "."
    #
    # Convolute Spectra
    #
    data = np.loadtxt(Opts['InputFile'])
    widths = Opts['LineWidth'] * len(data)
    specden = calcspecden(data, ls=Opts['LineShape'], widths=widths)

    #
    # Convoluted Spectral Density
    #
    titles = ("E (cm^-1)", "E (eV)", "Spectral Density (cm^-1)")
    header = ("\nCalculated Spectral Density\n"
              "Lineshape: %s\n"
              "Linewidths (cm^-1): %s \n\n" % (Opts['LineShape'], widths[0]))

    header1 = "%9s %10s %25s\n" % (titles)
    fmt = "%11.4f %17.8e"
    np.savetxt(Opts['OutPref'] + "specden.txt", specden, fmt=fmt, header=header+header1)

    if Opts['FigExt']:

         figname = Opts['OutPref'] + "specden." + Opts['FigExt'] 
         plot(specden, unit=Opts['SpecUnit'], figname=figname)

    pass
