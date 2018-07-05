# This Python file uses the following encoding: utf-8
import sys
sys.path.insert(0, "../../")
sys.path.insert(0, "../")
sys.path.insert(0, "./")

from scipy.io import netcdf
import numpy as np
import pytest
import subprocess
import math

#import matplotlib
#matplotlib.use('Agg')
import matplotlib.pyplot as plt

from parcel import parcel

def plot_profiles(data, output_folder):

    plt.figure()
    plots    = []
    legend_l = []
    for i in range(3):
        plots.append(plt.subplot(1,3,i+1))

    plots[0].set_xlabel('T [K]')
    plots[1].set_xlabel('rv [g/kg]')
    plots[2].set_xlabel('RH')

    for ax in plots:
        ax.set_ylabel('z [m]')

    z = data.variables["z"][:]
    plots[0].plot(data.variables["T"][:]          , z)
    plots[1].plot(data.variables["r_v"][:] * 1000 , z)
    plots[2].plot(
        data.variables["RH"][:]                     , z,
        [data.variables["RH"][:].max()] * z.shape[0], z,
        [1.] * z.shape[0], z
        )

    plt.tight_layout()
    plt.savefig(output_folder + "acnv_profiles.pdf")

def plot_spectrum_m0(data, output_folder):

    ymax = 1e12
    ymin = 1

    rr = data.variables["rradii_r_wet"][:] * 1e6
    cr = data.variables["cradii_r_wet"][:] * 1e6
    ar = data.variables["aradii_r_wet"][:] * 1e6

    drr = data.variables["rradii_dr_wet"][:] * 1e6
    dcr = data.variables["cradii_dr_wet"][:] * 1e6
    dar = data.variables["aradii_dr_wet"][:] * 1e6

    for t in range(data.variables['t'].shape[0]):

        plt.clf()
        f, ax = plt.subplots()
        ax.set_xscale("log")
        ax.set_yscale("log")
        ax.set_xlabel("particle radius [$\mu m$]")
        ax.set_ylabel("concentration [$kg^{-1} \mu m^{-1}$]")
        ax.set_ylim(ymin, ymax)

        nr = data.variables['rradii_m0'][t,:] / drr
        nc = data.variables['cradii_m0'][t,:] / dcr
        na = data.variables['aradii_m0'][t,:] / dar

        plt.step(rr, nr, color='magenta', label="rain",   lw = 2.5, rasterized=False)
        plt.step(cr, nc, color='blue',    label="cloud",  lw = 2.5, rasterized=False)
        plt.step(ar, na, color='orange',  label="aerosol",lw = 2.5, rasterized=False)

        plt.axvline(x=1, color='gray')
        plt.axvline(x=25, color='gray')

        plt.legend(frameon=False)

        plt.savefig(output_folder + 'acnv_size_distr_m0_' + str("%03d" % t) + '.pdf')

def plot_spectrum_m3(data, output_folder):

    ymax = 1e-0
    ymin = 1e-6

    rr = data.variables["rradii_r_wet"][:] * 1e6
    cr = data.variables["cradii_r_wet"][:] * 1e6
    ar = data.variables["aradii_r_wet"][:] * 1e6

    drr = data.variables["rradii_dr_wet"][:] * 1e6
    dcr = data.variables["cradii_dr_wet"][:] * 1e6
    dar = data.variables["aradii_dr_wet"][:] * 1e6

    to_mass = 4./3 * math.pi * 1e3                 * 1e3
    #                          ^density of water   ^ convert to g
    for t in range(data.variables['t'].shape[0]):

        plt.clf()
        f, ax = plt.subplots()
        ax.set_xscale("log")
        ax.set_yscale("log")
        ax.set_xlabel("particle radius [$\mu m$]")
        ax.set_ylabel("mass [$g_{water} \; kg^{-1} \mu m^{-1}$]")
        ax.set_ylim(ymin, ymax)

        nr = data.variables['rradii_m3'][t,:] / drr * to_mass
        nc = data.variables['cradii_m3'][t,:] / dcr * to_mass
        na = data.variables['aradii_m3'][t,:] / dar * to_mass

        plt.step(rr, nr, color='magenta', label="rain",   lw = 2.5, rasterized=False)
        plt.step(cr, nc, color='blue',    label="cloud",  lw = 2.5, rasterized=False)
        plt.step(ar, na, color='orange',  label="aerosol",lw = 2.5, rasterized=False)

        plt.axvline(x=1, color='gray')
        plt.axvline(x=25, color='gray')

        plt.legend(frameon=False)

        plt.savefig(output_folder + 'acnv_size_distr_m3_' + str("%03d" % t) + '.pdf')

def plot_spectrum_m6(data, output_folder):

    ymax = 10
    ymin = -200
    xmin = 1e-3
    xmax = 1e3

    rr = data.variables["rradii_r_wet"][:] * 1e6
    cr = data.variables["cradii_r_wet"][:] * 1e6
    ar = data.variables["aradii_r_wet"][:] * 1e6

    drr = data.variables["rradii_dr_wet"][:] * 1e6
    dcr = data.variables["cradii_dr_wet"][:] * 1e6
    dar = data.variables["aradii_dr_wet"][:] * 1e6

    for t in range(data.variables['t'].shape[0]):

        plt.clf()
        f, ax = plt.subplots()
        ax.set_xscale("log")
        #ax.set_yscale("log")
        ax.set_xlabel("particle radius [$\mu m$]")
        ax.set_ylabel("TODO [$kg^{-1} \mu m^{-1}$]")
        ax.set_ylim(ymin, ymax)
        ax.set_xlim(xmin, xmax)

        nr = 10 * np.log10(data.variables['rradii_m6'][t,:] / 1e-18) / drr
        nc = 10 * np.log10(data.variables['cradii_m6'][t,:] / 1e-18) / dcr
        na = 10 * np.log10(data.variables['aradii_m6'][t,:] / 1e-18) / dar

        plt.step(rr, nr, color='magenta', label="rain",   lw = 2.5, rasterized=False)
        plt.step(cr, nc, color='blue',    label="cloud",  lw = 2.5, rasterized=False)
        plt.step(ar, na, color='orange',  label="aerosol",lw = 2.5, rasterized=False)

        plt.axvline(x=1, color='gray')
        plt.axvline(x=25, color='gray')

        plt.legend(frameon=False)

        plt.savefig(output_folder + 'acnv_size_distr_m6_' + str("%03d" % t) + '.pdf')


def main():

    outfile = "test_acnv.nc"
    aerosol = '{"ammonium_sulfate": {"kappa": 0.61, "mean_r": [0.08e-6], "gstdev": [2.4], "n_tot": [500.0e6]}}'
    out_bin = '{"rradii": {"rght": 1000e-6, "left": 25e-6,  "drwt": "wet", "lnli": "log", "nbin": 51, "moms": [0]},\
                "cradii": {"rght": 25e-6,   "left": 1e-6,   "drwt": "wet", "lnli": "log", "nbin": 26, "moms": [0]},\
                "aradii": {"rght": 1e-6,    "left": 1e-9,   "drwt": "wet", "lnli": "log", "nbin": 26, "moms": [0]}}'


    # run parcel run!
    parcel(dt = .5, w = .5, sd_conc = 8192, outfreq = 400, z_max = 2000., RH_0 = 0.9, T_0 = 300, p_0 = 90000.,\
           outfile = outfile, out_bin = out_bin, aerosol = aerosol, wait = 0)

    data = netcdf.netcdf_file(outfile, "r")

    # plotting
    plot_spectrum(data, output_folder="../outputs/")

    # cleanup
    #subprocess.call(["rm", outfile])

if __name__ == '__main__':
    main()

