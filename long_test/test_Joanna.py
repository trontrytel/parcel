import sys
sys.path.insert(0, "../")
sys.path.insert(0, "./")
sys.path.insert(0, "plots/one_simulat/")

from scipy.io import netcdf
import numpy as np
import math
import subprocess
import pytest
import copy
import pprint as pp

from parcel import parcel
from libcloudphxx import common as cm
from autoconversion import plot_spectrum_m0, plot_profiles, plot_spectrum_m3, plot_spectrum_m6

#from open_nc import open_nc

import functions as fn

@pytest.fixture(scope="module")
def data(request):

    p_dict = {}
    p_dict['outfile']  = "Joanna_test_acnv.nc"
    p_dict['sd_conc']  = 30000
    p_dict['outfreq']  = 100
    p_dict['dt']       = 1        # dt = 5 unstable   dt = 1 stable
    p_dict['w']        = 0.5
    p_dict['z_max']    = 1000.
    p_dict['RH_0']     = 0.99     # 0.9
    p_dict['T_0']      = 290      # 350K that is 170F. That is a too high to be reasonable atmospheric condition
    p_dict['p_0']      = 90000.
    p_dict['wait']     = 0
    p_dict['coal']     = True
    p_dict['coal_kernel'] = "hall_davis_no_waals"
    p_dict['terminal_vel'] = "khvorostyanov_spherical"
    p_dict['aerosol'] = '{"ammonium_sulfate": {"kappa": 0.61, "mean_r": [0.06e-6, 0.12e-6], "gstdev": [1.6, 2.5], "n_tot": [500.0e6, 320.0e6]}}'
    #p_dict['out_bin'] = '{\
    #                "rradii": {"rght": 1000e-6, "left": 25e-6,  "drwt": "wet", "lnli": "log", "nbin": 1, "moms": [6]},\
    #                "cradii": {"rght": 25e-6,   "left": 1e-6,   "drwt": "wet", "lnli": "log", "nbin": 1, "moms": [6]},\
    #                "aradii": {"rght": 1e-6,    "left": 1e-9,   "drwt": "wet", "lnli": "log", "nbin": 1, "moms": [6]}}'

    p_dict['out_bin'] = '{\
                "rradii": {"rght": 1000e-6, "left": 25e-6,  "drwt": "wet", "lnli": "log", "nbin": 51, "moms": [0, 3, 6]},\
                "cradii": {"rght": 25e-6,   "left": 1e-6,   "drwt": "wet", "lnli": "log", "nbin": 26, "moms": [0, 3, 6]},\
                "aradii": {"rght": 1e-6,    "left": 1e-9,   "drwt": "wet", "lnli": "log", "nbin": 26, "moms": [0, 3, 6]}}'
    #open_nc(p_dict)

    # run parcel
    #parcel(**p_dict)
    #simulation results
    data = netcdf.netcdf_file(p_dict['outfile'],   "r")
    # removing all netcdf files after all tests
    def removing_files():
        subprocess.call(["rm", p_dict['outfile']])
    #request.addfinalizer(removing_files)
    return data

def test_acnv_plot(data):
    """
    quicklook for spectrum
    """
    plot_spectrum_m0(data, output_folder="plots/outputs/")
    plot_spectrum_m3(data, output_folder="plots/outputs/")
    #plot_spectrum_m6(data, output_folder="plots/outputs/")

def test_profiles(data):
    """
    quicklook for profiles
    """
    plot_profiles(data, output_folder="plots/outputs/")
