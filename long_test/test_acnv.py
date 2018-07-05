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
from autoconversion import plot_spectrum_m0

import functions as fn

@pytest.fixture(scope="module")
def data(request):

    p_dict = {}
    p_dict['outfile']  = "test_acnv.nc"
    p_dict['sd_conc']  = 1000
    p_dict['outfreq']  = 300
    p_dict['dt']       = 1.
    p_dict['w']        = 0.5
    p_dict['z_max']    = 1500.
    p_dict['RH_0']     = 0.9
    p_dict['T_0']      = 10
    p_dict['p_0']      = 90000.
    p_dict['wait']     = 0

    p_dict['coal']     = True
    p_dict['coal_kernel'] = "onishi_hall"
    p_dict['terminal_vel'] = "khvorostyanov_spherical"
    p_dict['coal_dissipation_rate'] = 0.01  # or 0.04
    p_dict['coal_Reynolds_number'] = 100.

    # initial aerosol: 2-mode ammonium sulfate
    p_dict['aerosol'] = '{"ammonium_sulfate": {"kappa":  0.61,\
                                                 "mean_r": [0.02e-6,  0.07e-7],\
                                                 "gstdev": [1.4,      1.2],\
                                                 "n_tot":  [120.0e6,  80.0e6]}}'

    # additional giant aerosols
    p_dict['gccn_flag'] = False
    p_dict['large_tail'] = True
    p_dict['gccn'] = '{"kappa" : 1.28, "r_d" : [10e-7, 80e-7], "n_tot" : [10, 5]}'

    p_dict['out_bin'] = '{\
                "rradii": {"rght": 1000e-6, "left": 25e-6,  "drwt": "wet", "lnli": "log", "nbin": 51, "moms": [0]},\
                "cradii": {"rght": 25e-6,   "left": 1e-6,   "drwt": "wet", "lnli": "log", "nbin": 26, "moms": [0]},\
                "aradii": {"rght": 1e-6,    "left": 1e-9,   "drwt": "wet", "lnli": "log", "nbin": 26, "moms": [0]}}'

    #pp.pprint(p_dict)

    # run parcel
    parcel(**p_dict)

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
