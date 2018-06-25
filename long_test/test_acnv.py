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
from autoconversion import plot_spectrum

import functions as fn

@pytest.fixture(scope="module")
def data(request):

    p_dict = {}
    p_dict['outfile']  = "test_acnv.nc"
    p_dict['sd_conc']  = 10000
    p_dict['outfreq']  = 200
    p_dict['dt']       = 0.5
    p_dict['w']        = 1.
    p_dict['z_max']    = 1500.
    p_dict['RH_0']     = 0.9
    p_dict['T_0']      = 300
    p_dict['p_0']      = 90000.
    p_dict['wait']     = 0
    p_dict['coal']     = True
    p_dict['coal_kernel'] = "hall"
    p_dict['terminal_vel'] = "khvorostyanov_nonspherical"

    p_dict['aerosol'] = '{"ammonium_sulfate": {"kappa": 0.61, "mean_r": [0.08e-6], "gstdev": [2.4], "n_tot": [500.0e6]}}'
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

def test_chem_plot(data):
    """
    quicklook for spectrum
    """
    plot_spectrum(data, output_folder="plots/outputs/")
