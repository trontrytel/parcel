import sys
import os
sys.path.insert(0, "/home/pracownicy/ajaruga/devel/libcloudphxx/build/bindings/python")

from libcloudphxx import common as cm

def sat_pressure(T):
    return cm.p_vs(T)

def epsilon():
    return cm.eps
