#__init__.py
#GENERAL
from numpy import sqrt
import numpy as np
import ast
from numpy import unique
import pandas as pd
from scipy.spatial.distance import pdist,cdist
import os
import glob
import math as m
import itertools

# ASE
import ase.io.vasp
import ase
from ase.visualize import view
from ase import atoms,Atoms
from ase.calculators.emt import EMT
from ase.constraints import FixAtoms
from ase.optimize import QuasiNewton
from ase.build import bulk,fcc111, add_adsorbate

from ase.spacegroup.symmetrize import *
from ase.spacegroup.xtal import crystal
from ase.spacegroup import get_basis
from ase.spacegroup import get_spacegroup
from ase.spacegroup import Spacegroup

# HETEROSTRCTURES
import supercell_core as sc

# SLABCUT
from catkit.gen.surface import SlabGenerator

#PLOTTING
from matplotlib.cm import get_cmap
import matplotlib.pyplot as plt
from matplotlib.ticker import ScalarFormatter,AutoMinorLocator
import matplotlib as mpl

from matplotlib import rc
rc('font',**{'family':'sans-serif','sans-serif':['Helvetica']})
rc('font',**{'family':'serif','serif':['Garamond']})
rc('text', usetex=True)


import unittest
#from catkit.gen.surface import align_crystal
