#from mpi4py import MPI
import numpy as np
import glob
import pandas as pd
import re
import time
import os as os
import re
import AA_functions as ff


LOCAL = os.path.abspath('../')

#....................dopo avere creato cartelle e sotto cartelle
#....................restituisce in un vettori le varie possibili combinazioni L,D

path_file = LOCAL+'/DATI*/*/*/'
directory = sorted(glob.glob(path_file))	
namesLD = ff.folder_crea(LOCAL,directory)


for kk,name in enumerate(namesLD):
	
	ff.Io_fascio_tuto(name,LOCAL)
	uga = 1

	print('YO', name)



