from mpi4py import MPI
import numpy as np
import glob
import pandas as pd
import re
import time
import os as os
import re
import AA_functions as ff


LOCAL = os.path.abspath('.')

#....................dopo avere creato cartelle e sotto cartelle
#....................restituisce in un vettori le varie possibili combinazioni L,D

namesLD = ff.folder_crea()

dir_nameALL = '../**/corr.prp'
All_files   = glob.glob(dir_nameALL, recursive=True)

comm = MPI.COMM_WORLD
size = comm.Get_size()
rank = comm.Get_rank()


for kk,name in enumerate(All_files):
	if kk%size!=rank: continue

	ff.CorCon_exp(name)


