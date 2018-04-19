from mpi4py import MPI
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

comm = MPI.COMM_WORLD
size = comm.Get_size()
rank = comm.Get_rank()

if rank == 0:
	path_file = LOCAL+'/DATI*/*/*/'
	directory = sorted(glob.glob(path_file))	
	namesLD = ff.folder_crea(LOCAL,directory)


namesLD = comm.scatter(data, root=0)
print 'rank',rank,'has data:',data



#for kk,name in enumerate(namesLD):
#	if kk%size!=rank: continue
	
	#ff.Io_fascio_tuto(name,LOCAL)
#	uga = 1

print('YO', name)



