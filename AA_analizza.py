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

directory = sorted(glob.glob('../DATI*/*/*/'))	
namesLD = ff.folder_crea(directory)


comm = MPI.COMM_WORLD
size = comm.Get_size()
rank = comm.Get_rank()

t_max=10002


for kk,name in enumerate(namesLD):
	if kk%size!=rank: continue
	
	L = name[0]
	D = name[1]

	L_int = int(L)+1

	dir_nameALL = '../**/L_'+L+'/D_'+D+'/*/corr.prp'
	All_files = glob.glob(dir_nameALL, recursive=True)
	
	n_rel=len(All_files)

	Big_Mat = np.zeros((n_rel,t_max,L_int), dtype=np.float)

	i=0

	for filename in All_files:

		dir_path  = os.path.dirname(filename)
		corrcon_file = dir_path+'/corr_con.prp'

		if os.path.isfile(corrcon_file) :
			DD 		  = pd.read_csv(corrcon_file, header=None, sep=r"\s+")
			corr_conn = pd.DataFrame.as_matrix(DD)

		else:
			corr_conn = ff.CorCon_exp(filename)

		corr_conn_nan = ff.putnan(t_max,corr_conn)

		Big_Mat[i] = corr_conn_nan
		i+=1

	dirAV_path = '../average/L_'+L+'/D_'+D
	
	media_0 = np.nanmean(Big_Mat,axis=0)
	media   = np.hstack((media_0,np.array([media_0[:,0]]).T))
	np.savetxt(dirAV_path+'/corr_aver.prp', media, fmt='%.9f')

	std_0   = np.nanstd(Big_Mat,axis=0)
	std     = np.hstack((std_0,np.array([std_0[:,0]]).T))
	np.savetxt(dirAV_path+'/corr_std.prp', std, fmt='%.9f')
	
	print(L, D, n_rel, "FATTO")



