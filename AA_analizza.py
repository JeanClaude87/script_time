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


t_max=10002

#comm = MPI.COMM_WORLD
#size = comm.Get_size()
#rank = comm.Get_rank()


for kk,name in enumerate(namesLD):
#	if kk%size!=rank: continue
	
	L = name[0]
	D = name[1]

	L_int = int(L)

	dir_nameALL = LOCAL+'/**/L_'+L+'/D_'+D+'/*/corr.prp'
	All_files = glob.glob(dir_nameALL, recursive=True)
	
	n_rel=len(All_files)

	Big_Mat = np.zeros((n_rel,t_max,L_int+1), dtype=np.float)
#	print(L_int, 'lung')

	i=0

	for filename in All_files:

		dir_path  = os.path.dirname(filename)
		corrcon_file = dir_path+'/corr_con.prp'

		print(filename,Lint)

		if os.path.isfile(corrcon_file) :
			DD 		  = pd.read_csv(corrcon_file, header=None, sep=r"\s+")
			corr_conn = pd.DataFrame.as_matrix(DD)

		else:
			corr_conn = ff.CorCon_exp(L_int,LOCAL,filename)

		corr_conn_nan = ff.putnan(L_int+1,t_max,corr_conn)

		Big_Mat[i] = corr_conn_nan
		i+=1

		print(corrcon_file)

	dirAV_path = LOCAL+'/average/L_'+L+'/D_'+D
	
	media_0 = np.nanmean(Big_Mat,axis=0)
	media   = np.hstack((media_0,np.array([media_0[:,0]]).T))
	np.savetxt(dirAV_path+'/corr_aver.prp', media, fmt='%.9f')

	std_0   = np.nanstd(Big_Mat,axis=0)
	std     = np.hstack((std_0,np.array([std_0[:,0]]).T))
	np.savetxt(dirAV_path+'/corr_std.prp', std, fmt='%.9f')
	
	print(L, D, n_rel, "FATTO")



