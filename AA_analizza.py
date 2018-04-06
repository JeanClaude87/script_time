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


directory = sorted(glob.glob('../DATI1/*/*/'))

t_max  = 10002
t_step = 0.1


names = ff.folder_crea(directory)

#print(enumerate(names))

comm = MPI.COMM_WORLD
size = comm.Get_size()
rank = comm.Get_rank()


for kk,cose in enumerate(names):
	if kk%size!=rank: continue
	
	L = cose[0]
	D = cose[1]

	L_int  = int(L)

	dir_nameALL = '../**/L_'+L+'/D_'+D+'/*/corr.prp'
	All_files = glob.glob(dir_nameALL, recursive=True)
	
	n_rel=len(All_files)
	print(L, D, n_rel)

	Big_Mat = np.zeros((n_rel,t_max,L_int), dtype=np.float)

	i=0

	for filename in All_files:

		dir_path  = os.path.dirname(filename)
		corrcon_file = dir_path+'/corr_con.prp'

		if os.path.isfile(corrcon_file) :
			corr_conn = np.loadtxt(corrcon_file, dtype=np.float)
			print(corr_conn.shape)
		else:
			corr_conn = ff.CorCon_exp(filename,L_int,t_max)
			print(filename)

		Big_Mat[i] = corr_conn
		i+=1

	dirAV_path = '../average/L_'+L+'/D_'+D
	
	media_0 = np.nanmean(Big_Mat,axis=0)
	media   = np.hstack((media_0,np.array([media_0[:,0]]).T))
	np.savetxt(dirAV_path+'/corr_aver.prp', media, fmt='%.9f')

	std_0   = np.nanstd(Big_Mat,axis=0)
	std     = np.hstack((std_0,np.array([std_0[:,0]]).T))
	np.savetxt(dirAV_path+'/corr_std.prp', std, fmt='%.9f')
	print(L, D, n_rel, "FATTO")
	



















