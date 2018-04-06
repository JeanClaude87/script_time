import numpy as np
import math
import scipy.linalg as _la
import scipy.special as special
from math import factorial
import time 
from datetime import datetime
import time
import re
import os as os
from glob import glob
from numpy import inf
import re
import fnmatch


#.............................................Traslations MEAN
def Trasl_Mean(A):
	a = A.shape
	B = np.zeros((a[1],a[1]), dtype=np.float)
	for i in range(a[1]):
		B[i] = np.roll(A[i],-i)
	return np.mean(B, axis=0)

#.................Export connected function
def	CorCon_exp(filename,L_int,t_max):

	dir_path  = os.path.dirname(filename)
	corr_file = dir_path+'/corr.prp'
	dens_file = dir_path+'/dens.prp'
	corrcon_file = dir_path+'/corr_con.prp'

	#...........find Size
	corr = np.loadtxt(corr_file, dtype=np.float)	
	dens = np.loadtxt(dens_file, dtype=np.float)

	corr_conn  = CorCon_TrAver(corr,dens,L_int,t_max)

	np.savetxt(corrcon_file, corr_conn, fmt='%.9f')

	return corr_conn

#.................Average over Trasl of the connected function
def CorCon_TrAver(A,B,lx,t):
	#...A->corr
	#...B->dens	
	ltA = int(np.shape(A)[0]/(lx*lx))
	ltB = int(np.shape(B)[0]/(lx))

	corr_tab  = np.reshape(A[:,2],(ltA,lx,lx))
	dens_tab  = np.reshape(B[:,2],(ltB,lx))

	corr_aver = np.empty((int(t),lx))
	corr_aver[:] = np.nan

	for x in range(0,np.amin([ltA,ltB])):
		dens_dens = np.tensordot(dens_tab[x],dens_tab[x],0)
		data_tab  = corr_tab[x]-dens_dens
		corr_aver[x] = Trasl_Mean(data_tab)
	return corr_aver


#.............................................creation of folder
def	folder_crea(directory):
	names = [[0 for x in range(2)] for y in range(len(directory))]
	j=0
	for dir_name in directory:

		str_spl = re.split('/|_',dir_name)
		D=str_spl[-2]
		L=str_spl[-4]

		names[j][0] = L
		names[j][1] = D
		dirAV_path = '../average/L_'+L+'/D_'+D

		if not os.path.exists(dirAV_path):
			os.makedirs(dirAV_path)

		j+=1

	xx = list(set(tuple(element) for element in names))
	return xx


#.............................................FileName_abstime
def generate_filename(basename):
	unix_timestamp = int(time.time())
	local_time = str(int(round(time.time() * 1000)))
	xx = basename + local_time + ".dat"
	if os.path.isfile(xx):
		time.sleep(1)
		return generate_filename(basename)
	return xx		

#.....................................................FindName
def find_files(treeroot, pattern):
	results = []
	for base, dirs, files in os.walk(treeroot):
		goodfiles = fnmatch.filter(files, pattern)
		results.extend(os.path.join(base, f) for f in goodfiles)
	return results


