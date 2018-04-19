import csv
import numpy as np
import math
import scipy.linalg as _la
import scipy.special as special
from math import factorial 
from datetime import datetime
import time
import re
import os as os
import glob
from numpy import inf
import re
import fnmatch
import pandas as pd
import shutil

#.............................................Export connected function
def	Io_fascio_tuto(name,LOCAL):

	L = name[0]
	D = name[1]

	t_max=10002
	L_int = int(L)

	dir_nameALL = LOCAL+'/**/L_'+L+'/D_'+D+'/*/corr.prp'
	All_files = glob.glob(dir_nameALL, recursive=True)
	
	n_rel=len(All_files)

	Big_Mat = np.zeros((n_rel,t_max,L_int+1), dtype=np.float)
	#print(L_int, 'lung')

	i=0

	for filename in All_files:

		dir_path  = os.path.dirname(filename)
		corrcon_file = dir_path+'/corr_con.prp'
		print(corrcon_file, (i+1), '/', n_rel)

	#	print(filename,L_int)

		if os.path.isfile(corrcon_file) :
			DD 		  = pd.read_csv(corrcon_file, header=None, sep=r"\s+")
			corr_conn = pd.DataFrame.as_matrix(DD)

		else:
			corr_conn = CorCon_exp(L_int,LOCAL,filename)

		corr_conn_nan = putnan(L_int+1,t_max,corr_conn)

	#	print(corr_conn.shape)

		Big_Mat[i] = corr_conn_nan
		i+=1


	dirAV_path = LOCAL+'/average/L_'+L+'/D_'+D
	
	media_0 = np.nanmean(Big_Mat,axis=0)
	media   = np.hstack((media_0,np.array([media_0[:,0]]).T))
	np.savetxt(dirAV_path+'/corr_aver.prp', media, fmt='%.9f')

	std_0   = np.nanstd(Big_Mat,axis=0)
	std     = np.hstack((std_0,np.array([std_0[:,0]]).T))
	np.savetxt(dirAV_path+'/corr_std.prp', std, fmt='%.9f')

	print(L, D, n_rel, "FATTO")	

	return 1

#.............................................Export connected function
def	CorCon_exp(Lint,LOCAL,filename):
	
	str_spl_tmp = re.split('L',filename)
	str_spl     = re.split('_|/',str_spl_tmp[1])

	L = str_spl[1]
	D = str_spl[3]

	dir_path  = os.path.dirname(filename)
	corr_file = dir_path+'/corr.prp'
	dens_file = dir_path+'/dens.prp'
	corrcon_file = dir_path+'/corr_con.prp'

	CC   = pd.read_csv(corr_file, header=None, sep=r"\s+")
	corr = pd.DataFrame.as_matrix(CC)

	DD   = pd.read_csv(dens_file, header=None, sep=r"\s+")
	dens = pd.DataFrame.as_matrix(DD)


	#Lint = int(np.amax(dens[:,1]))+1

	lt_corr = int(np.shape(corr)[0]/(Lint*Lint))-2
	lt_dens = int(np.shape(dens)[0]/(Lint))-2
	t_max = int(np.amin([lt_corr,lt_dens]))-1
	
	corr_tab  = np.reshape(corr[:t_max*Lint*Lint,2],(lt_corr,Lint,Lint))
	dens_tab  = np.reshape(dens[:t_max*Lint     ,2],(lt_dens,Lint))

	corr_aver = np.zeros((t_max,Lint+1))

	for x in range(0,t_max):
		dens_dens = np.tensordot(dens_tab[x],dens_tab[x],0)
		data_tab0 = corr_tab[x]-dens_dens
		data_tab  = Trasl_Mean(data_tab0)
		plus_last_el = np.append(data_tab,data_tab[0])
		corr_aver[x] = plus_last_el

	namefold = '/L_'+L+'/D_'+D+'/corr_con'
	dat_fold = LOCAL+'/datas'
	dirdat   = os.path.abspath(glob.glob(dat_fold)[0])
	
	namegen  = generate_filename(dirdat+os.sep+namefold)

	np.savetxt(namegen, corr_aver, fmt='%.9f')
	np.savetxt(corrcon_file, corr_aver, fmt='%.9f')
	print('bella YO')
	return corr_aver


#......................................................Traslations MEAN
def Trasl_Mean(A):
	a = A.shape
	B = np.zeros((a[1],a[1]), dtype=np.float)
	for i in range(a[1]):
		B[i] = np.roll(A[i],-i)
	return np.mean(B, axis=0)


#......................................................put nan
def putnan(Space,t,A):
	Time  = A.shape[0]
	Space = A.shape[1]	
	nantime = t-Time

#	print(Space, 'putnan')

	B = np.empty((nantime,Space,))
	B[:] = np.nan

#	print(A.shape,B.shape)

	xx = np.concatenate((A,B), axis=0)
	return xx

#....................................................creation of folder
def	folder_crea(LOCAL,directory):	
	names = [[0 for x in range(2)] for y in range(len(directory))]
	j=0
	dat_fold = LOCAL+'/datas'
	ave_fold = LOCAL+'/average'

	print(dat_fold, ave_fold)

	if os.path.exists(dat_fold):
		shutil.rmtree(dat_fold, ignore_errors=True)

	if os.path.exists(ave_fold):
		shutil.rmtree(ave_fold, ignore_errors=True)

	for dir_name in directory:

		str_spl = re.split('/|_',dir_name)
		D=str_spl[-2]
		L=str_spl[-4]

		names[j][0] = L
		names[j][1] = D
		dirDT_path = dat_fold+'/L_'+L+'/D_'+D

		if not os.path.exists(dirDT_path):
			os.makedirs(dirDT_path)

		dirAV_path = ave_fold+'/L_'+L+'/D_'+D

		if not os.path.exists(dirAV_path):
			os.makedirs(dirAV_path)


		j+=1

	xx = list(set(tuple(element) for element in names))
	return xx


#......................................................FileName_abstime
def generate_filename(basename):
	unix_timestamp = int(time.time())
	local_time = str(int(round(time.time() * 1000)))
	xx = basename + local_time + ".dat"
	if os.path.isfile(xx):
		time.sleep(1)
		return generate_filename(basename)
	return xx		

#..............................................................FindName
def find_files(treeroot, pattern):
	results = []
	for base, dirs, files in os.walk(treeroot):
		goodfiles = fnmatch.filter(files, pattern)
		results.extend(os.path.join(base, f) for f in goodfiles)
	return results


