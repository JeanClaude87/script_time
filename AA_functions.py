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


#.............................................Export connected function
def	CorCon_exp(filename):

	
	str_spl_tmp = re.split('L',filename)
	str_spl     = re.split('_|/',str_spl_tmp[1])

	L = str_spl[1]
	D = str_spl[3]

	dir_path  = os.path.dirname(filename)
	corr_file = dir_path+'/corr.prp'
	dens_file = dir_path+'/dens.prp'
	corrcon_file = dir_path+'/corr_con.prp'

	print(corrcon_file)

	#if os.path.isfile(corrcon_file) :
	#	ciao=0
	#	print(corrcon_file, "FATTO")
	
	#else:

	CC   = pd.read_csv(corr_file, header=None, sep=r"\s+")
	corr = pd.DataFrame.as_matrix(CC)

	DD   = pd.read_csv(dens_file, header=None, sep=r"\s+")
	dens = pd.DataFrame.as_matrix(DD)


	Lint = int(np.amax(dens[:,1]))+1

	lt_corr = int(np.shape(corr)[0]/(Lint*Lint))
	lt_dens = int(np.shape(dens)[0]/(Lint))
	
	corr_tab  = np.reshape(corr[:,2],(lt_corr,Lint,Lint))
	dens_tab  = np.reshape(dens[:,2],(lt_dens,Lint))

	t_max = int(np.amin([lt_corr,lt_dens]))

	corr_aver = np.zeros((t_max,Lint+1))

	for x in range(0,t_max):
		dens_dens = np.tensordot(dens_tab[x],dens_tab[x],0)
		data_tab0 = corr_tab[x]-dens_dens
		data_tab  = Trasl_Mean(data_tab0)

		corr_aver[x] = np.append(data_tab,data_tab[0])

	print(corrcon_file)

	namefold = '/L_'+L+'/D_'+D+'/corr_con'
	dirdat   = os.path.abspath(glob.glob('../datas')[0])
	
	namegen  = generate_filename(dirdat+os.sep+namefold)

	np.savetxt(namegen, corr_aver, fmt='%.9f')

	return 1



#......................................................Traslations MEAN
def Trasl_Mean(A):
	a = A.shape
	B = np.zeros((a[1],a[1]), dtype=np.float)
	for i in range(a[1]):
		B[i] = np.roll(A[i],-i)
	return np.mean(B, axis=0)

#....................................................creation of folder
def	folder_crea(directory):	
	names = [[0 for x in range(2)] for y in range(len(directory))]
	j=0
	for dir_name in directory:

		str_spl = re.split('/|_',dir_name)
		D=str_spl[-2]
		L=str_spl[-4]

		names[j][0] = L
		names[j][1] = D
		dirDT_path = '../datas/L_'+L+'/D_'+D

		if not os.path.exists(dirDT_path):
			os.makedirs(dirDT_path)

		dirAV_path = '../average/L_'+L+'/D_'+D

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


