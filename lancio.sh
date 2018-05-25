#!/bin/bash

### nom du job (a changer)
#$ -N nomino

### parallel environment & nb cpu (NSLOTS)

##$ -pe mpi8_debian 8
###$ -q "E5-2667v4deb256A,E5-2670deb128B,E5-2670deb128C,E5-2670deb128D,E5-2670deb128E,E5-2670deb128F,E5-2670deb128A,E5-2670deb128nl"

#$ -pe openmp8 8
#$ -q "h6-E5-2667v4deb128"

module load icc/2017.4b

### exporter les variables d'environnement sur tous les noeuds d'execution
#$ -V

###$ -e /dev/null
###$ -o /dev/null
