#!/bin/bash

### nom du job (a changer)
#$ -N nomino

###$ -q "E5-2670deb128A,E5-2670deb128B,E5-2670deb128C,E5-2670deb128D,E5-2670deb128E,E5-2670deb128F,E5-2667v4deb128nl,E5-2667v4deb256A,E5-2667v2d2deb128,E5-2667v2h6deb128"

### exporter les variables d'environnement sur tous les noeuds d'execution
#$ -V

###$ -e /dev/null
###$ -o /dev/null