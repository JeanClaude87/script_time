#!/bin/bash

rm -f script_wait.sh


printf "%s\n" "#!/bin/bash" >> script_wait.sh

for gamma in $(seq 10 2 16)
        do

for beta in $(seq 2.0 0.5 5.0)
        do

for alpha in $(seq 1 1 100)
        do

for uga in $(seq 1 1 8)
        do

basetime=$(date +%s%N)



        if [  ! -d DATI2/L_$gamma/D_$beta/rel_$alpha ]
                then mkdir -p DATI2/L_$gamma/D_$beta/$basetime
        fi

                IN=L_$gamma-D_$beta-$basetime.inp
                OUT=L_$gamma-D_$beta-$basetime.out

                CODIR=$(pwd)
                WODIR=$CODIR/DATI2/L_$gamma/D_$beta/$basetime

                        awk -v Beta=$beta -v Gamma=$gamma ' 
                        {if (($1=="define") && ($2=="L"))       $3="     " Gamma;       }
                        {if (($1=="define") && ($2=="Delta"))   $3="     " Beta;        } { print $0;} 
                        '  rand.inp  > $WODIR/$IN;

                        
# crea e metti il lancio.sh & dmrg
############################################################

cp DMRG_ST/dmrgtev $WODIR
cp lancio.sh $WODIR/lancio.tp

printf "

WORKDIR=\"$WODIR\"
CURRENTDIR=\"$CODIR\"

cd \${WORKDIR}

nohup ./dmrgtev < $IN > $OUT 

cd \${CODIR}

" >> $WODIR/lancio.tp

nome=L_$gamma-dis_$beta-rel_$alpha

sed -e "s/nomino/$nome/g" < $WODIR/lancio.tp > $WODIR/lancio.sh
rm $WODIR/lancio.tp

done

# scrivi i lancio su script_wait.sh
############################################################
printf "%s\n" "

cd $WODIR 
qsub lancio.sh
cd $CODIR

" >> script_wait.sh

done
done
done

wait
chmod +x script_wait.sh
#./script_wait.sh

