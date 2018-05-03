#!/bin/bash

rm -f ../script_wait.sh
printf "%s\n" "#!/bin/bash" >> ../script_wait.sh
CODIR=$(pwd)

for gamma in $(seq -w 10 2 12)
   do

   for beta in $(seq -w 2.0 0.5 3.0)
      do

      for alpha in $(seq -w 1 1 2)
         do

         basetime=$(date +%s%N)
         WODIR=$CODIR/../DATI/L_$gamma/D_$beta/$basetime

            if [  ! -d $WODIR ]
               then mkdir -p $WODIR
            fi

         cp $CODIR/lancio.sh $WODIR/lancio.tp

         for x in $(seq 1 1 8)
            do

            WODIR_core=$WODIR/core_$x

            if [  ! -d $WODIR_core ]
               then mkdir -p $WODIR_core
            fi

            cp DMRG_ST/dmrgtev-12.5/dmrgtev $WODIR_core

            IN=L_$gamma-D_$beta-$basetime-core_$uga.inp
            OUT=L_$gamma-D_$beta-$basetime-core_$uga.out

            awk -v Beta=$beta -v Gamma=$gamma ' 
            {if (($1=="define") && ($2=="L"))       $3="     " Gamma;       }
            {if (($1=="define") && ($2=="Delta"))   $3="     " Beta;        } { print $0;} 
            '  rand.inp  > $WODIR_core/$IN;

                                 
            # crea e metti il lancio.sh & dmrg
            ############################################################

            printf "

            #.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-

            WORKDIR=\"$WODIR_core\"
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

