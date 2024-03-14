#!/bin/tcsh

module load topas/3.2.test
module load geant4/10.5.1

set ITER = $1
set ADDITION = $2
if ($ITER == "") then
  set ITER = 1
endif

set OPTION = `cat inputfiles.txt`
foreach LINE ( $OPTION )
  set COUNT = 0
  while ($COUNT < $ITER)
    echo $LINE $COUNT
    set INFILE = `basename $LINE .txt`
    set USER = `whoami`
    set CURRENTPATH = `pwd`
    set DATEDAY  = `date | awk '{print $3}'`
    set DATEMONTH = `date | awk '{print $2}'`
    set DATEYEAR = `date | awk '{print $6}'`
    set DATEHOUR  = `date | awk '{print $4}' | awk -F: '{print $1}'`
    set DATEMIN   = `date | awk '{print $4}' | awk -F: '{print $2}'`
    set DATE = $DATEYEAR$DATEMONTH$DATEDAY
    set UNAME = `uname`

    set DIR = $CURRENTPATH"/run/"$ADDITION$DATE/$INFILE/$COUNT
    if ( -d $DIR ) then
       echo Directory exists, removing and recreating $DIR
       rm -rf $DIR
    endif

    mkdir -p $DIR
    cp ParameterFiles/depFile*txt $DIR
    cp $LINE $DIR
   
    set SEED = `bash -c 'echo $RANDOM'` 
    echo Ts/Seed = $SEED >> $DIR/$INFILE.txt    


    set SCRIPT=$DIR/$INFILE-$COUNT".csh"

    cat - << EOF > $SCRIPT

#!/bin/bash
#BSUB -J dbscan
#BSUB -q long
#BSUB -r
#BSUB -C 0
#BSUB -n 1
#BSUB -R "rusage[mem=2500]"
#BSUB -Q "140"
cd $DIR

time topas $INFILE.txt 

EOF

    chmod +x $SCRIPT

    bsub  -e $DIR/log.err -o $DIR/log.out  < $SCRIPT
    
    @ COUNT = $COUNT + 1
  end
end

