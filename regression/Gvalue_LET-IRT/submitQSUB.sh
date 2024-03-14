#!/bin/tcsh

#module load topas/3.2.test
#module load geant4/10.5.1


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
    set INFILE = `basename $LINE .py`
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
    echo Ts/Seed = $SEED >> $DIR/depFile1.txt


    set SCRIPT=$DIR/$INFILE-$COUNT".csh"

    cat - << EOF > $SCRIPT

#!/bin/env bash
cd $DIR

export TOPAS_G4_DATA_DIR=/wynton/home/faddegonlab/jramos/software/G4Data
export LD_LIBRARY_PATH=/wynton/home/faddegonlab/jramos/software/topas.3.5/lib

/wynton/home/faddegonlab/jramos/software/topas.3.5/bin/topas $INFILE.txt  > log.out

EOF

    chmod +x $SCRIPT

    qsub -cwd -pe smp 1 -l mem_free=2G < $SCRIPT
    
    @ COUNT = $COUNT + 1
  end
end

