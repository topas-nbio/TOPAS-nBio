#!/bin/tcsh

set ITER = $1
set ITER2 = $2
set ADDITION = $3
set DOSE = $4
set NPULSES = $5
set TPULSES = $6
set MEAN = $7
set FWHM = $8

if ($ITER == "") then
  set ITER = 1
endif

set OPTION = `cat inputfiles.txt`
foreach LINE ( $OPTION )
  set COUNT = $ITER 
  while ($COUNT < $ITER2)
    echo $LINE $COUNT
    set INFILE = `basename $LINE .py`
    set USER   = `whoami`
    set CURRENTPATH = `pwd`
    set DATEDAY   = `date | awk '{print $3}'`
    set DATEMONTH = `date | awk '{print $2}'`
    set DATEYEAR  = `date | awk '{print $6}'`
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
    cp ParameterFiles/TOPAS*tps $DIR
    cp ParameterFiles/pUC* $DIR
    cp ParameterFiles/Plasmid_50ugg_Sphere_1um_diameter_envelopes.xyz $DIR
    cp ParameterFiles/TOPAS*Reactions.txt $DIR

    cp $LINE $DIR
   
    set SEED = `bash -c 'echo $RANDOM'` 
    echo Ts/Seed = $SEED >> $DIR/depFile1.txt
    echo d:Sc/nbOfMol/PrescribedDose    = $DOSE Gy >> $DIR/depFile1.txt
    echo i:Sc/nbOfMol/NumberOfPulses    = $NPULSES  >> $DIR/depFile1.txt
    echo d:Sc/nbOfMol/TimeUpper                 = $TPULSES s  >> $DIR/depFile1.txt
    echo d:Sc/nbOfMol/PulseSimulationTime       = $TPULSES s  >> $DIR/depFile1.txt 
    echo d:Sc/nbOfMol/PulseTimeMean     = $MEAN s >> $DIR/depFile1.txt 
    echo d:Sc/nbOfMol/PulseTimeFWHM     = $FWHM s >> $DIR/depFile1.txt 

    set SCRIPT=$DIR/$INFILE-$COUNT".csh"

    cat - << EOF > $SCRIPT

#!/bin/bash
cd $DIR
nohup time python3 $INFILE.py nBio_dev > log.out &
EOF
    chmod +x $SCRIPT
    bash $SCRIPT 
    @ COUNT = $COUNT + 1
  end
  #python3 analysis/analysis.py $CURRENTPATH"/run/"$ADDITION$DATE/$INFILE $CURRENTPATH"/run/"$ADDITION$DATE/$INFILE
end
