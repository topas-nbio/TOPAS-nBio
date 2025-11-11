#!/bin/csh
# TOPAS-nBio demo runner for C-shell users.
# This script locates the repository root so it can be invoked via
#   source rundemos.csh, ./rundemos.csh, or /path/to/rundemos.csh
# from any directory.

set start_dir = `pwd`
set repo_root = ""
set script_path = "$0"

# When sourced, $0 may just be the shell name. Fall back to env vars.
if ("$script_path" == "csh" || "$script_path" == "-csh" || "$script_path" == "tcsh") then
    set script_path = ""
endif

if ("$script_path" != "") then
    if ("$script_path" !~ */*) then
        set resolved = `which "$script_path" 2> /dev/null`
        if ($status == 0 && "$resolved" != "") then
            set script_path = "$resolved"
        endif
    endif
    set script_dir = "$script_path:h"
    if ("$script_dir" == "" || "$script_dir" == "$script_path") then
        set script_dir = "."
    endif
    if (-d "$script_dir") then
        pushd "$script_dir" > /dev/null
        set repo_root = `pwd`
        popd > /dev/null
    endif
endif

if ("$repo_root" == "" && $?TOPAS_NBIO_HOME) then
    set repo_root = "$TOPAS_NBIO_HOME"
endif

if ("$repo_root" == "") then
    set repo_root = "$start_dir"
endif

if (! -d "$repo_root/examples") then
    echo "Unable to locate TOPAS-nBio repository root from '$repo_root'."
    echo "Set TOPAS_NBIO_HOME or run the script via an explicit path."
    exit 1
endif

set _demo_had_nonomatch = $?nonomatch
if (! $_demo_had_nonomatch) then
    set nonomatch
endif

set geometry_cases = ( \
    dna/CharltonDNA.txt \
    dna/LinearDNA.txt \
    dna/CircularPlasmid.txt \
    cells/EllipsoidCell.txt \
    cells/FibroblastCell.txt \
    cells/SphericalCell.txt \
    cells/HiC/Parameters.txt \
    cells/Neuron/Neuron.txt \
    other/generateRandomCylinders.txt \
    other/readBackRandomCylindres.txt )

set process_cases = ( \
    ActiveChemistryDefault.txt \
    ActiveChemistryRevised.txt \
    ActiveCustomizablePhysics.txt \
    G4DNAModelPerRegion.txt \
    GoldNanoParticle.txt \
    RemoveChemicalSpeciesInVolume.txt )

set scorer_cases = ( \
    DBSCAN/DBSCAN.txt \
    Fricke/FrickeIRT.txt \
    IonizationDetail/IonizationDetailInRandomCylinders.txt \
    IRTGetGValue/TsIRTGvalue.txt \
    IRTInterTrack/TsIRTInterTrack.txt \
    IRTCummulative/TsIRTCummulative.txt \
    IRTTemperature/TemperatureExample_90C.txt \
    IRTAndGillespieContinuous/TsIRTAndGillespieGvalue.txt \
    ProteinDataBank/PDB4DNA.txt \
    ParticleTuple/particleTuple.txt \
    SBSDamageToDNAPlasmid/FullDNADamageInPlasmid.txt \
    SBSDamageToDNAPlasmid/SSBandDSbWithDBSCAN.txt \
    SBSMoleculesAtATime/TsSpeciesAtTime.txt \
    SBSGetGValue/GvalueG4DNADefault.txt \
    SBSGetGValue/GvalueRevisedPhysicsChemistry.txt )

set green = `tput setaf 2` 
set red   = `tput setaf 1`
set reset = `tput sgr0` 

# --- Detect UTF-8 support and choose symbols ---
if (`locale charmap` == "UTF-8") then
    set check = "✅"
    set cross = "❌"
else
    set check = "[ Success ]"
    set cross = "[ FAIL ]"
endif

set exec_log = out

set start_time = `date +%s` 
echo "Running TOPAS-nBio demos using repository root: $repo_root"

cd "$repo_root/examples/geometry"
echo ""
echo "======================================"
echo "** Geometry demos **"
echo "======================================"
foreach case ($geometry_cases)
    set case_dir = $case:h
    set case_file = $case:t
    if ("$case_dir" == "" || "$case_dir" == "$case") then
        set case_dir = "."
    endif
    pushd "$case_dir" > /dev/null
    echo "includeFile = $case_file" > run.txt
    echo 'b:Ts/UseQt = "False"' >> run.txt
    echo 'b:Ts/PauseBeforeQuit = "False"' >> run.txt
    echo 'b:Gr/Enable = "False"' >> run.txt
    echo "🚀 Running $case " 
    topas run.txt >&! $exec_log
    set matches = `awk '/Execution/ {count++} END {print count+0}' $exec_log`
    if ("$matches" == "1") then
        printf "   %s%s [ Success ]%s %s\n" "$green" "$check" "$reset" "$case" 
    else
        printf "   %s%s [ Fail ]%s %s\n" "$red" "$cross" "$reset" "$case" 
    endif
    rm -f $exec_log run.txt *.xyz *.phsp *.header *.csv *.sdd
    popd > /dev/null
end
rm -f *.xyz *.phsp *.header *.csv $exec_log

cd "$repo_root/examples/processes"
echo ""
echo "======================================"
echo "** Process demos **"
echo "======================================"
foreach case ($process_cases)
    set case_dir = $case:h
    set case_file = $case:t
    if ("$case_dir" == "" || "$case_dir" == "$case") then
        set case_dir = "."
    endif
    pushd "$case_dir" > /dev/null
    echo "includeFile = $case_file" > run.txt
    echo 'b:Ts/UseQt = "False"' >> run.txt
    echo 'b:Ts/PauseBeforeQuit = "False"' >> run.txt
    echo 'b:Gr/Enable = "False"' >> run.txt
    echo "🚀 Running $case " 
    topas run.txt >&! $exec_log
    set matches = `awk '/Execution/ {count++} END {print count+0}' $exec_log`
    if ("$matches" == "1") then
        printf "   %s%s [ Success ]%s %s\n" "$green" "$check" "$reset" "$case" 
    else
        printf "   %s%s [ Fail ]%s %s\n" "$red" "$cross" "$reset" "$case" 
    endif
    rm -f $exec_log run.txt
    popd > /dev/null
end
rm -f $exec_log

cd "$repo_root/examples/scorers"
echo ""
echo "======================================"
echo "** Scorer demos **"
echo "======================================"
foreach case ($scorer_cases)
    set case_dir = $case:h
    set case_file = $case:t
    if ("$case_dir" == "" || "$case_dir" == "$case") then
        set case_dir = "."
    endif
    pushd "$case_dir" > /dev/null
    echo "includeFile = $case_file" > run.txt
    echo 'b:Ts/UseQt = "False"' >> run.txt
    echo 'b:Ts/PauseBeforeQuit = "False"' >> run.txt
    echo 'b:Gr/Enable = "False"' >> run.txt
    echo "🚀 Running $case " 
    topas run.txt >&! $exec_log
    set matches = `awk '/Execution/ {count++} END {print count+0}' $exec_log`
    if ("$matches" == "1") then
        printf "   %s%s [ Success ]%s %s\n" "$green" "$check" "$reset" "$case" 
    else
        printf "   %s%s [ Fail ]%s %s\n" "$red" "$cross" "$reset" "$case" 
    endif
    rm -f $exec_log run.txt *.root *.phsp *.header *.csv *.bin
    popd > /dev/null
end
rm -f *.root *.phsp *.header *csv $exec_log

cd "$start_dir"

set end_time = `date +%s`
@ elapsed = $end_time - $start_time

printf "\n📊  Regression completed in %d seconds (%.2f minutes)\n" $elapsed `echo "$elapsed / 60.0" | bc -l`

if (! $_demo_had_nonomatch) then
    unset nonomatch
endif

