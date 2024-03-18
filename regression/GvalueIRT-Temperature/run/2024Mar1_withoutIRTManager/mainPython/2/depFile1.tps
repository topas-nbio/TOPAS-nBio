# Setup to help to reproduce IRT temperature-dependent G-values

################# TOPAS PARAM #################

i:Ts/ShowHistoryCountAtInterval = 100
i:Ts/MaxInterruptedHistories  = So/Demo/NumberOfHistoriesInRun
Ts/NumberOfThreads = 0
b:Ts/ShowCPUTime = "True"

################## MODEL LISTS ##################

sv:Ph/Default/Modules = 2 "TsEmDNAPhysics" "TsEmDNAChemistry"
s:Ph/Default/Electron/SetElasticScatteringModel  = "ELSEPA" 

#Solvated Electron options are: ritchie, terrisol, meesungnoen, meesungnoensolid, krepl
s:Ph/Default/SolvatedElectronThermalizationModel = "Meesungnoen" 

################### CHEMISTRY ###################

s:Ch/ChemistryName = "TOPASChemistry"
b:Ch/TOPASChemistry/ChemicalStageTransportActive = "True"
dv:Ch/TOPASChemistry/ChemicalStageTimeStepsHighEdges    = 1 999999 ps
dv:Ch/TOPASChemistry/ChemicalStageTimeStepsResolutions  = 1 1e-6 ps
d:Ch/TOPASChemistry/ChemicalStageTimeEnd    = 1.00001 ps
u:Ch/TOPASChemistry/Temperature             = fTemperature
b:Ch/TOPASChemistry/TestForContactReactions = "True"
b:Ch/TOPASChemistry/ApplyCorrectionScalingForTemperature = "True"

includeFile = depFile2.tps

################# Water Density #################

s:Ma/G4_WATER_MODIFIED/CloneFromMaterial = "G4_WATER"
d:Ma/G4_WATER_MODIFIED/CloneWithDensity  = fDensity g/cm3

#################### GEOMETRY ###################

d:Ge/World/HLX= 10 cm
d:Ge/World/HLY= 10 cm
d:Ge/World/HLZ= 10 cm
s:Ge/World/Material ="G4_WATER"

d:Ge/Target/HLX= Ge/World/HLX cm 
d:Ge/Target/HLY= Ge/World/HLY cm
d:Ge/Target/HLZ= Ge/World/HLZ cm 
s:Ge/Target/Material = "G4_WATER"
s:Ge/Target/Type     = "TsBox"
s:Ge/Target/Parent   = "World"

################# SOURCE #################

Ge/BeamPosition/TransZ = 0.0 um
Ge/BeamPosition/RotX = 0 deg

So/Demo/BeamPositionDistribution = "None"
So/Demo/BeamEnergySpread = 0
So/Demo/BeamParticle = "e-"
So/Demo/BeamEnergy   = 999.999 keV
i:So/Demo/NumberOfHistoriesInRun = 1000

################# SCORER #################

u:perCent = 0.01 

s:Sc/nbOfMol/Quantity   = "TsIRTGvalue"
s:Sc/nbOfMol/Component  = "Target"
s:Sc/nbOfMol/OutputType = "ascii"
i:Sc/nbOfMol/OutputBufferSize = 1
d:Sc/nbOfMol/KillPrimaryIfEnergyLossExceeds = So/Demo/BeamEnergy keV * perCent
s:Sc/nbOfMol/IfOutputFileAlreadyExists      = "Overwrite"
b:Sc/nbOfMol/KillPrimaryBasedOnEnergyLoss   = "true"
d:Sc/nbOfMol/KillPrimaryBasedOnTrackLength  = 100 km 
d:Sc/nbOfMol/AbortEventIfPrimaryEnergyLossExceeds = 1.001 * Sc/nbOfMol/KillPrimaryIfEnergyLossExceeds keV
d:Sc/nbOfMol/CutOffTime = 1 ps
d:Sc/nbOfMol/TimeLower  = 1.0 ps
d:Sc/nbOfMol/TimeUpper  = 1e7 ps 
i:Sc/nbOfMol/TimeBins   = 100
s:Sc/nbOfMol/OutputFile = "electron_Gvalue_Corrected_fTemperature"

##########################################
Ts/Seed = 6046
