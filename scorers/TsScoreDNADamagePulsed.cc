// Scorer for TsIRTStrandBreaksPulsed
//
// ********************************************************************
// *                                                                  *
// * This file is part of the TOPAS-nBio extensions to the            *
// *   TOPAS Simulation Toolkit.                                      *
// * The TOPAS-nBio extensions are freely available under the license *
// *   agreement set forth at: https://topas-nbio.readthedocs.io/     *
// *                                                                  *
// ********************************************************************
//

#include "TsScoreDNADamagePulsed.hh"
#include "TsIRTManager.hh"
#include "TsIRTConfiguration.hh"
#include "TsIRTUtils.hh"

#include "G4ITTrackHolder.hh"
#include "G4EventManager.hh"
#include "G4RunManager.hh"
#include "G4SystemOfUnits.hh"
#include "G4DNAChemistryManager.hh"
#include "G4Scheduler.hh"
#include "G4MolecularConfiguration.hh"
#include "G4Molecule.hh"
#include "G4H2O.hh"
#include "G4Electron_aq.hh"
#include "G4UnitsTable.hh"
#include "Randomize.hh"
#include "TsGeometryManager.hh"
#include "TsIRTPlasmidSupercoiled.hh"

TsScoreDNADamagePulsed::TsScoreDNADamagePulsed(TsParameterManager* pM, TsMaterialManager* mM, TsGeometryManager* gM, TsScoringManager* scM, TsExtensionManager* eM,
                                                     G4String scorerName, G4String quantity, G4String outFileName, G4bool isSubScorer)
: TsVNtupleScorer(pM, mM, gM, scM, eM, scorerName, quantity, outFileName, isSubScorer),
fPm(pM), fEnergyDepositPerEvent(0), fEnergyDepositPerEventEverywhere(0), fName(scorerName), fOldEvent(-1)
{
    SetUnit("");

    fIRT       = new TsIRTManager(fPm, fName);
    fUtils     = fIRT->GetUtils();
    fStepTimes = fIRT->GetStepTimes();
    
    fNtuple->RegisterColumnD(&fGValue,      "GValue", "");
    fNtuple->RegisterColumnD(&fGValueError, "GValue Error", "");
    fNtuple->RegisterColumnD(&fNumberOfMolecules, "Total Number of molecules","");
    fNtuple->RegisterColumnD(&fNumberOfMoleculesError, "Error of Number of molecules","");
    fNtuple->RegisterColumnD(&fEnergy,      "EnergyDeposit", "");
    fNtuple->RegisterColumnD(&fTime,        "Time", "picosecond");
    fNtuple->RegisterColumnS(&fMoleculeName, "MoleculeName");
    
    fPrescribedDose = fPm->GetDoubleParameter(GetFullParmName("PrescribedDose"),"Dose");

    fTimeDistribution = fPm->GetStringParameter(GetFullParmName("PulseDistribution"));

    G4StrUtil::to_lower(fTimeDistribution);

    fNumberOfPulses  = 1;
    fPulsesTimeDelay = 0;
    if ( fPm->ParameterExists(GetFullParmName("NumberOfPulses")) ) {
        fNumberOfPulses = fPm->GetIntegerParameter(GetFullParmName("NumberOfPulses"));
        G4double pulsesFrequency = fPm->GetDoubleParameter(GetFullParmName("PulsesFrequency"),"perTime");
        if ( fNumberOfPulses > 1 && pulsesFrequency == 0 ) {
            G4cerr << "TOPAS is exiting due to an error in parameter" << G4endl;
            G4cerr << GetFullParmName("PulsesFrequency") << G4endl;
            G4cerr << "The number of pulses is bigger than 1, therefore, the frequency must be larger than zero" << G4endl;
        }
        fPulsesTimeDelay = 1./pulsesFrequency;
    }
    
    fLowLimitTo1ps = false;
    if ( fPm->ParameterExists(GetFullParmName("ForceLowTimeCutTo1ps")) )
        fLowLimitTo1ps = fPm->GetBooleanParameter(GetFullParmName("ForceLowTimeCutTo1ps"));
    
    if ( fTimeDistribution == "gaussian" ) {
        fTimeDistributionType = 1;
        fTimeMean = fPm->GetDoubleParameter(GetFullParmName("PulseTimeMean"),"Time");
        fTimeFWHM = fPm->GetDoubleParameter(GetFullParmName("PulseTimeFWHM"),"Time");
        fTimeStdv = fTimeFWHM/2.354820045;
    } else if ( fTimeDistribution == "uniform" ) {
        fTimeDistributionType = 2;
        fTimeMean = fPm->GetDoubleParameter(GetFullParmName("PulseTimeMean"),"Time");
        fTimeFWHM = fPm->GetDoubleParameter(GetFullParmName("PulseTimeFWHM"),"Time");
        fTimeStdv = fTimeFWHM/2.354820045;
    } else if ( fTimeDistribution == "none" ) {
        fTimeDistributionType = 0;
    } else if ( fTimeDistribution == "discrete" ) {
        fTimeDistributionType = 3;
        fTimeValues = fPm->GetDoubleVector(GetFullParmName("TimeValues"),"Time");
        fTimeWeights = fPm->GetUnitlessVector(GetFullParmName("TimeWeights"));
        fNbOfTimes = fPm->GetVectorLength(GetFullParmName("TimeValues"));
        G4double sum = 0.0;
        for ( int i = 0; i < fNbOfTimes; i++ )
            sum += fTimeWeights[i];
        
        if ( sum != 1.0 )
            for ( int i = 0; i < fNbOfTimes; i++ )
                fTimeWeights[i] /= sum;
        
        fTimeTops.push_back(fTimeWeights[0]);
        for ( int i = 1; i < fNbOfTimes; i++ )
            fTimeTops.push_back(fTimeTops[i-1] + fTimeWeights[i]);
        
    } else if ( fTimeDistribution == "exponential" ) {
        fTimeDistributionType = 4;
        fTimeMean = fPm->GetDoubleParameter(GetFullParmName("PulseTimeMean"),"Time");
    } else {
        G4cerr << "TOPAS is exiting due to an error in parameter" << G4endl;
        G4cerr << GetFullParmName("PulseDistribution") << G4endl;
        G4cerr << "Distribution: " << fTimeDistribution << " does not found" << G4endl;
        G4cerr << "Available distributions: Gaussian, Poisson and Uniform" << G4endl;
    }
    
    // To filter only physical interactions within a virtual TsBox.
    fTestIsInside = false;
    if ( fPm->ParameterExists(GetFullParmName("SensitiveVolumeName"))) {
        fTestIsInside = true;
        fSensitiveVolume = fPm->GetStringParameter(GetFullParmName("SensitiveVolumeName"));
    }

    fNumberOfRepetitions = 1;
    if (fPm->ParameterExists(GetFullParmName("NumberOfRepetitions"))) {
        fNumberOfRepetitions = fPm->GetIntegerParameter(GetFullParmName("NumberOfRepetitions"));
    }

    fNumberOfThreads = fPm->GetIntegerParameter("Ts/NumberOfThreads");
    if (fNumberOfThreads <= 0)
        fNumberOfThreads = G4Threading::G4GetNumberOfCores() + fNumberOfThreads;
    
    fTCut = 1.0 * ps;

    fSavePulseInformation = true;
    if (fPm->ParameterExists(GetFullParmName("SavePulseInformation"))) {
        fSavePulseInformation = fPm->GetBooleanParameter(GetFullParmName("SavePulseInformation"));
    }

    fRootFileName = fPm->GetStringParameter("Sc/RootFileName") + ".bin";
    remove(fRootFileName);
    
    fNbOfScoredEvents = 0;
    fNbOfScoredEventsEverywhere = 0;
    fNbOfIRTRuns = 0;

    fTotalDose = 0.0;
    fDosePerPulse = 0.0;
    fPulseTimeShift = 0.0; //fPulsesTimeDelay;
    
    //SampleShiftTime();
    fMinShiftTime = 1.0 * s;
    fCurrentPulse = 1;

    fPulseTimeLimits[fCurrentPulse].first  =  1E150;
    fPulseTimeLimits[fCurrentPulse].second = 1E-150;

    fSimulationFinished = false;
    fTimeSampled        = false;
    fDNAHasBeenInserted = false;

    fVolume = 0;

    fComponent = fPm->GetStringParameter(GetFullParmName("Component"));

    for ( size_t u = 0; u < fStepTimes.size(); u++ ) {
        fVEnergyDeposit.push_back(0.0);
        fAccumulatedEnergy.push_back(0.0);
        fAccumulatedEnergy2.push_back(0.0);
    }

    fOutputFile = "TOPAS";
    if ( fPm->ParameterExists(GetFullParmName("OutputFile"))) {
        fOutputFile = fPm->GetStringParameter(GetFullParmName("OutputFile"));
    }

    fPulseRecycle = false;
    if ( fPm->ParameterExists(GetFullParmName("RecyclePulses"))) {
        fPulseRecycle = fPm->GetBooleanParameter(GetFullParmName("RecyclePulses"));
    }

    fTimeSimulationForPulse = 1 * us;
    if ( fPm->ParameterExists(GetFullParmName("PulseSimulationTime"))) {
        fTimeSimulationForPulse = fPm->GetDoubleParameter(GetFullParmName("PulseSimulationTime"),"Time");

    }

    G4int nbBreakMolecules   = fPm->GetVectorLength(GetFullParmName("StrandBreakMolecules"));
    G4String* breakMolecules = fPm->GetStringVector(GetFullParmName("StrandBreakMolecules"));

    for(G4int i = 0; i < nbBreakMolecules; i++) {
        fStrandBreakMolecules.push_back(breakMolecules[i]);
    }

    fMass = 0;

    if (fPm->ParameterExists(GetFullParmName("DNAMoleculesNames"))) {
        G4int NbOfDNANames   = fPm->GetVectorLength(GetFullParmName("DNAMoleculesNames"));
        G4int NbOfDNAWeights = fPm->GetVectorLength(GetFullParmName("DNAMoleculesWeights"));

        G4String* DNAMoleculesNames   = fPm->GetStringVector(GetFullParmName("DNAMoleculesNames"));
        G4double* DNAMoleculesWeights = fPm->GetUnitlessVector(GetFullParmName("DNAMoleculesWeights"));

        if (NbOfDNANames != NbOfDNAWeights) {
            G4cout << "Please provide as many DNA molecules as Weights!!" << G4endl;
            exit(1);
        }

        for (G4int i = 0; i < NbOfDNANames; i++) {
            fDNAMoleculesToInsert.push_back(DNAMoleculesNames[i]);
            fDNAMoleculesWeights.push_back(DNAMoleculesWeights[i]);
        }
    }
    else {
        fDNAMoleculesToInsert.push_back("TsDeoxyribose^0");
        fDNAMoleculesWeights.push_back(1);
    }

    if (fDNAMoleculesToInsert.size() != fDNAMoleculesWeights.size()) {
        G4cout << "Please provide as many DNA molecules as Weights!!" << G4endl;
        exit(1);
    }

    GetDNAInformation();
}


TsScoreDNADamagePulsed::~TsScoreDNADamagePulsed()
{
    delete fIRT;
}


G4bool TsScoreDNADamagePulsed::ProcessHits(G4Step* aStep, G4TouchableHistory*)
{
    if (!fIsActive) {
        fSkippedWhileInactive++;
        return false;
    }

    if (!fSimulationFinished) {
        G4int TrckID  = aStep->GetTrack()->GetTrackID();

        if ( -1 < TrckID ) {
            ResolveSolid(aStep);
            G4double edep = aStep->GetTotalEnergyDeposit() ;
            if ( edep > 0 ) {
                edep *= aStep->GetPreStepPoint()->GetWeight();
                
                G4double density = aStep->GetPreStepPoint()->GetMaterial()->GetDensity();
                fEnergyDepositPerEventEverywhere += edep ;

                if ( fTestIsInside ) {
                    G4TouchableHistory* touchable = (G4TouchableHistory*)(aStep->GetPreStepPoint()->GetTouchable());
                    const G4String& volumeName = touchable->GetVolume()->GetName();
                    if( !G4StrUtil::contains(volumeName,fSensitiveVolume)) {
                        return false;
                    }
                }

                G4TouchableHistory* touchable = (G4TouchableHistory*)(aStep->GetPreStepPoint()->GetTouchable());
                const G4String& vName = touchable->GetVolume()->GetName();

                if (G4StrUtil::contains(vName,"deoxyribose") || G4StrUtil::contains(vName,"phosphate")) {
                    G4int PlasmidID;
                    G4int BaseID = touchable->GetVolume(0)->GetCopyNo();
                    G4int StrandID;

                    if (G4StrUtil::contains(vName,"water")) {PlasmidID = touchable->GetVolume(1)->GetCopyNo();}
                    else {PlasmidID = touchable->GetVolume(2)->GetCopyNo();}

                    if (G4StrUtil::contains(vName,"1")) {StrandID = 1;}
                    else {StrandID = 2;}

                    if (fDNAEnergyDeposit.find(PlasmidID) == fDNAEnergyDeposit.end()) {
                        fDNAEnergyDeposit[PlasmidID][BaseID][StrandID] = edep;
                    } else {
                        if (fDNAEnergyDeposit[PlasmidID].find(BaseID) == fDNAEnergyDeposit[PlasmidID].end()) {
                            fDNAEnergyDeposit[PlasmidID][BaseID][StrandID] = edep;
                        } else {
                            if (fDNAEnergyDeposit[PlasmidID][BaseID].find(StrandID) == fDNAEnergyDeposit[PlasmidID][BaseID].end()) {
                                fDNAEnergyDeposit[PlasmidID][BaseID][StrandID] = edep;
                            } else {
                                fDNAEnergyDeposit[PlasmidID][BaseID][StrandID] += edep;
                            }
                        }
                    }
                }

                //else if (!G4StrUtil::contains(vName,"VoxelStraight")) {
                if (vName == "Plasmids") {
                    fMass = (density * fSolid->GetCubicVolume())/kg;
                }

                if (fMass == 0) {fMass = fVolume * density;}
            
                G4double dose = 0;

                if (fMass > 0) {
                    dose = (edep/eV)*1.60218e-19 / fMass;
                }
                //if (dose > 20*gray) {return false;} 

                fTotalDose += dose;
                fEnergyDepositPerEvent += edep ;
                fDosePerPulse += dose;

                if (fTotalDose >= fPrescribedDose/gray) {
                    aStep->GetTrack()->SetTrackStatus(fStopAndKill);
                    return false;
                }

                if (fDosePerPulse >= (fPrescribedDose/fNumberOfPulses)/gray) {
                    aStep->GetTrack()->SetTrackStatus(fStopAndKill);
                    //return false;
                }
                return true;
            }
        } 
    }

    else
        aStep->GetTrack()->SetTrackStatus(fStopAndKill);

    return false;
}


void TsScoreDNADamagePulsed::UserHookForPreTimeStepAction() {
    if (!G4EventManager::GetEventManager()->GetConstCurrentEvent()->IsAborted()) {

        FilterMoleculesInsideDNAVolumes();

        G4TrackManyList* trackList = G4ITTrackHolder::Instance()->GetMainList();
        G4ManyFastLists<G4Track>::iterator it_begin = trackList->begin();
        G4ManyFastLists<G4Track>::iterator it_end   = trackList->end();

        for(;it_begin!=it_end;++it_begin) {
            G4Molecule* Molecule = GetMolecule(*it_begin);
            const G4MoleculeDefinition* MolDef = Molecule->GetDefinition();
            if (MolDef == G4H2O::Definition()){return;}
        }

        SampleShiftTime();
        it_begin = trackList->begin();
        if (fShiftTime < fPulseTimeLimits[fCurrentPulse].first)
            fPulseTimeLimits[fCurrentPulse].first  = fShiftTime;

        if (fShiftTime > fPulseTimeLimits[fCurrentPulse].second)
            fPulseTimeLimits[fCurrentPulse].second = fShiftTime;

        for(;it_begin!=it_end;++it_begin) {
            TsIRTConfiguration::TsMolecule aMol = fIRT->ConstructMolecule(*it_begin, fShiftTime, it_begin->GetTrackID(), G4ThreeVector());
            aMol.parentID = fCurrentPulse;
            fMolecules[fCurrentPulse].push_back(aMol);
        }

        G4Scheduler::Instance()->Stop();
    }
}


void TsScoreDNADamagePulsed::UserHookForEndOfEvent() {
    if (fSimulationFinished) {return;}
    fTimeSampled = false;

    if ( fEnergyDepositPerEvent > 0 ) {
        fTotalEnergyDeposit += fEnergyDepositPerEvent;
        if (fPulseRecycle) {
            for (G4int i = 0; i < fNumberOfPulses; i++) {
                fVEnergyDepositPerEvent.push_back(std::make_pair(fShiftTime+(fPulsesTimeDelay*i),fEnergyDepositPerEvent));
            }
        }
        else
            fVEnergyDepositPerEvent.push_back(std::make_pair(fShiftTime,fEnergyDepositPerEvent));
        fNbOfScoredEvents++;
    }
    
    if ( fEnergyDepositPerEventEverywhere > 0 ) {
        fVEnergyDepositPerEventEverywhere.push_back(std::make_pair(fShiftTime,fEnergyDepositPerEventEverywhere));
        fNbOfScoredEventsEverywhere++;
    }

    fEnergyDepositPerEvent = 0.0;
    fEnergyDepositPerEventEverywhere = 0.0;

    if(fNumberOfPulses > 1 && fDosePerPulse >= (fPrescribedDose/fNumberOfPulses)/gray && fCurrentPulse < fNumberOfPulses && !fPulseRecycle) {
        fDosePerPulse = 0.0;
        fPulseTimeShift += fPulsesTimeDelay;
        fCurrentPulse++;
        fPulseTimeLimits[fCurrentPulse].first  =  1E150;
        fPulseTimeLimits[fCurrentPulse].second = 1E-150;
    }

    if(fTotalDose > fPrescribedDose/gray)
        RunAndSaveInfo();
    if (fPulseRecycle && fDosePerPulse >= (fPrescribedDose/fNumberOfPulses)/gray)
        RunAndSaveInfo();
}


void TsScoreDNADamagePulsed::AbsorbResultsFromWorkerScorer(TsVScorer* workerScorer) {
    TsVNtupleScorer::AbsorbResultsFromWorkerScorer(workerScorer);
    
    TsScoreDNADamagePulsed* workerGvalueScorer = dynamic_cast<TsScoreDNADamagePulsed*>(workerScorer);
    
    std::map<G4String, std::map<G4double,G4double>>::iterator wIter;
    std::map<G4String, std::map<G4double,G4double>>::iterator mIter;

    for(wIter = workerGvalueScorer->fGValuePerSpeciePerTime.begin(); wIter != workerGvalueScorer->fGValuePerSpeciePerTime.end(); wIter++) {
        mIter = fGValuePerSpeciePerTime.find(wIter->first);
        if (mIter == fGValuePerSpeciePerTime.end()) {
            fGValuePerSpeciePerTime.insert(std::pair<G4String,std::map<G4double,G4double>> (wIter->first,wIter->second));
        }
        else {
            std::map<G4double, G4double>::iterator witer;
            std::map<G4double, G4double>::iterator miter;
            for (witer = (wIter->second).begin(); witer != (wIter->second).end();witer++){
                miter = (mIter->second).find(witer->first);
                miter->second += witer->second;
            }
        }
    }

    for(wIter = workerGvalueScorer->fGValuePerSpeciePerTime2.begin(); wIter != workerGvalueScorer->fGValuePerSpeciePerTime2.end(); wIter++) {
        mIter = fGValuePerSpeciePerTime2.find(wIter->first);
        if (mIter == fGValuePerSpeciePerTime2.end()) {
            fGValuePerSpeciePerTime2.insert(std::pair<G4String,std::map<G4double,G4double>> (wIter->first,wIter->second));
        }
        else {
            std::map<G4double, G4double>::iterator witer;
            std::map<G4double, G4double>::iterator miter;
            for (witer = (wIter->second).begin(); witer != (wIter->second).end();witer++){
                miter = (mIter->second).find(witer->first);
                miter->second += witer->second;
            }
        }
    }

    for(wIter = workerGvalueScorer->fSpeciePerTime.begin(); wIter != workerGvalueScorer->fSpeciePerTime.end(); wIter++) {
        mIter = fSpeciePerTime.find(wIter->first);
        if (mIter == fSpeciePerTime.end()) {
            fSpeciePerTime.insert(std::pair<G4String,std::map<G4double,G4double>> (wIter->first,wIter->second));
        }
        else {
            std::map<G4double, G4double>::iterator witer;
            std::map<G4double, G4double>::iterator miter;
            for (witer = (wIter->second).begin(); witer != (wIter->second).end();witer++){
                miter = (mIter->second).find(witer->first);
                miter->second += witer->second;
            }
        }
    }

    for(wIter = workerGvalueScorer->fSpeciePerTime2.begin(); wIter != workerGvalueScorer->fSpeciePerTime2.end(); wIter++) {
        mIter = fSpeciePerTime2.find(wIter->first);
        if (mIter == fSpeciePerTime2.end()) {
            fSpeciePerTime2.insert(std::pair<G4String,std::map<G4double,G4double>> (wIter->first,wIter->second));
        }
        else {
            std::map<G4double, G4double>::iterator witer;
            std::map<G4double, G4double>::iterator miter;
            for (witer = (wIter->second).begin(); witer != (wIter->second).end();witer++){
                miter = (mIter->second).find(witer->first);
                miter->second += witer->second;
            }
        }
    }

    for (size_t i = 0; i < workerGvalueScorer->fAccumulatedEnergy.size();i++) {
        fAccumulatedEnergy[i]  += workerGvalueScorer->fAccumulatedEnergy[i];
        fAccumulatedEnergy2[i] += workerGvalueScorer->fAccumulatedEnergy2[i];
    }

    for (size_t i = 0; i < workerGvalueScorer->fPulseInformation.size(); i++) {
        fPulseInformation.push_back(workerGvalueScorer->fPulseInformation[i]);
    }

    for (auto& IndexTimeAndDelta:workerGvalueScorer->fDeltaGValues) {
        G4int   Index = IndexTimeAndDelta.first;

        for (auto& TimeAndDelta:IndexTimeAndDelta.second) {
            G4double Time = TimeAndDelta.first;
            G4double Delta = TimeAndDelta.second;
            fDeltaGValues[Index][Time]  += Delta;
        }
    }

    for (auto& IndexTimeAndDelta:workerGvalueScorer->fDeltaGValues2) {
        G4int   Index = IndexTimeAndDelta.first;

        for (auto& TimeAndDelta:IndexTimeAndDelta.second) {
            G4double Time = TimeAndDelta.first;
            G4double Delta = TimeAndDelta.second;
            fDeltaGValues2[Index][Time]  += Delta;
        }
    }

    for (auto& EventAndDirectSB : workerGvalueScorer->fDirectStrandBreaks) {
        G4int Event = EventAndDirectSB.first;
        std::vector<std::vector<G4int>> DirectSBs = EventAndDirectSB.second;

        for (size_t i = 0 ; i < DirectSBs.size(); i++)
            fDirectStrandBreaks[Event+fNbOfIRTRuns].push_back(DirectSBs[i]);
    }

    for (auto& EventAndIndirectSB : workerGvalueScorer->fIndirectStrandBreaks) {
        G4int Event = EventAndIndirectSB.first;
        std::vector<std::vector<G4int>> IndirectSBs = EventAndIndirectSB.second;

        for (size_t i = 0; i < IndirectSBs.size(); i++) {
            fIndirectStrandBreaks[Event+fNbOfIRTRuns].push_back(IndirectSBs[i]);
        }
    }

    fNbOfIRTRuns += workerGvalueScorer->fNbOfIRTRuns;

    workerGvalueScorer->fNbOfIRTRuns = 0;
    workerGvalueScorer->fGValuePerSpeciePerTime.clear();
    workerGvalueScorer->fGValuePerSpeciePerTime2.clear();
    workerGvalueScorer->fSpeciePerTime.clear();
    workerGvalueScorer->fSpeciePerTime2.clear();
    workerGvalueScorer->fPulseInformation.clear();
    workerGvalueScorer->fVEnergyDeposit.clear();
    workerGvalueScorer->fAccumulatedEnergy.clear();
    workerGvalueScorer->fAccumulatedEnergy2.clear();
    workerGvalueScorer->fDeltaGValues.clear();
    workerGvalueScorer->fDeltaGValues2.clear();
    workerGvalueScorer->fDirectStrandBreaks.clear();
    workerGvalueScorer->fIndirectStrandBreaks.clear();
}


void TsScoreDNADamagePulsed::Output() {
    if (fNbOfIRTRuns == 0) {
        return;
    }

    for (auto& NameAndTimeGValuePerTime: fGValuePerSpeciePerTime) {
        G4String Name = NameAndTimeGValuePerTime.first;

        G4int tBin = 0;
        //for (iter = (wIter->second).begin(); iter != (wIter->second).end(); iter++) {
        for (auto& TimeAndGValue:NameAndTimeGValuePerTime.second) {
            G4double Time = TimeAndGValue.first;
            G4double GVal = TimeAndGValue.second;
            G4double Mols = fSpeciePerTime[Name][Time];

            fGValue = GVal/fNbOfIRTRuns;
            fNumberOfMolecules = Mols/fNbOfIRTRuns;
            fEnergy  = fAccumulatedEnergy[tBin]/fNbOfIRTRuns;

            if (fNbOfIRTRuns > 1) {
                fGValueError = sqrt( (1.0/(fNbOfIRTRuns-1)) * ((fGValuePerSpeciePerTime2[Name][Time])/fNbOfIRTRuns - fGValue*fGValue));
                fNumberOfMoleculesError = sqrt( (1.0/(fNbOfIRTRuns-1)) * ((fSpeciePerTime2[Name][Time])/fNbOfIRTRuns - fNumberOfMolecules*fNumberOfMolecules));
                fEnergyError = sqrt( (1.0/(fNbOfIRTRuns-1)) * ((fAccumulatedEnergy2[tBin])/fNbOfIRTRuns - fEnergy*fEnergy));
            }
            else {
                fGValueError = 1;
                fNumberOfMoleculesError = 1;
                fEnergyError = 1;
            }
            fTime = Time;
            fMoleculeName = Name;
            fNtuple->Fill();
            tBin++;
        }
    }

    fNtuple->Write();

    if (fSavePulseInformation) {
        fTimeOutFile.open(fRootFileName, std::ios::binary | std::ios::app);
        for(size_t i = 0; i < fPulseInformation.size(); i++) {
            G4double saveTime = fPulseInformation[i].first/ps;
            G4double saveEdep = fPulseInformation[i].second/eV;
            fTimeOutFile.write(reinterpret_cast<char*>(&saveTime), sizeof saveTime);
            fTimeOutFile.write(reinterpret_cast<char*>(&saveEdep), sizeof saveEdep);
        }
        fTimeOutFile.close();
    }

    std::ofstream StrandBreakOut;
    StrandBreakOut.open(fOutputFile + "_StrandBreak.phsp");

    StrandBreakOut << "# TOPAS-nBio IRT: DNA Strand Break Map File"   << std::endl;
    StrandBreakOut << "# BreakID = 1 if Direct, 2 if Indirect SB"     << std::endl;
    StrandBreakOut << "# | EventID | VolumeID | BaseID | StrandID | MoleculeID | BreakID |" << std::endl;

    for (G4int i = 0; i < fNbOfIRTRuns; i++) {
        if (fDirectStrandBreaks.find(i) != fDirectStrandBreaks.end()) {
            for (size_t j = 0; j  < fDirectStrandBreaks[i].size(); j++) {
                StrandBreakOut << std::setw(12) << i
                               << std::setw(11) << fDirectStrandBreaks[i][j][0]
                               << std::setw(9)  << fDirectStrandBreaks[i][j][1]
                               << std::setw(11) << fDirectStrandBreaks[i][j][2]
                               << std::setw(12) << 0
                               << std::setw(10) << 1 << std::endl;          
            }
        }

        if (fIndirectStrandBreaks.find(i) != fIndirectStrandBreaks.end()) {
            for (size_t j = 0; j  < fIndirectStrandBreaks[i].size(); j++) {
                StrandBreakOut << std::setw(12) << i
                               << std::setw(11) << fIndirectStrandBreaks[i][j][0]
                               << std::setw(9)  << fIndirectStrandBreaks[i][j][1]
                               << std::setw(11) << fIndirectStrandBreaks[i][j][2]
                               << std::setw(12) << fIndirectStrandBreaks[i][j][3]
                               << std::setw(10) << 2 << std::endl;          
            }
        }
    }
    StrandBreakOut.close();
}


void TsScoreDNADamagePulsed::SampleShiftTime() {
    if ( fTimeDistributionType == 1 ) {
        while ( true ) {
            fShiftTime = G4RandGauss::shoot(fTimeMean+fPulseTimeShift, fTimeStdv);
            if ( fShiftTime > 1*ps ) {break;}
        }
    } else if ( fTimeDistributionType == 2 ) {
        while (true) {
            fShiftTime = G4RandFlat::shoot(fTimeMean - 0.5*fTimeFWHM + fPulseTimeShift, fTimeMean + 0.5*fTimeFWHM + fPulseTimeShift);
            if ( fShiftTime > 1*ps ) {break;}
        }
        fShiftTime += fPulseTimeShift;
    } else if ( fTimeDistributionType == 3 ) {
        G4double aRandom = G4UniformRand();
        G4int j = fNbOfTimes - 1;
        while ((fTimeTops[j] >= aRandom) && (j >= 0)) {
            j--;
        }
        fShiftTime = fTimeValues[j];
        fShiftTime += fPulseTimeShift;
    } else if ( fTimeDistributionType == 4 ) {
        while (true) {
            fShiftTime = G4RandExponential::shoot(fTimeMean + fPulseTimeShift);
            if ( fShiftTime >= 1.0*ps ) {break;}
        }
    } else {
        fShiftTime = 1.0 * ps;
    }
    
    if ( fLowLimitTo1ps && fTimeDistributionType != 3 && (fNbOfScoredEvents == 0 && fNbOfScoredEventsEverywhere == 0)) {
        fShiftTime = 1.0 * ps;
    }
}


void TsScoreDNADamagePulsed::RunAndSaveInfo() {
    G4cout << "Starting Chemistry after " << fNbOfScoredEvents << " primaries | Dose = " << fTotalDose << " Gy" << G4endl;
    RunChemistry();
    G4int tBin;
    for(size_t t = 0; t < fVEnergyDepositPerEvent.size();t++) {
        tBin = fUtils->FindBin(fVEnergyDepositPerEvent[t].first, fStepTimes);

        for (int i = tBin; i < (int)fStepTimes.size(); i++) {
            fVEnergyDeposit[i] += fVEnergyDepositPerEvent[t].second;
        }

        if (fSavePulseInformation) {
            fPulseInformation.push_back(fVEnergyDepositPerEvent[t]);
        }
    }

    for (size_t t = 0; t < fVEnergyDeposit.size();t++){
        fAccumulatedEnergy[t]  += fVEnergyDeposit[t]/MeV;
        fAccumulatedEnergy2[t] += (fVEnergyDeposit[t]/MeV)*(fVEnergyDeposit[t]/MeV);
    }

    for ( auto& nameTimeAndGvalue : fGValuesInRun ) {
        G4String name = nameTimeAndGvalue.first;
        for ( auto& timeAndGvalue : nameTimeAndGvalue.second ) {
            G4double time      = timeAndGvalue.first;
            G4double gvalue    = timeAndGvalue.second;
            G4double molecules = timeAndGvalue.second;
            tBin = fUtils->FindBin(time, fStepTimes);
            G4double Energy = fVEnergyDeposit[tBin];

            fSpeciePerTime[name][time]  += molecules;
            fSpeciePerTime2[name][time] += molecules*molecules;

            if (fVEnergyDeposit[tBin] > 0) {
                gvalue = 100*gvalue/(Energy/eV);
            }

            fGValuePerSpeciePerTime[name][time]  += gvalue;
            fGValuePerSpeciePerTime2[name][time] += gvalue*gvalue;
        }
    }

    for ( auto& indexTimeAndDeltaG : fDeltaGInRun ) {
        G4int index = indexTimeAndDeltaG.first;
        for ( auto& timeAndDeltaG : indexTimeAndDeltaG.second ) {
            G4double time      = timeAndDeltaG.first;
            G4double deltaG    = timeAndDeltaG.second;

            fDeltaGValues[index][time]  += deltaG;
            fDeltaGValues2[index][time] += deltaG*deltaG;
        }
    }
    
    fIRT->Clean();
    fMolecules.clear();
    fPulseTimeLimits.clear();
    fDNAHasBeenInserted = false;

    fTotalDose = 0;

    if ((fNumberOfRepetitions / fNumberOfThreads) <= fNbOfIRTRuns) {
        fSimulationFinished = true;
        G4cout << "Work At thread Finished" << G4endl;
    }

    fPulseTimeShift   = 0.0;
    fNbOfScoredEvents = 0;
    fNbOfScoredEventsEverywhere = 0;
    fVEnergyDepositPerEvent.clear();
    fVEnergyDepositPerEventEverywhere.clear();
    fCurrentPulse = 1;
    fDosePerPulse = 0.0;
    fTotalEnergyDeposit = 0;

    fPulseTimeLimits[fCurrentPulse].first  =  1E150;
    fPulseTimeLimits[fCurrentPulse].second = 1E-150;

    for (size_t enI = 0; enI < fVEnergyDeposit.size(); enI++) 
        fVEnergyDeposit[enI] = 0;    
}


void TsScoreDNADamagePulsed::RunChemistry() {
    fGValuesInRun.clear();
    fDeltaGInRun.clear();
    CheckForDirectDNABreaks();
    InsertDNAMolecules();
    G4double LastMinPulseTime = -1;
    G4double LastMaxPulseTime = -1;

    if (fNumberOfPulses == 1) {
        fTimeSimulationForPulse = fStepTimes[fStepTimes.size()-1];
    }

    std::map<G4int,G4String> MoleculeNames = fIRT->GetIRTConfiguration()->GetMoleculeNames();
    for (auto& PulseIndexAndTimeLimits:fPulseTimeLimits) {
        G4int cPulse = PulseIndexAndTimeLimits.first;
        G4double MinPulseTime = fPulseTimeLimits[cPulse].first;
        G4double MaxPulseTime = fPulseTimeLimits[cPulse].second;

        // Add the Molecules from the current Pulse to the IRT
        for (size_t MolIndex = 0; MolIndex < fMolecules[cPulse].size(); MolIndex++) {
            fIRT->AddMolecule(fMolecules[cPulse][MolIndex]);
        }
        G4cout << "---- IRT Starts for Pulse [" << cPulse << "] from " << MinPulseTime/us << " us to " << (MaxPulseTime+fTimeSimulationForPulse)/us << " us" << G4endl;

        if (LastMinPulseTime < 0)
            fIRT->runIRT(MinPulseTime, MaxPulseTime+(fTimeSimulationForPulse));
        else if (LastMinPulseTime < MinPulseTime) {
            fIRT->runIRT(MinPulseTime, MaxPulseTime+(fTimeSimulationForPulse));
        }
        else {
            fIRT->runIRT(LastMinPulseTime/(cPulse-1),(MaxPulseTime/(cPulse-1))+(fTimeSimulationForPulse));
        }

        G4cout << "---- IRT Ends for Pulse [" << cPulse << "]" << G4endl;
        LastMinPulseTime = MinPulseTime;
        LastMaxPulseTime = MaxPulseTime;
    }
    fGValuesInRun = fIRT->GetGValues();
    fDeltaGInRun  = fIRT->GetDeltaGValues();
 
    CheckForIndirectDNABreaks();
    fNbOfIRTRuns++;
}


void TsScoreDNADamagePulsed::GetDNAInformation() {
    G4String Component = fComponent;
    G4StrUtil::to_lower(Component);
    TsVGeometryComponent* DNAComponent     = fGm->GetComponent(Component);
    TsIRTPlasmidSupercoiled* PlasmidComponent1 = dynamic_cast<TsIRTPlasmidSupercoiled*>(DNAComponent);

    if (PlasmidComponent1 != 0) {
        fPHSPMolecules = PlasmidComponent1->GetDNANames();
        fPHSPPosition  = PlasmidComponent1->GetDNAPositions();
        fPHSPTime      = PlasmidComponent1->GetDNATimes();
        fDNADetails    = PlasmidComponent1->GetDNADetails();
    }

    else {
        G4cout << "No DNA Compatible scorer found!!!" << G4endl;
        exit(1);
    }

    for (size_t i = 0; i < fPHSPTime.size(); i++) {
        fPHSPIsAlive.push_back(true);
    }
}


void TsScoreDNADamagePulsed::InsertDNAMolecules() {
    if (fDNAHasBeenInserted) {return;}

    std::map<G4String,G4int> MolIDs = fIRT->GetIRTConfiguration()->GetMoleculeIDs();
    
    for (size_t i = 0; i < fPHSPMolecules.size(); i++) {
        if (fPHSPIsAlive[i]) {
            G4String MolName = SampleDNAMolecule(); 
            G4int    MolID   = MolIDs[MolName];
            TsIRTConfiguration::TsMolecule aMol;
            aMol.id       = MolID;
            aMol.position = fPHSPPosition[i];
            aMol.time     = 1 * ps;
            aMol.spin     = 0;
            aMol.trackID  = 100;
            aMol.parentID = 0;
            aMol.reacted  = false;
            aMol.isDNA    = false;
            aMol.isNew    = true;
            //aMol.exinfo   = fDNADetails[i];
            aMol.volumeID = fDNADetails[i][0];
            aMol.baseID   = fDNADetails[i][1];
            aMol.strandID = fDNADetails[i][2];
            fMolecules[1].push_back(aMol);
        }
        else
            fPHSPIsAlive[i] = true;
    }
    fDNAHasBeenInserted = true;
}


void TsScoreDNADamagePulsed::FilterMoleculesInsideDNAVolumes() {
    G4TrackManyList* trackList = G4ITTrackHolder::Instance()->GetMainList();
    G4ManyFastLists<G4Track>::iterator it_begin = trackList->begin();
    G4ManyFastLists<G4Track>::iterator it_end   = trackList->end();
    for(;it_begin!=it_end;++it_begin){
        if (it_begin->GetMaterial()->GetName() != "G4_WATER") {
            it_begin->SetTrackStatus(fStopAndKill);
        }
    }
}


void TsScoreDNADamagePulsed::CheckForDirectDNABreaks() {
    // Check for Direct Strand Breaks
    for ( auto& PlasmidAndElse : fDNAEnergyDeposit ) {
        for ( auto& BaseAndElse : PlasmidAndElse.second ) {
            for ( auto& StrandAndEnergy : BaseAndElse.second) {
                if ( StrandAndEnergy.second >= 17.5 * eV ) {
                    G4int PlasmidID = PlasmidAndElse.first;
                    G4int BaseID    = BaseAndElse.first;
                    G4int StrandID  = StrandAndEnergy.first;
                    std::vector<G4int> tempbreak;
                    tempbreak.push_back(PlasmidID);
                    tempbreak.push_back(BaseID);
                    tempbreak.push_back(StrandID);
                    fDirectStrandBreaksInEvent.push_back(tempbreak);
                }
            }
        }
    }

    if (fDNAEnergyDeposit.size() != 0)
        fDNAEnergyDeposit.clear();

    // Remove DNA Molecule of Strand Break from the Chemical Stage
    for (size_t i = 0; i < fDirectStrandBreaksInEvent.size(); i++) {
        for (size_t j = 0; j < fDNADetails.size(); j++) {
            if ((fDNADetails[j][0] == fDirectStrandBreaksInEvent[i][0]) &&
                (fDNADetails[j][1] == fDirectStrandBreaksInEvent[i][1]) &&
                (fDNADetails[j][2] == fDirectStrandBreaksInEvent[i][2])) {
                fPHSPIsAlive[j] = false;
                break;
            }
        }
    }
}


void TsScoreDNADamagePulsed::CheckForIndirectDNABreaks() {
    // Score the indirect strand breaks
    for (size_t i = 0; i < fStrandBreakMolecules.size(); i++) {
        G4int BreakMolID = (fIRT->GetIRTConfiguration()->GetMoleculeIDs())[fStrandBreakMolecules[i]];
        std::vector<TsIRTConfiguration::TsMolecule> SurvivingMolecules = fIRT->GetSurvivingMoleculesWithMolID(BreakMolID);
        for (size_t j = 0; j < SurvivingMolecules.size(); j++) {
            if (fIndirectStrandBreaks.count(fNbOfIRTRuns)){
                //std::vector<G4int> thisBreak = SurvivingMolecules[j].exinfo;
                std::vector<G4int> thisBreak;
                thisBreak.push_back(SurvivingMolecules[j].volumeID);
                thisBreak.push_back(SurvivingMolecules[j].baseID);
                thisBreak.push_back(SurvivingMolecules[j].strandID);
                thisBreak.push_back(int(i));
                fIndirectStrandBreaks[fNbOfIRTRuns].push_back(thisBreak);
            }
            else {
                //std::vector<G4int> thisBreak = SurvivingMolecules[j].exinfo;
                std::vector<G4int> thisBreak;
                thisBreak.push_back(SurvivingMolecules[j].volumeID);
                thisBreak.push_back(SurvivingMolecules[j].baseID);
                thisBreak.push_back(SurvivingMolecules[j].strandID);
                thisBreak.push_back(int(i));
                std::vector<std::vector<G4int>> tempBreak;
                tempBreak.push_back(thisBreak);
                fIndirectStrandBreaks[fNbOfIRTRuns] = tempBreak;
            }
        }
    }

    for (size_t i = 0; i < fDirectStrandBreaksInEvent.size(); i++)
        if (fDirectStrandBreaks.count(fNbOfIRTRuns))
            fDirectStrandBreaks[fNbOfIRTRuns].push_back(fDirectStrandBreaksInEvent[i]);
        else {
            std::vector<std::vector<G4int>> tempBreak;
            tempBreak.push_back(fDirectStrandBreaksInEvent[i]);
            fDirectStrandBreaks[fNbOfIRTRuns]= tempBreak;
        }

    if (fDirectStrandBreaksInEvent.size() != 0)
        fDirectStrandBreaksInEvent.clear();
}

G4String TsScoreDNADamagePulsed::SampleDNAMolecule() {
    if (fDNAMoleculesToInsert.size() == 1)
        return fDNAMoleculesToInsert[0];

    G4double Random = G4UniformRand();
    G4double Sum    = 0;

    for (size_t i = 0; i < fDNAMoleculesWeights.size(); i++) {
        Sum += fDNAMoleculesWeights[i];
        if (Random <= Sum) {
            return fDNAMoleculesToInsert[i];
        }
    }

    return fDNAMoleculesToInsert[fDNAMoleculesToInsert.size()-1];
}
