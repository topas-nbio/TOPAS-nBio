// Scorer for IRTDoseInTOPAS2
//
// ********************************************************************
// *                                                                  *
// * This  code  implementation is the  intellectual property  of the *
// * TOPAS collaboration.                                             *
// * Use or redistribution of this code is not permitted without the  *
// * explicit approval of the TOPAS collaboration.                    *
// * Contact: Joseph Perl, perl@slac.stanford.edu                     *
// *                                                                  *
// ********************************************************************
//

#include "TsScoreWithIRTCummulativeInTOPAS2.hh"
#include "TsHybridIRT.hh"
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
//#include <algorithm>
//#include <stdlib.h>
// #include <algorithm>

TsScoreWithIRTCummulativeInTOPAS2::TsScoreWithIRTCummulativeInTOPAS2(TsParameterManager* pM, TsMaterialManager* mM, TsGeometryManager* gM, TsScoringManager* scM, TsExtensionManager* eM,
                                                     G4String scorerName, G4String quantity, G4String outFileName, G4bool isSubScorer)
: TsVNtupleScorer(pM, mM, gM, scM, eM, scorerName, quantity, outFileName, isSubScorer),
fPm(pM), fEnergyDepositPerEvent(0), fEnergyDepositPerEventEverywhere(0), fName(scorerName), fOldEvent(-1)
{
    SetUnit("");
    
    fIRT       = new TsHybridIRT(fPm, fName);
    fUtils     = fIRT->GetUtils();
    fStepTimes = fIRT->GetStepTimes();
    
    fNtuple->RegisterColumnD(&fGValue,      "GValue", "");
    fNtuple->RegisterColumnD(&fGValueError, "GValue Error", "");
    fNtuple->RegisterColumnD(&fNumberOfMolecules, "Total Number of molecules","");
    fNtuple->RegisterColumnD(&fNumberOfMoleculesError, "Error of Number of molecules","");
    fNtuple->RegisterColumnD(&fEnergy,      "EnergyDeposit", "");
    //fNtuple->RegisterColumnD(&fEnergyError, "EnergyDeposit Error", "");
    fNtuple->RegisterColumnD(&fTime,        "Time", "picosecond");
    fNtuple->RegisterColumnS(&fMoleculeName, "MoleculeName");
    
    fPrescribedDose = fPm->GetDoubleParameter(GetFullParmName("PrescribedDose"),"Dose");

    fTimeDistribution = fPm->GetStringParameter(GetFullParmName("PulseDistribution"));
    fTimeDistribution.toLower();

    fNumberOfPulses  = 1;
    fPulsesTimeDelay = 0;
    if ( fPm->ParameterExists(GetFullParmName("NumberOfPulses")) ) {
        fNumberOfPulses = fPm->GetIntegerParameter(GetFullParmName("NumberOfPulses"));
        //if (fNumberOfPulses > 1) {
            G4double pulsesFrequency = fPm->GetDoubleParameter(GetFullParmName("PulsesFrequency"),"perTime");
            if ( fNumberOfPulses > 1 && pulsesFrequency == 0 ) {
                G4cerr << "TOPAS is exiting due to an error in parameter" << G4endl;
                G4cerr << GetFullParmName("PulsesFrequency") << G4endl;
                G4cerr << "The number of pulses is bigger than 1, therefore, the frequency must be larger than zero" << G4endl;
            }
            fPulsesTimeDelay = 1./pulsesFrequency;
        //}
        //else 
        //fPulsesTimeDelay = 1 * ps;
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
    fPulseTimeShift = fPulsesTimeDelay;
    
    SampleShiftTime();
    fMinShiftTime = 1.0 * s;
    fCurrentPulse = 1;

    fPulseTimeLimits[fCurrentPulse].first  =  1E150;
    fPulseTimeLimits[fCurrentPulse].second = 1E-150;

    fSimulationFinished = false;
    fTimeSampled        = false;

    fVolume = 0;

    G4double TimeStart = fStepTimes[0];
    G4double TimeEnd   = fStepTimes[fStepTimes.size()-1];
    G4int Bins         = fStepTimes.size();

    for ( size_t u = 0; u < fStepTimes.size(); u++ ) {
        fVEnergyDeposit.push_back(0.0);
        fAccumulatedEnergy.push_back(0.0);
        fAccumulatedEnergy2.push_back(0.0);
    }

    fOutputFile = "TOPAS";
    if ( fPm->ParameterExists(GetFullParmName("OutputFile"))) {
        fOutputFile = fPm->GetStringParameter(GetFullParmName("OutputFile"));
    }

    fReportEachHistory = false;
    if ( fPm->ParameterExists(GetFullParmName("ReportAtEachHistory"))) {
        fReportEachHistory = fPm->GetBooleanParameter(GetFullParmName("ReportAtEachHistory"));
    }

    fReportDelta = false;
    if ( fPm->ParameterExists(GetFullParmName("ReportDeltaGValues"))) {
        fReportDelta = fPm->GetBooleanParameter(GetFullParmName("ReportDeltaGValues"));
    }
}


TsScoreWithIRTCummulativeInTOPAS2::~TsScoreWithIRTCummulativeInTOPAS2()
{
    delete fIRT;
}


G4bool TsScoreWithIRTCummulativeInTOPAS2::ProcessHits(G4Step* aStep, G4TouchableHistory*)
{
    if (!fIsActive) {
        fSkippedWhileInactive++;
        return false;
    }

    if (!fSimulationFinished) {
        G4int TrckID  = aStep->GetTrack()->GetTrackID();

        if (!fTimeSampled) {
            fTimeSampled = true;
            SampleShiftTime();
        }

        if ( -1 < TrckID ) {
            ResolveSolid(aStep);
            fVolume = fSolid->GetCubicVolume();
            G4double edep = aStep->GetTotalEnergyDeposit() ;
            if ( edep > 0 ) {
                ResolveSolid(aStep);
                edep *= aStep->GetPreStepPoint()->GetWeight();
                
                G4double mass = 0.0;
                G4double density = aStep->GetPreStepPoint()->GetMaterial()->GetDensity();
                fEnergyDepositPerEventEverywhere += edep ;

                if ( fTestIsInside ) {
                    G4TouchableHistory* touchable = (G4TouchableHistory*)(aStep->GetPreStepPoint()->GetTouchable());
                    const G4String& volumeName = touchable->GetVolume()->GetName();
                    if( !volumeName.contains(fSensitiveVolume)) {
                        return false;
                    }
                }
                
                mass = (density * fSolid->GetCubicVolume())/kg;
                G4double dose = (edep/eV)*1.60218e-19 / mass;
                fTotalDose += dose;
                fEnergyDepositPerEvent += edep ;
                fDosePerPulse += dose;

                if (fTotalDose >= fPrescribedDose/gray) {
                    aStep->GetTrack()->SetTrackStatus(fStopAndKill);
                }

                if (fDosePerPulse >= (fPrescribedDose/fNumberOfPulses)/gray) {
                    aStep->GetTrack()->SetTrackStatus(fStopAndKill);
                }
                
                return true;
            }
        } 
    }

    else
        aStep->GetTrack()->SetTrackStatus(fStopAndKill);

    return false;
}


void TsScoreWithIRTCummulativeInTOPAS2::UserHookForPreTimeStepAction() {
    if (!G4EventManager::GetEventManager()->GetConstCurrentEvent()->IsAborted()) {

        if (G4Scheduler::Instance()->GetStatus() == 0) {
            return;
        }

        G4TrackManyList* trackList = G4ITTrackHolder::Instance()->GetMainList();
        G4ManyFastLists<G4Track>::iterator it_begin = trackList->begin();
        G4ManyFastLists<G4Track>::iterator it_end   = trackList->end();

        if (fShiftTime < fPulseTimeLimits[fCurrentPulse].first)
            fPulseTimeLimits[fCurrentPulse].first  = fShiftTime;

        if (fShiftTime > fPulseTimeLimits[fCurrentPulse].second)
            fPulseTimeLimits[fCurrentPulse].second = fShiftTime;

        for(;it_begin!=it_end;++it_begin) {
            TsIRTConfiguration::TsMolecule aMol = fIRT->ConstructMolecule(*it_begin, fShiftTime, fCurrentPulse, G4ThreeVector());
            fMolecules[fCurrentPulse].push_back(aMol);
        }
        G4Scheduler::Instance()->Stop();
    }
}


void TsScoreWithIRTCummulativeInTOPAS2::UserHookForEndOfEvent() {
    if (fSimulationFinished) {return;}
    fTimeSampled = false;
    //G4cout << "End Of Event: " << fNbOfScoredEvents << "  " << fTotalDose << "  " << fDosePerPulse << "  " << fCurrentPulse << "  " << (fPrescribedDose/fNumberOfPulses)/gray << "  " << fShiftTime << G4endl;

    if ( fEnergyDepositPerEvent > 0 ) {
        fTotalEnergyDeposit += fEnergyDepositPerEvent;
        fVEnergyDepositPerEvent.push_back(std::make_pair(fShiftTime,fEnergyDepositPerEvent));
        fNbOfScoredEvents++;
    }
    
    if ( fEnergyDepositPerEventEverywhere > 0 ) {
        fVEnergyDepositPerEventEverywhere.push_back(std::make_pair(fShiftTime,fEnergyDepositPerEventEverywhere));
        fNbOfScoredEventsEverywhere++;
    }

    fEnergyDepositPerEvent = 0.0;
    fEnergyDepositPerEventEverywhere = 0.0;
    G4double prescribedDosePerPulse = (fPrescribedDose/fNumberOfPulses)/gray;
    G4double doseDifference = fabs(fDosePerPulse - prescribedDosePerPulse);

    if(fNumberOfPulses > 1 && fDosePerPulse >= prescribedDosePerPulse && fCurrentPulse < fNumberOfPulses) {
        fDosePerPulse = 0.0;
        fPulseTimeShift += fPulsesTimeDelay;
        fCurrentPulse++;
        fPulseTimeLimits[fCurrentPulse].first  =  DBL_MAX;
        fPulseTimeLimits[fCurrentPulse].second = -DBL_MAX;
    }
    else if (fNumberOfPulses > 1 && fDosePerPulse >= prescribedDosePerPulse && fCurrentPulse == fNumberOfPulses) {
        RunAndSaveInfo();
    }

    else if(fNumberOfPulses == 1 && fDosePerPulse >= prescribedDosePerPulse) {
        RunAndSaveInfo();
    }
}


void TsScoreWithIRTCummulativeInTOPAS2::AbsorbResultsFromWorkerScorer(TsVScorer* workerScorer) {
    TsVNtupleScorer::AbsorbResultsFromWorkerScorer(workerScorer);
    
    TsScoreWithIRTCummulativeInTOPAS2* workerGvalueScorer = dynamic_cast<TsScoreWithIRTCummulativeInTOPAS2*>(workerScorer);
    
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
}


void TsScoreWithIRTCummulativeInTOPAS2::OutputForThread(G4int Thread) {
    std::ofstream PartialFile;
    G4String str = std::to_string(Thread);
    G4String FileName = fOutputFile + "_Partial_ThreadID" + str + ".phsp";
    remove(FileName.c_str());
    PartialFile.open(FileName, std::ofstream::app);

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
            PartialFile << "  " << fGValue << "     " << fGValueError << "  " << fNumberOfMolecules << "  " << fNumberOfMoleculesError << "  " << fEnergy << "  " << fTime / ps << "  " << fMoleculeName << G4endl;
            tBin++;
        }
    }
    PartialFile.close();
}


void TsScoreWithIRTCummulativeInTOPAS2::Output() {
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

    if (!fReportDelta) {return;}

    std::ofstream DeltaGFile;

    std::map<G4int, std::map<G4double, G4double> >::iterator wDeltaIter;
    std::map<G4int, std::map<G4double, G4double> >::iterator wDeltaIter2;
    std::map<G4double, G4double>::iterator deltaiter;
    std::map<G4double, G4double>::iterator deltaiter2;

    G4int ReactionIndex;
    G4double ReactionTime;
    G4double DeltaReaction;
    G4double DeltaError;
    G4String ReactA;
    G4String ReactB;

    G4String OutputFileName = fPm->GetStringParameter(GetFullParmName("OutputFile")) + "_DeltaG.phsp";
    remove(OutputFileName.c_str());
    DeltaGFile.open(OutputFileName, std::ofstream::app);

    for ( wDeltaIter = fDeltaGValues.begin(); wDeltaIter != fDeltaGValues.end(); wDeltaIter++ ) {
        wDeltaIter2 = fDeltaGValues.find(wDeltaIter->first );

        for ( deltaiter = (wDeltaIter->second).begin(); deltaiter != (wDeltaIter->second).end(); deltaiter++) {
            deltaiter2 = (wDeltaIter2->second).find( deltaiter->first );

            DeltaReaction = deltaiter->second/fNbOfIRTRuns;
            if ( fNbOfScoredEvents > 1 ) {
                DeltaError = sqrt( (1.0/(fNbOfIRTRuns-1)) * ( (deltaiter2->second)/fNbOfIRTRuns - DeltaReaction*DeltaReaction));
            } else {
                DeltaError = 1.0;
            }

            ReactionTime  = deltaiter->first;
            ReactionIndex = wDeltaIter->first;
            ReactA = (fIRT->GetReactants(ReactionIndex)).first;
            ReactB = (fIRT->GetReactants(ReactionIndex)).second;
            std::vector<G4String> Products = fIRT->GetProducts(ReactionIndex);
            DeltaGFile << ReactionIndex+1 << "    "  << ReactA << "  "  << ReactB << "  ";

                while (Products.size() < 3) {
                    Products.push_back("None");
                }

                for (size_t prod = 0; prod < Products.size(); prod++) {
                    DeltaGFile << Products[prod] << "  ";
                }

            DeltaGFile << "  " << ReactionTime * 1000 << "     " << DeltaReaction << "    " << DeltaError << std::endl;
        }
    }

    DeltaGFile.close();
}


void TsScoreWithIRTCummulativeInTOPAS2::SampleShiftTime() {
    if ( fTimeDistributionType == 1 ) {
        while ( true ) {
            fShiftTime = G4RandGauss::shoot(fTimeMean+fPulseTimeShift, fTimeStdv);
            if ( fShiftTime > 1*ps ) {break;}
        }
    } else if ( fTimeDistributionType == 2 ) {
        while (true) {
            //fShiftTime = G4RandFlat::shoot(fTimeMean- fTimeFWHM*0.5 + fPulseTimeShift, fTimeMean + fTimeFWHM*0.5 + fPulseTimeShift);
            fShiftTime = G4RandFlat::shoot(fTimeMean + fPulseTimeShift, fTimeMean + fTimeFWHM + fPulseTimeShift);
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
    
    //if ( fLowLimitTo1ps && fTimeDistributionType != 3 && (fNbOfScoredEvents == 0 && fNbOfScoredEventsEverywhere == 0)) {
    //    fShiftTime = 1.0 * ps;
    //}
}


void TsScoreWithIRTCummulativeInTOPAS2::RunAndSaveInfo() {
    G4cout << "Starting Chemistry after " << fNbOfScoredEvents << " primaries" << G4endl;
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

    if (fReportEachHistory) {
        OutputForThread(G4Threading::G4GetThreadId());
    }

    fTotalDose = 0;

    if ((fNumberOfRepetitions / fNumberOfThreads) <= fNbOfIRTRuns) {
        fSimulationFinished = true;
        G4cout << "Work At thread Finished" << G4endl;
    }
    fPulseTimeShift   = fPulsesTimeDelay;
    fNbOfScoredEvents = 0;
    fNbOfScoredEventsEverywhere = 0;
    fVEnergyDepositPerEvent.clear();
    fVEnergyDepositPerEventEverywhere.clear();
    fCurrentPulse = 1;
    fDosePerPulse = 0.0;
    fTotalEnergyDeposit = 0;

    fPulseTimeLimits[fCurrentPulse].first  =  DBL_MAX;
    fPulseTimeLimits[fCurrentPulse].second = -DBL_MAX;

    for (size_t enI = 0; enI < fVEnergyDeposit.size(); enI++) 
        fVEnergyDeposit[enI] = 0;    
}


void TsScoreWithIRTCummulativeInTOPAS2::RunChemistry() {
    fGValuesInRun.clear();
    fDeltaGInRun.clear();
    std::map<G4int,G4String> MoleculeNames = fIRT->GetIRTConfiguration()->GetMoleculeNames();

    for (auto& PulseIndexAndTimeLimits:fPulseTimeLimits) {
        G4int cPulse = PulseIndexAndTimeLimits.first;
        G4double MinPulseTime = fPulseTimeLimits[cPulse].first;
        G4double MaxPulseTime = fPulseTimeLimits[cPulse].second;

        if (fNumberOfPulses == cPulse) {
            MaxPulseTime = fStepTimes[fStepTimes.size()-1];
        }
        else {
            MaxPulseTime = fPulseTimeLimits[cPulse+1].first;
        }

        // Add the Molecules from the current Pulse to the IRT
        for (size_t MolIndex = 0; MolIndex < fMolecules[cPulse].size(); MolIndex++) {
            fIRT->AddMolecule(fMolecules[cPulse][MolIndex]);
        }
        G4cout << "---- IRT Starts for Pulse [" << cPulse << "] from " << MinPulseTime/us << " us to " << MaxPulseTime/us << " us | Pulse contributed with " << fMolecules[cPulse].size() << " Species" << G4endl;

        fIRT->runIRT(MinPulseTime, MaxPulseTime, fPulseTimeLimits[cPulse].second + (100*us));
        fIRT->SetContainersForNextPulse();
        G4cout << "---- IRT Ends for Pulse [" << cPulse << "]" << G4endl;
    }

    fGValuesInRun = fIRT->GetGValues();
    fDeltaGInRun =  fIRT->GetDeltaGValues();
    fIRT->Clean();
    fNbOfIRTRuns++;
}