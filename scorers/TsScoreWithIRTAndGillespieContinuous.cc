// Scorer for IRTContinuous
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

#include "TsScoreWithIRTAndGillespieContinuous.hh"

#include "TsIRTManager.hh"
#include "TsVIRTProcedure.hh"
#include "TsIRTContinuous.hh"
#include "TsIRTConfiguration.hh"

#include "TsGillespie.hh"
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

TsScoreWithIRTAndGillespieContinuous::TsScoreWithIRTAndGillespieContinuous(TsParameterManager* pM, TsMaterialManager* mM, TsGeometryManager* gM, TsScoringManager* scM, TsExtensionManager* eM,
                                                     G4String scorerName, G4String quantity, G4String outFileName, G4bool isSubScorer)
: TsVNtupleScorer(pM, mM, gM, scM, eM, scorerName, quantity, outFileName, isSubScorer),
fPm(pM), fEnergyDepositPerEvent(0), fName(scorerName), fOldEvent(-1)
{
    SetUnit("");
    
    auto manager = new TsIRTManager(fPm, fName);
    fIRT       = (TsIRTContinuous*) manager->GetIRTProcedure();
    fUtils     = fIRT->GetUtils();
    fStepTimes = fIRT->GetStepTimes();
    
    fNtuple->RegisterColumnD(&fGValue,      "GValue", "");
    fNtuple->RegisterColumnD(&fGValueError, "GValue Error", "");
//    fNtuple->RegisterColumnD(&fEnergy,      "EnergyDeposit", "");
    //fNtuple->RegisterColumnD(&fEnergyError, "EnergyDeposit Error", "");
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
        /*
        fTimeTops.push_back(fTimeWeights[0]);
        for ( int i = 1; i < fNbOfTimes; i++ )
            fTimeTops.push_back(fTimeTops[i-1] + fTimeWeights[i]);*/
        
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
    fNbOfIRTRuns = 0;

    fTotalDose = 0.0;
    fDosePerPulse = 0.0;
    fPulseTimeShift = fPulsesTimeDelay;
    
    SampleShiftTime();
    fMinShiftTime = 1.0 * s;
    fCurrentPulse = 0;

    fSimulationFinished = false;

    fOutputFile = "TOPAS";
    if ( fPm->ParameterExists(GetFullParmName("OutputFile"))) {
        fOutputFile = fPm->GetStringParameter(GetFullParmName("OutputFile"));
    }

    fReportDelta = false;
    if ( fPm->ParameterExists(GetFullParmName("ReportDeltaGValues"))) {
        fReportDelta = fPm->GetBooleanParameter(GetFullParmName("ReportDeltaGValues"));
    }

    fPulseRecycle = false;
    if ( fPm->ParameterExists(GetFullParmName("RecyclePulses"))) {
        fPulseRecycle = fPm->GetBooleanParameter(GetFullParmName("RecyclePulses"));
    }

    if ( fPm->ParameterExists(GetFullParmName("IrradiationTime"))) {
    	fIrradiationTime = fPm->GetDoubleParameter(GetFullParmName("IrradiationTime"),"Time");
    }else{
    	fIrradiationTime = fStepTimes[fStepTimes.size()-1];
    }

    if ( fPm->ParameterExists("Ge/WorldLength")) {
    	fSideLength = fPm->GetDoubleParameter("Ge/WorldLength","Length");
    }

    fGillespieWorld = std::pow(fSideLength, 3);
}


TsScoreWithIRTAndGillespieContinuous::~TsScoreWithIRTAndGillespieContinuous()
{
    delete fIRT;
    delete fGillespie;
}


G4bool TsScoreWithIRTAndGillespieContinuous::ProcessHits(G4Step* aStep, G4TouchableHistory*)
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
                ResolveSolid(aStep);
                edep *= aStep->GetPreStepPoint()->GetWeight();
                
                G4double mass = 0.0;
                G4double density = aStep->GetPreStepPoint()->GetMaterial()->GetDensity();

                mass = (density * fSolid->GetCubicVolume())/kg;
                G4double dose = (edep/eV)*1.60218e-19 / mass;

                fEnergyDepositPerEvent += edep ;
                fDosePerPulse += dose;
                
                return true;
            }
        } 
    }

    else
        aStep->GetTrack()->SetTrackStatus(fStopAndKill);

    return false;
}


void TsScoreWithIRTAndGillespieContinuous::UserHookForPreTimeStepAction() {
    if (!G4EventManager::GetEventManager()->GetConstCurrentEvent()->IsAborted()) {

        G4TrackManyList* trackList = G4ITTrackHolder::Instance()->GetMainList();
        G4ManyFastLists<G4Track>::iterator it_begin = trackList->begin();
        G4ManyFastLists<G4Track>::iterator it_end   = trackList->end();

		for(;it_begin!=it_end;++it_begin) {
			G4Molecule* Molecule = GetMolecule(*it_begin);
			const G4MoleculeDefinition* MolDef = Molecule->GetDefinition();
	     	if (MolDef == G4H2O::Definition()) {
	     		return;
	     	}
		}

        it_begin = trackList->begin();

        for(;it_begin!=it_end;++it_begin) {
            TsIRTConfiguration::TsMolecule aMol = fIRT->ConstructMolecule(*it_begin, fShiftTime, fCurrentPulse, G4ThreeVector());
            fMolecules[fCurrentPulse].push_back(aMol);
        }
        G4Scheduler::Instance()->Stop();
    }
}


void TsScoreWithIRTAndGillespieContinuous::UserHookForEndOfEvent() {
    if (fSimulationFinished) {return;}

    if (fEnergyDepositPerEvent == 0 && fMolecules[fCurrentPulse].size() == 0){
        fDosePerPulse = 0.0;
        fEnergyDepositPerEvent = 0;
    	return;
    }

    fEnergyDep.push_back(std::make_pair(fCurrentPulse, fEnergyDepositPerEvent));
    fCurrentPulse++;
    fTotalDose += fDosePerPulse;

    fNbOfScoredEvents++;

    fEnergyDepositPerEvent = 0;
    fDosePerPulse = 0;

	if(fTotalDose > fPrescribedDose/gray){
		RunAndSaveInfo();
	}
}


void TsScoreWithIRTAndGillespieContinuous::AbsorbResultsFromWorkerScorer(TsVScorer* workerScorer) {
    TsVNtupleScorer::AbsorbResultsFromWorkerScorer(workerScorer);
    
    TsScoreWithIRTAndGillespieContinuous* workerGvalueScorer = dynamic_cast<TsScoreWithIRTAndGillespieContinuous*>(workerScorer);
    
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

    for (auto& IndexTimeAndDelta:workerGvalueScorer->fDeltaGValues_nh) {
        G4int   Index = IndexTimeAndDelta.first;

        for (auto& TimeAndDelta:IndexTimeAndDelta.second) {
            G4double Time = TimeAndDelta.first;
            G4double Delta = TimeAndDelta.second;
            fDeltaGValues_nh[Index][Time]  += Delta;
        }
    }

    for (auto& IndexTimeAndDelta:workerGvalueScorer->fDeltaGValues2_nh) {
        G4int   Index = IndexTimeAndDelta.first;

        for (auto& TimeAndDelta:IndexTimeAndDelta.second) {
            G4double Time = TimeAndDelta.first;
            G4double Delta = TimeAndDelta.second;
            fDeltaGValues2_nh[Index][Time]  += Delta;
        }
    }

    fNbOfIRTRuns += workerGvalueScorer->fNbOfIRTRuns;
    fTotalEnergyDeposit += workerGvalueScorer->fTotalEnergyDeposit;
    fNbOfScoredEvents += workerGvalueScorer->fNbOfScoredEvents;

    workerGvalueScorer->fNbOfIRTRuns = 0;
    workerGvalueScorer->fEnergyDep.clear();
    workerGvalueScorer->fGValuePerSpeciePerTime.clear();
    workerGvalueScorer->fGValuePerSpeciePerTime2.clear();
    workerGvalueScorer->fDeltaGValues.clear();
    workerGvalueScorer->fDeltaGValues2.clear();
    workerGvalueScorer->fDeltaGValues_nh.clear();
    workerGvalueScorer->fDeltaGValues2_nh.clear();
}

void TsScoreWithIRTAndGillespieContinuous::Output() {
    if (fNbOfScoredEvents == 0) {
        return;
    }

    G4double mass = fGillespieWorld * g/cm3;
    for (auto& NameAndTimeGValuePerTime: fGValuePerSpeciePerTime) {
        G4String Name = NameAndTimeGValuePerTime.first;

        for (auto& TimeAndGValue:NameAndTimeGValuePerTime.second) {
            G4double Time = TimeAndGValue.first;
            G4double GVal = TimeAndGValue.second;

            GVal = (G4double) GVal / (fPrescribedDose * mass) * (100 * eV);

            fGValue = GVal/fNumberOfThreads;
            G4double Gerror = fGValuePerSpeciePerTime2[Name][Time] / pow(fPrescribedDose * mass / (100 * eV), 2);
            fEnergy = 1;

            if (fNbOfScoredEvents > 1) {
                fGValueError = sqrt( (1.0/(fNumberOfThreads-1)) * (Gerror/fNumberOfThreads - fGValue*fGValue));
            }
            else {
                fGValueError = 1;
                fNumberOfMoleculesError = 1;
                fEnergyError = 1;
            }
            fTime = Time;
            fMoleculeName = Name;
            fNtuple->Fill();
        }
    }

    fNtuple->Write();

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

            DeltaReaction = deltaiter->second/fNumberOfThreads / (fPrescribedDose * mass) * (100 * eV);


            if ( fNbOfScoredEvents > 1 ) {
                DeltaError = sqrt( (1.0/(fNumberOfThreads-1)) * ( (deltaiter2->second)/fNumberOfThreads - DeltaReaction*DeltaReaction));
            } else {
                DeltaError = 1.0;
            }

            ReactionTime  = deltaiter->first;
            ReactionIndex = wDeltaIter->first;

            if(ReactionIndex >= 0){
                ReactA = (fIRT->GetReactants(ReactionIndex)).first;
                ReactB = (fIRT->GetReactants(ReactionIndex)).second;
                DeltaGFile << ReactionIndex+1 << "    "  << ReactA << "  "  << ReactB << "  ";
            }else{
            	ReactA = fIRT->GetIRTConfiguration()->GetMoleculeNames()[-ReactionIndex];
                DeltaGFile << ReactionIndex << "    "  << ReactA << "  "  << "IRT" << "  ";
            }

            std::vector<G4String> Products = fIRT->GetProducts(ReactionIndex);

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


    OutputFileName = fPm->GetStringParameter(GetFullParmName("OutputFile")) + "_DeltaG_NH.phsp";
    remove(OutputFileName.c_str());
    DeltaGFile.open(OutputFileName, std::ofstream::app);

    for ( wDeltaIter = fDeltaGValues_nh.begin(); wDeltaIter != fDeltaGValues_nh.end(); wDeltaIter++ ) {
        wDeltaIter2 = fDeltaGValues_nh.find(wDeltaIter->first );

        for ( deltaiter = (wDeltaIter->second).begin(); deltaiter != (wDeltaIter->second).end(); deltaiter++) {
            deltaiter2 = (wDeltaIter2->second).find( deltaiter->first );

            DeltaReaction = deltaiter->second/fNumberOfThreads / (fPrescribedDose * mass) * (100 * eV);


            if ( fNbOfScoredEvents > 1 ) {
                DeltaError = sqrt( (1.0/(fNumberOfThreads-1)) * ( (deltaiter2->second)/fNumberOfThreads - DeltaReaction*DeltaReaction));
            } else {
                DeltaError = 1.0;
            }

            ReactionTime  = deltaiter->first;
            ReactionIndex = wDeltaIter->first;

            if(ReactionIndex >= 0){
                ReactA = (fIRT->GetReactants(ReactionIndex)).first;
                ReactB = (fIRT->GetReactants(ReactionIndex)).second;
                DeltaGFile << ReactionIndex+1 << "    "  << ReactA << "  "  << ReactB << "  ";
            }else{
            	ReactA = fIRT->GetIRTConfiguration()->GetMoleculeNames()[-ReactionIndex];
                DeltaGFile << ReactionIndex << "    "  << ReactA << "  "  << "IRT" << "  ";
            }

            std::vector<G4String> Products = fIRT->GetProducts(ReactionIndex);

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


void TsScoreWithIRTAndGillespieContinuous::SampleShiftTime() {
    if ( fTimeDistributionType == 1 ) {
        while ( true ) {
            fShiftTime = G4RandGauss::shoot(fTimeMean, fTimeStdv);
            if ( fShiftTime > 1*ps ) {break;}
        }
    } else if ( fTimeDistributionType == 2 ) {
        while (true) {
            fShiftTime = G4RandFlat::shoot(fTimeMean- fTimeFWHM*0.5, fTimeMean + fTimeFWHM*0.5);
            if ( fShiftTime > 1*ps ) {break;}
        }
    }
    /*else if ( fTimeDistributionType == 3 ) {
        G4double aRandom = G4UniformRand();
        G4int j = fNbOfTimes - 1;
        while ((fTimeTops[j] >= aRandom) && (j >= 0)) {
            j--;
        }
        fShiftTime = fTimeValues[j];
    }*/
    else if ( fTimeDistributionType == 4 ) {
        while (true) {
            fShiftTime = G4RandExponential::shoot(fTimeMean);
            if ( fShiftTime >= 1.0*ps ) {break;}
        }
    } else {
        fShiftTime = 1.0 * ps;
    }

    if ( fLowLimitTo1ps && fTimeDistributionType != 3 && (fNbOfScoredEvents == 0 && fNbOfScoredEventsEverywhere == 0)) {
        fShiftTime = 1.0 * ps;
    }

}


void TsScoreWithIRTAndGillespieContinuous::RunAndSaveInfo() {
    G4cout << "Starting Chemistry after " << fNbOfScoredEvents << " primaries" << G4endl;
    RunChemistry();
    G4int tBin;

    for ( auto& nameTimeAndGvalue : fGValuesInRun ) {
        G4String name = nameTimeAndGvalue.first;
        G4int i = 0;
        for ( auto& timeAndGvalue : nameTimeAndGvalue.second ) {
            G4double time      = timeAndGvalue.first;
            G4double gvalue    = timeAndGvalue.second;
            G4double molecules = timeAndGvalue.second;
            tBin = fUtils->FindBin(time, fStepTimes);
            G4double Energy = fEnergyDep[time].second;

            if (Energy > 0) {
                gvalue = 100*gvalue/(Energy/eV);
            }

            fGValuePerSpeciePerTime[name][fStepTimes[0]]  += gvalue;
            fGValuePerSpeciePerTime2[name][fStepTimes[0]] += gvalue*gvalue;

            i++;
        }
    }

    fIRT->Clean();
    fMolecules.clear();

    fTotalDose = 0;

    if ((fNumberOfRepetitions / fNumberOfThreads) <= fNbOfIRTRuns) {
        fSimulationFinished = true;

        G4double DR = fPrescribedDose / fIrradiationTime;

        G4cout<<"Dose: "<<fPrescribedDose/gray<<"Gy Irradiation time: "<<fIrradiationTime/s<<"s Dose rate: "<<DR<<G4endl;
        G4cout<<"1 s: "<<s<<'\t'<<"1 Gy: "<<gray<<'\t'<<"1 kg: "<<kg<<G4endl;

        G4double TimeStart = fStepTimes[0];
        G4double TimeEnd   = fStepTimes[fStepTimes.size()-1];
        G4int Bins         = fStepTimes.size();

        fGillespie = new TsGillespie(fPm,fName,fIRT,TimeStart,TimeEnd,Bins,fGillespieWorld, DR, fIrradiationTime);

        G4cout<<"Gillespie start"<<G4endl;

        for ( auto& nameTimeAndGvalue : fGValuePerSpeciePerTime ) {
        	G4String name = nameTimeAndGvalue.first;
            G4int tBin = 0;
            for (auto& TimeAndGValue:nameTimeAndGvalue.second) {
                G4double Time = TimeAndGValue.first;
                G4double GVal = TimeAndGValue.second;

                G4double gvalue = GVal/fNbOfScoredEvents;
            	fGillespie->SetMolPerTime(fIRT->GetIRTConfiguration()->GetMoleculeIDs()[name], gvalue);
            }
        }

        fGillespie->RunAlt(TimeStart, TimeEnd);

        fGValuesInRun = fGillespie->GetGValues();
        fDeltaGInRun = fGillespie->GetDeltaGValues();


        for ( auto& indexTimeAndDeltaG : fDeltaGInRun ) {
            G4int index = indexTimeAndDeltaG.first;
            for ( auto& timeAndDeltaG : indexTimeAndDeltaG.second ) {
                G4double time      = timeAndDeltaG.first;
                G4double deltaG    = timeAndDeltaG.second;

                fDeltaGValues[index][time]  += deltaG;
                fDeltaGValues2[index][time] += deltaG*deltaG;
            }
        }


        for ( auto& indexTimeAndDeltaG : fDeltaGInRun_nh ) {
            G4int index = indexTimeAndDeltaG.first;
            for ( auto& timeAndDeltaG : indexTimeAndDeltaG.second ) {
                G4double time      = timeAndDeltaG.first;
                G4double deltaG    = timeAndDeltaG.second;

                fDeltaGValues_nh[index][time]  += deltaG;
                fDeltaGValues2_nh[index][time] += deltaG*deltaG;
            }
        }



        fGValuePerSpeciePerTime.clear();
        fGValuePerSpeciePerTime2.clear();

        for ( auto it = fGValuesInRun.begin(); it != fGValuesInRun.end(); it++ ) {
            G4String name = it->first;
            G4int i = 0;
            for ( auto& timeAndGvalue : it->second ) {
                G4double time      = timeAndGvalue.first;
                G4double gvalue    = timeAndGvalue.second;
                G4double molecules = timeAndGvalue.second;
                tBin = fUtils->FindBin(time, fStepTimes);

                fGValuePerSpeciePerTime[name][time]  += gvalue;
                fGValuePerSpeciePerTime2[name][time] += gvalue*gvalue;

                i++;
            }
        }


        G4cout << "Work At thread Finished" << G4endl;
    }
    fPulseTimeShift   = fPulsesTimeDelay;
    fCurrentPulse = 0;
    fDosePerPulse = 0.0;
    fTotalEnergyDeposit = 0.0;
    fEnergyDep.clear();
}


void TsScoreWithIRTAndGillespieContinuous::RunChemistry() {

    fGValuesInRun.clear();
    fDeltaGInRun.clear();
    fDeltaGInRun_nh.clear();
    G4double LastMinPulseTime = -1;
    G4double LastMaxPulseTime = -1;

	auto end = fEnergyDep.end();
	end--;
    std::vector<G4double> pulse_time;
	for(G4int i=0;i<=end->first;i++){

		SampleShiftTime();
		G4double time = floor(fIrradiationTime / fPulsesTimeDelay * G4UniformRand()) * fPulsesTimeDelay;

		pulse_time.push_back(time + fShiftTime);
	}
	std::sort(pulse_time.begin(), pulse_time.end());

	for (auto& PulseIndex:fEnergyDep) {
		G4int cPulse = PulseIndex.first;
		fEnergyDep[cPulse].first = pulse_time[cPulse];

		for (size_t MolIndex = 0; MolIndex < fMolecules[cPulse].size(); MolIndex++) {
			fMolecules[cPulse][MolIndex].time = pulse_time[cPulse];
			fIRT->AddMolecule(fMolecules[cPulse][MolIndex], cPulse);
		}
	}


	G4double TimeStart = fStepTimes[0];
	G4double TimeEnd   = fStepTimes[fStepTimes.size()-1];

	fIRT->runIRT(TimeStart, TimeEnd);

	fGValuesInRun = fIRT->GetGValues();
	fDeltaGInRun_nh = fIRT->GetDeltaGValues();

	fIRT->Clean();

	fNbOfIRTRuns++;
}
