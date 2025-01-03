// Scorer for TsIRTPulses
//
// ********************************************************************
// *                                                                  *
// * This file is part of the TOPAS-nBio extensions to the            *
// *   TOPAS Simulation Toolkit.                                      *
// * The TOPAS-nBio extensions are freely available under the license *
// *   agreement set forth at: https://topas-nbio.readthedocs.io/     *
// *                                                                  *
// ********************************************************************

#include "TsScoreWithIRTPulses.hh"
#include "TsIRTManager.hh"
#include "TsVIRTProcedure.hh"
#include "TsIRTUtils.hh"

#include "G4ITTrackHolder.hh"
#include "G4EventManager.hh"
#include "G4RunManager.hh"
#include "G4SystemOfUnits.hh"
#include "G4Scheduler.hh"
#include "G4Molecule.hh"
#include "G4H2O.hh"
#include "G4UnitsTable.hh"
#include "Randomize.hh"

TsScoreWithIRTPulses::TsScoreWithIRTPulses(TsParameterManager* pM, TsMaterialManager* mM, TsGeometryManager* gM, TsScoringManager* scM, TsExtensionManager* eM,
                                           G4String scorerName, G4String quantity, G4String outFileName, G4bool isSubScorer)
: TsVNtupleScorer(pM, mM, gM, scM, eM, scorerName, quantity, outFileName, isSubScorer),
fPm(pM), fEnergyDepositPerEvent(0)
{
    SetUnit("");
    
    if (fPm->GetIntegerParameter("Ts/NumberOfThreads") != 1) {
        G4cerr << "TOPAS is exiting due to an error in parameter." << G4endl;
        G4cerr << "Ts/NumberOfThreads has a value different than 1" << G4endl;
        G4cerr << "This example only can be run in single thread mode" << G4endl;
        G4cerr << "Set: i:Ts/NumberOfThreads = 1" << G4endl;
        fPm->AbortSession(1);
    }
    
    fNtuple->RegisterColumnD(&fGValue, "GValue: number of molecules per 100 eV of energy deposit", "");
    fNtuple->RegisterColumnD(&fTime,    "Time", "picosecond");
    fNtuple->RegisterColumnS(&fMoleculeName, "MoleculeName");
    fNtuple->RegisterColumnD(&fMolecules,"Number of molecules","");
    fNtuple->RegisterColumnD(&fEnergyDepositGvalue,"Energy deposited","eV");
    
    fIRT       = new TsIRTManager(fPm, scorerName);
    fUtils       = fIRT->GetUtils();
    fVStepTimes = fIRT->GetStepTimes();
    
    fEnergyDepositPerEvent = 0.0;
    fTotalEnergyDeposit = 0.0;
    
    fPrescribedDose = fPm->GetDoubleParameter(GetFullParmName("PrescribedDose"),"Dose");
    fDosePerPulse = 0*gray;
    fMass = 0*g;
    fDensity = 1 * g/(cm*cm*cm);
    fShiftTime = 1 * ps;
    fEventID = 0;
    fOldEvent = -1;
    
    fTimeMean = fPm->GetDoubleParameter(GetFullParmName("PulseTimeMean"),"Time");
    fTimeFWHM = fPm->GetDoubleParameter(GetFullParmName("PulseTimeFWHM"),"Time");
    fTimeStdv = fTimeFWHM/2.354820045;
    
    fNumberOfPulses = 1;
    fPulseCount = 0;
    fPulseTimeShift = 0*ps;
    fPulsesTimePeriod = 0*ps;
    if ( fPm->ParameterExists(GetFullParmName("NumberOfPulses")) ) {
        fNumberOfPulses = fPm->GetIntegerParameter(GetFullParmName("NumberOfPulses"));
        if ( fNumberOfPulses > 1  ) {
            G4double pulsesFrequency = fPm->GetDoubleParameter(GetFullParmName("PulsesFrequency"),"perTime");
            if ( pulsesFrequency <= 0 ) {
                G4cerr << "TOPAS is exiting due to an error in parameter" << G4endl;
                G4cerr << GetFullParmName("PulsesFrequency") << G4endl;
                G4cerr << "The number of pulses is greater than 1, therefore, the frequency must be larger than zero" << G4endl;
                fPm->AbortSession(1);
            }
            fPulsesTimePeriod = 1./pulsesFrequency;
        }
    }
    
    G4String fileName = fPm->GetStringParameter("Sc/RootFileName") + ".bin";
    remove(fileName);
    fTimeOutFile.open(fileName, std::ios::binary);
    
    for ( size_t u = 0; u < fVStepTimes.size(); u++ )
        fVEnergyDeposits.push_back(0.0);
}


TsScoreWithIRTPulses::~TsScoreWithIRTPulses()
{
    delete fIRT;
}


G4bool TsScoreWithIRTPulses::ProcessHits(G4Step* aStep, G4TouchableHistory*)
{
    if (!fIsActive) {
        fSkippedWhileInactive++;
        return false;
    }
    
    if ( -1 < aStep->GetTrack()->GetTrackID() ) {
        G4double edep = aStep->GetTotalEnergyDeposit();
        if ( edep > 0 ) {
            ResolveSolid(aStep);
            
            if (fIRT->GetIRTProcedure()->GetTestIsInside()) {
                if ( 0 == fMass )
                    fMass = fDensity * fIRT->GetIRTProcedure()->GetVirtualRegionCubicVolume();
                
                if (fIRT->GetIRTProcedure()->Inside(aStep->GetPreStepPoint()->GetPosition())) {
                    fEnergyDepositPerEvent += edep * aStep->GetPreStepPoint()->GetWeight();
                    return true;
                } else {
                    return false;
                }
                
            } else {
                if ( 0 == fMass )
                    fMass = fDensity * fSolid->GetCubicVolume();
                
                fEnergyDepositPerEvent += edep * aStep->GetPreStepPoint()->GetWeight();
                return true;
            }
        }
    }
    
    return false;
}


void TsScoreWithIRTPulses::UserHookForEndOfEvent() {
    if (fEnergyDepositPerEvent == 0 )
        return;
    
    if(fEnergyDepositPerEvent > 0) {
        SampleTimeShift();
        fVEnergyDepositPerSampledTime.push_back(std::make_pair(fShiftTime,fEnergyDepositPerEvent));
        fTotalEnergyDeposit += fEnergyDepositPerEvent;
        fDosePerPulse += fEnergyDepositPerEvent/fMass;
        fEnergyDepositPerEvent = 0;
        fEventID++; // only events that contributed to the dose are considered
    }
    
    // If several pulses, then dose is split in dose/numberOfPulses
    if(fNumberOfPulses > 1 && fDosePerPulse >= fPrescribedDose/fNumberOfPulses) {
        fPulseCount++;
        G4cout << "-- The pulse (" << fPulseCount << ") started at " << fPulseTimeShift/ps
        << " ps achieved dose per pulse " << fDosePerPulse/gray << " Gy. " << G4endl;
        fPulseTimeShift += fPulsesTimePeriod;
        fDosePerPulse = 0.0;
    }
    
    if ( fTotalEnergyDeposit/fMass >= fPrescribedDose ) {
        fPulseCount++;
        G4cout << "-- With pulse (" << fPulseCount << ") started at " << fPulseTimeShift/ps
        << " ps achieved a Total Dose " << (fTotalEnergyDeposit/fMass)/gray << " Gy after " << fEventID << " histories." << G4endl;
        
        for ( size_t t = 0; t < fVEnergyDepositPerSampledTime.size(); t++ ) {
            G4int tBin = fUtils->FindBin(fVEnergyDepositPerSampledTime[t].first, fVStepTimes);
            
            for ( int i = tBin; i < (int)fVStepTimes.size(); i++)
                fVEnergyDeposits[i] += fVEnergyDepositPerSampledTime[t].second;
        }
        
        fIRT->runIRT();
        
        std::map<G4String, std::map<G4double, G4int>> irt = fIRT->GetGValues();
        for ( auto& nameTimeAndGvalue : irt ) {
            G4String name = nameTimeAndGvalue.first;
            for ( auto& timeAndGvalue : (nameTimeAndGvalue.second) ) {
                G4double time = timeAndGvalue.first;
                G4double gvalue = timeAndGvalue.second;
                G4int tBin = fUtils->FindBin(time, fVStepTimes);
                
                fMoleculesPerSpeciePerTime[name][time] += gvalue;
                
                gvalue *= 100/(fVEnergyDeposits[tBin]/eV);
                fGValuePerSpeciePerTime[name][time] += gvalue;
            }
        }
        
        std::map<G4int, std::map<G4double, G4int>> DeltaG = fIRT->GetDeltaGValues();
        for ( auto& indexTimeAndDeltaG : DeltaG) {
            G4int reactionIndex = indexTimeAndDeltaG.first;
            for ( auto& timeAndDeltaG : (indexTimeAndDeltaG.second) ) {
                G4double time   = timeAndDeltaG.first;
                G4double deltaG = timeAndDeltaG.second;
                G4int tBin = fUtils->FindBin(time, fVStepTimes);
                
                fDeltaGPerReactionPerTime[reactionIndex][time] += deltaG;
            }
        }
        
        DeltaG.clear();
        irt.clear();
        fIRT->Clean();
        Output();
        
        G4RunManager::GetRunManager()->AbortRun(true);
    }
    
    return;
}


void TsScoreWithIRTPulses::UserHookForPreTimeStepAction() {
    if (0 == fEnergyDepositPerEvent) return;
    
    G4TrackManyList* trackList = G4ITTrackHolder::Instance()->GetMainList();
    G4ManyFastLists<G4Track>::iterator it_begin = trackList->begin();
    G4ManyFastLists<G4Track>::iterator it_end   = trackList->end();
    for(;it_begin!=it_end;++it_begin) {
        G4Molecule* Molecule = GetMolecule(*it_begin);
        const G4MoleculeDefinition* MolDef = Molecule->GetDefinition();
        if (MolDef == G4H2O::Definition())
            return;
    }
    
    it_begin = trackList->begin();
    for(;it_begin!=it_end;++it_begin){
        fIRT->AddMolecule(*it_begin, fShiftTime, 0, G4ThreeVector());
    }
    
    G4Scheduler::Instance()->Stop();
}


void TsScoreWithIRTPulses::SampleTimeShift() {
    //fEventID = GetEventID();
    if ( fEventID == 0 )
        return;
    
    if (fEventID != fOldEvent) { // sample time only once per history basis
        G4bool resample = true;
        while ( resample ) {
            fShiftTime = G4RandFlat::shoot(fTimeMean-fTimeFWHM*0.5, fTimeMean+fTimeFWHM*0.5);
            if ( fShiftTime >= 1*ps )
                resample = false;
        }
        
        fShiftTime += fPulseTimeShift; // fPulseTimeShift = 0 if numberOfPulses = 1
        fOldEvent = fEventID;
    }
}


void TsScoreWithIRTPulses::AbsorbResultsFromWorkerScorer(TsVScorer* workerScorer) {
    TsVNtupleScorer::AbsorbResultsFromWorkerScorer(workerScorer);
    
    TsScoreWithIRTPulses* workerGvalueScorer = dynamic_cast<TsScoreWithIRTPulses*>(workerScorer);
    
    // scalar quantities
    fTotalEnergyDeposit += workerGvalueScorer->fTotalEnergyDeposit;
    fEnergyDepositPerEvent += workerGvalueScorer->fEnergyDepositPerEvent;
    fDosePerPulse += workerGvalueScorer->fDosePerPulse;
    
    // G value accumulation
    for (auto& NameAndElse: workerGvalueScorer->fGValuePerSpeciePerTime) {
        G4String MoleculeName = NameAndElse.first;
        for (auto& TimeAndGvalue:NameAndElse.second) {
            G4double MoleculeTime = TimeAndGvalue.first;
            G4double MoleculesGvalue = TimeAndGvalue.second;
            fGValuePerSpeciePerTime[MoleculeName][MoleculeTime]+=MoleculesGvalue;
        }
    }
    
    // Merge the Number Of Molecules Per Time
    for (auto& NameAndElse: workerGvalueScorer->fMoleculesPerSpeciePerTime) {
        G4String MoleculeName = NameAndElse.first;
        for (auto& TimeAndMolecules:NameAndElse.second) {
            G4double MoleculeTime = TimeAndMolecules.first;
            G4double MoleculesNum = TimeAndMolecules.second;
            fMoleculesPerSpeciePerTime[MoleculeName][MoleculeTime]+=MoleculesNum;
        }
    }
    
    // Report Delta: Merge the DeltaG
    for (auto& indexTimeAndDelta : workerGvalueScorer->fDeltaGPerReactionPerTime){
        G4int index = indexTimeAndDelta.first;
        for (auto& timeAndDelta : indexTimeAndDelta.second) {
            G4double time  = timeAndDelta.first;
            G4double delta = timeAndDelta.second;
            fDeltaGPerReactionPerTime[index][time] += delta;
        }
    }
    
    // Energy deposits
    for ( size_t u = 0; u < (workerGvalueScorer->fVEnergyDeposits).size(); u++ ) {
        fVEnergyDeposits[u] += workerGvalueScorer->fVEnergyDeposits[u];
    }
    
    workerGvalueScorer->fGValuePerSpeciePerTime.clear();
    workerGvalueScorer->fDeltaGPerReactionPerTime.clear();
    workerGvalueScorer->fMoleculesPerSpeciePerTime.clear();
    workerGvalueScorer->fVEnergyDeposits.clear();
    workerGvalueScorer->fEnergyDepositPerEvent = 0.0;
    workerGvalueScorer->fDosePerPulse = 0.0;
}


void TsScoreWithIRTPulses::Output() {
    
    std::map<G4String, std::map<G4double, G4double> >::iterator wIter;
    std::map<G4double, G4double>::iterator iter;
    
    for ( wIter = fGValuePerSpeciePerTime.begin(); wIter != fGValuePerSpeciePerTime.end(); wIter++ ) {
        
        for ( iter = (wIter->second).begin(); iter != (wIter->second).end(); iter++ ) {
            fTime = iter->first;
            fMoleculeName = wIter->first;
            
            if (fMoleculeName == "")
                continue;
            
            fGValue    = iter->second;
            fMolecules = fMoleculesPerSpeciePerTime[fMoleculeName][fTime];
            fEnergyDepositGvalue = fVEnergyDeposits[fUtils->FindBin(fTime, fVStepTimes)];
            
            fNtuple->Fill();
        }
    }
    
    fNtuple->Write();
    
    for ( size_t u = 0; u < fVEnergyDepositPerSampledTime.size(); u++ )  {
        G4double saveTime = fVEnergyDepositPerSampledTime[u].first/ps;
        G4double saveEdep = fVEnergyDepositPerSampledTime[u].second/eV;
        fTimeOutFile.write(reinterpret_cast<char*>(&saveTime), sizeof saveTime);
        fTimeOutFile.write(reinterpret_cast<char*>(&saveEdep), sizeof saveEdep);
    }
    fTimeOutFile.close();
    
    std::map<G4int, std::map<G4double, G4double> >::iterator wDeltaIter;
    std::map<G4double, G4double>::iterator deltaiter;
    
    G4int ReactionIndex;
    G4double ReactionTime;
    G4double DeltaReaction;
    G4String ReactA;
    G4String ReactB;
    
    std::ofstream DeltaGFile(fOutFileName + "_DeltaG.phsp");
    DeltaGFile << "# ReactionID moleculeA moleculeB ProdA ProdB ProdC Time[ps] numberOfReactions accumulateEnergyDeposit[eV]" << G4endl;
    for ( wDeltaIter = fDeltaGPerReactionPerTime.begin(); wDeltaIter != fDeltaGPerReactionPerTime.end(); wDeltaIter++ ) {
        for ( deltaiter = (wDeltaIter->second).begin(); deltaiter != (wDeltaIter->second).end(); deltaiter++) {
            DeltaReaction = deltaiter->second;
            ReactionTime  = deltaiter->first;
            ReactionIndex = wDeltaIter->first;
            ReactA = (fIRT->GetReactants(ReactionIndex)).first;
            ReactB = (fIRT->GetReactants(ReactionIndex)).second;
            std::vector<G4String> Products = fIRT->GetProducts(ReactionIndex);
            DeltaGFile << ReactionIndex+1 << "    "  << ReactA << "  "  << ReactB << "  ";
            
            while (Products.size() < 3)
                Products.push_back("None");
            
            for (size_t prod = 0; prod < Products.size(); prod++)
                DeltaGFile << Products[prod] << "  ";
            
            G4int tBin = fUtils->FindBin(ReactionTime, fVStepTimes);
            DeltaGFile << "  " << ReactionTime/ps << "     " << DeltaReaction << "    " << fVEnergyDeposits[tBin]/eV << std::endl;
        }
    }
    DeltaGFile.close();
}


void TsScoreWithIRTPulses::Clear() {
    fGValuePerSpeciePerTime.clear();
    fMoleculesPerSpeciePerTime.clear();
    fDeltaGPerReactionPerTime.clear();
    
    UpdateFileNameForUpcomingRun();
}



