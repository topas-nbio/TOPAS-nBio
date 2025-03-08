// Scorer for TsDNADamageWithIRTPulses
//
// ********************************************************************
// *                                                                  *
// * This file is part of the TOPAS-nBio extensions to the            *
// *   TOPAS Simulation Toolkit.                                      *
// * The TOPAS-nBio extensions are freely available under the license *
// *   agreement set forth at: https://topas-nbio.readthedocs.io/     *
// *                                                                  *
// ********************************************************************

#include "TsScoreDNADamagePulsed2.hh"
#include "TsIRTPlasmidSupercoiled.hh"
#include "TsGeometryManager.hh"
#include "TsVGeometryComponent.hh"

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

TsScoreDNADamagePulsed2::TsScoreDNADamagePulsed2(TsParameterManager* pM, TsMaterialManager* mM, TsGeometryManager* gM, TsScoringManager* scM, TsExtensionManager* eM,
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
    
    for ( size_t u = 0; u < fVStepTimes.size(); u++ )
        fVEnergyDeposits.push_back(0.0);
    
    // For DNA information
    fComponentName = fPm->GetStringParameter(GetFullParmName("Component"));
    fStrandBreakMolecules = fPm->GetStringVector(GetFullParmName("StrandBreakMolecules"));
    fNumberOfStrandBreakMolecules = fPm->GetVectorLength(GetFullParmName("StrandBreakMolecules"));
    
    fDNAMoleculesToInsert = fPm->GetStringVector(GetFullParmName("DNAMoleculesNames"));
    fDNAMoleculesWeights = fPm->GetUnitlessVector(GetFullParmName("DNAMoleculesWeights"));
    fNumberOfMoleculesToInsert = fPm->GetVectorLength(GetFullParmName("DNAMoleculesWeights"));
    
    fDNAHasBeenInserted = false;
    GetDNAInformation();
}


TsScoreDNADamagePulsed2::~TsScoreDNADamagePulsed2()
{
    delete fIRT;
}


G4bool TsScoreDNADamagePulsed2::ProcessHits(G4Step* aStep, G4TouchableHistory*)
{
    if (!fIsActive) {
        fSkippedWhileInactive++;
        return false;
    }
    
    if ( -1 < aStep->GetTrack()->GetTrackID() ) {
        G4double edep = aStep->GetTotalEnergyDeposit();
        if ( edep > 0 ) {
            ResolveSolid(aStep);
            
            // If used, it test if the hits fell within a virtual box.
            // It assumes all density is unity, thus gets the mass from rho x cubic_volume
            if (fIRT->GetIRTProcedure()->GetTestIsInside()) {
                if ( 0 == fMass )
                    fMass = fDensity * fIRT->GetIRTProcedure()->GetVirtualRegionCubicVolume();
                
                if (fIRT->GetIRTProcedure()->Inside(aStep->GetPreStepPoint()->GetPosition())) {
                    fEnergyDepositPerEvent += edep * aStep->GetPreStepPoint()->GetWeight();
                    ScoreEnergyDepositInDNA(aStep);
                    return true;
                } else {
                    return false;
                }
                
            } else {
                if ( 0 == fMass )
                    fMass = fDensity * fSolid->GetCubicVolume();
                
                fEnergyDepositPerEvent += edep * aStep->GetPreStepPoint()->GetWeight();
                ScoreEnergyDepositInDNA(aStep);
                return true;
            }
        }
    }
    
    return false;
}

void TsScoreDNADamagePulsed2::ScoreEnergyDepositInDNA(G4Step* aStep) {
    G4TouchableHistory* touchable = (G4TouchableHistory*)(aStep->GetPreStepPoint()->GetTouchable());
    const G4String& vName = touchable->GetVolume()->GetName();
    
    G4double edep = aStep->GetTotalEnergyDeposit() * aStep->GetPreStepPoint()->GetWeight();
    // Comparing strings is expensive, but the chemistry takes much longer in comparison.
    
    if (G4StrUtil::contains(vName,"deoxyribose") || G4StrUtil::contains(vName,"phosphate")) {

        G4int plasmidID, strandID;
        G4int baseID = touchable->GetVolume(0)->GetCopyNo();
        
        if (G4StrUtil::contains(vName,"water")) {plasmidID = touchable->GetVolume(1)->GetCopyNo();}
        else {plasmidID = touchable->GetVolume(2)->GetCopyNo();}
        
        if (G4StrUtil::contains(vName,"1")) {strandID = 1;}
        else {strandID = 2;}
        
        if (fDNAEnergyDeposit.find(plasmidID) == fDNAEnergyDeposit.end()) {
            fDNAEnergyDeposit[plasmidID][baseID][strandID] = edep;
        } else {
            if (fDNAEnergyDeposit[plasmidID].find(baseID) == fDNAEnergyDeposit[plasmidID].end()) {
                fDNAEnergyDeposit[plasmidID][baseID][strandID] = edep;
            } else {
                if (fDNAEnergyDeposit[plasmidID][baseID].find(strandID) == fDNAEnergyDeposit[plasmidID][baseID].end()) {
                    fDNAEnergyDeposit[plasmidID][baseID][strandID] = edep;
                } else {
                    fDNAEnergyDeposit[plasmidID][baseID][strandID] += edep;
                }
            }
        }
    } else {
        return;
    }
    
    return;
}

void TsScoreDNADamagePulsed2::UserHookForEndOfEvent() {
    if(fEnergyDepositPerEvent == 0) return;

    if(fEnergyDepositPerEvent > 0) {
        SampleTimeShift();
        fVEnergyDepositPerSampledTime.push_back(std::make_pair(fShiftTime,fEnergyDepositPerEvent));
        fTotalEnergyDeposit += fEnergyDepositPerEvent;
        fDosePerPulse += fEnergyDepositPerEvent/fMass;
        fEnergyDepositPerEvent = 0;
        fEventID++;
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
        << " ps achieved a Total Dose " << (fTotalEnergyDeposit/fMass)/gray << " Gy with " << fEventID << " histories. " << G4endl;
        
        for ( size_t t = 0; t < fVEnergyDepositPerSampledTime.size(); t++ ) {
            G4int tBin = fUtils->FindBin(fVEnergyDepositPerSampledTime[t].first, fVStepTimes);
            
            for ( int i = tBin; i < (int)fVStepTimes.size(); i++)
                fVEnergyDeposits[i] += fVEnergyDepositPerSampledTime[t].second;
        }
        
        // Total energy deposit per history
        std::ofstream timeOutFile(fOutFileName + ".bin", std::ios::binary);
        for ( size_t u = 0; u < fVEnergyDepositPerSampledTime.size(); u++ )  {
            G4double saveTime = fVEnergyDepositPerSampledTime[u].first/ps;
            G4double saveEdep = fVEnergyDepositPerSampledTime[u].second/eV;
            timeOutFile.write(reinterpret_cast<char*>(&saveTime), sizeof saveTime);
            timeOutFile.write(reinterpret_cast<char*>(&saveEdep), sizeof saveEdep);
        }
        timeOutFile.close();
        
        //------------------------------------------------
        CheckForDirectDNABreaks();
        InsertDNAMolecules();
        fIRT->runIRT();
        CheckForIndirectDNABreaks();
        //------------------------------------------------
        
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
        
        //Output();

        DeltaG.clear();
        irt.clear();
        fIRT->Clean();
        
        G4RunManager::GetRunManager()->AbortRun(true);
    }
    
    return;
}


void TsScoreDNADamagePulsed2::UserHookForPreTimeStepAction() {
    if(fEnergyDepositPerEvent == 0) return;
    
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
    for(;it_begin!=it_end;++it_begin)
        fIRT->AddMolecule(*it_begin, fShiftTime, 0, G4ThreeVector());
    
    G4Scheduler::Instance()->Stop();
}


void TsScoreDNADamagePulsed2::SampleTimeShift() {
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


void TsScoreDNADamagePulsed2::AbsorbResultsFromWorkerScorer(TsVScorer* workerScorer) {
    TsVNtupleScorer::AbsorbResultsFromWorkerScorer(workerScorer);
    
    TsScoreDNADamagePulsed2* workerGvalueScorer = dynamic_cast<TsScoreDNADamagePulsed2*>(workerScorer);
    
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
    for ( size_t u = 0; u < (workerGvalueScorer->fVEnergyDeposits).size(); u++ )
        fVEnergyDeposits[u] += workerGvalueScorer->fVEnergyDeposits[u];
    
    for ( size_t u = 0; u < (workerGvalueScorer->fDirectStrandBreaks).size(); u++ )
        fDirectStrandBreaks.push_back(workerGvalueScorer->fDirectStrandBreaks[u]);
        
    for ( size_t u = 0; u < (workerGvalueScorer->fIndirectStrandBreaks).size(); u++ )
        fIndirectStrandBreaks.push_back(workerGvalueScorer->fIndirectStrandBreaks[u]);
            
    workerGvalueScorer->fGValuePerSpeciePerTime.clear();
    workerGvalueScorer->fDeltaGPerReactionPerTime.clear();
    workerGvalueScorer->fMoleculesPerSpeciePerTime.clear();
    workerGvalueScorer->fVEnergyDeposits.clear();
    workerGvalueScorer->fIndirectStrandBreaks.clear();
    workerGvalueScorer->fDirectStrandBreaks.clear();
    workerGvalueScorer->fEnergyDepositPerEvent = 0.0;
    workerGvalueScorer->fDosePerPulse = 0.0;
    return;
}


void TsScoreDNADamagePulsed2::Output() {
    // G values
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
    
    // Delta G
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
    
    // DNA damage
    std::ofstream StrandBreakOut(fOutFileName + "_StrandBreak.phsp");
    
    StrandBreakOut << "# TOPAS-nBio IRT: DNA Strand Break Map File"   << std::endl;
    StrandBreakOut << "# BreakID = 1 if Direct SB, 2 if indirect SB from IRT, 3 if indirect SB from Gillespie"     << std::endl;
    StrandBreakOut << "# plasmidID | BaseID | StrandID | MoleculeID | BreakID |" << std::endl;
    
    for (size_t j = 0; j  < fDirectStrandBreaks.size(); j++) {
        StrandBreakOut
        << std::setw(11) << fDirectStrandBreaks[j][0]
        << std::setw(9)  << fDirectStrandBreaks[j][1]
        << std::setw(11) << fDirectStrandBreaks[j][2]
        << std::setw(30) << 0
        << std::setw(10) << 1 << std::endl;
    }
    
    for (size_t j = 0; j  < fIndirectStrandBreaks.size(); j++) {
        int algoCheck = 2;
        if (fIndirectStrandBreaks[j][4] == 0) { // Checking the algorithm used
            algoCheck = 3;
        }
        
        StrandBreakOut
        << std::setw(11) << fIndirectStrandBreaks[j][0]
        << std::setw(9)  << fIndirectStrandBreaks[j][1]
        << std::setw(11) << fIndirectStrandBreaks[j][2]
        << std::setw(30) << fIRT->GetIRTConfiguration()->GetMoleculeNameFromMoleculeID(fIndirectStrandBreaks[j][3])
        << std::setw(10) << algoCheck << std::endl;
    }
            
    StrandBreakOut.close();
}


void TsScoreDNADamagePulsed2::GetDNAInformation() {
    G4String Component = fComponentName;
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

void TsScoreDNADamagePulsed2::InsertDNAMolecules() {
    if (fDNAHasBeenInserted)
        return;
    
    std::map<G4String,G4int> molIDs = fIRT->GetIRTConfiguration()->GetMoleculeIDs();
    G4int accepted = 0;
    for (size_t i = 0; i < fPHSPMolecules.size(); i++) {
        if (fPHSPIsAlive[i]) {
            G4String molName = GetSampledDNAMolecule();
            G4int    molID   = molIDs[molName];
            
            TsIRTConfiguration::TsMolecule aMol;
            aMol.id       = molID;
            aMol.position = fPHSPPosition[i];
            aMol.time     = 1 * ps;
            aMol.spin     = 0;
            aMol.trackID  = 100;
            aMol.isDNA    = true;
            aMol.isNew    = true;
            aMol.volumeID = fDNADetails[i][0];
            aMol.baseID   = fDNADetails[i][1];
            aMol.strandID = fDNADetails[i][2];
            
            fIRT->AddMolecule(aMol);
            accepted++;
        } else {
            fPHSPIsAlive[i] = true;
        }
    }
    
    G4cout << " -- DNA molecules inserted " << accepted << G4endl;
    
    fDNAHasBeenInserted = true;
}

G4String TsScoreDNADamagePulsed2::GetSampledDNAMolecule() {
    G4double random = G4UniformRand();
    G4double sum    = 0;
    
    for (int i = 0; i < fNumberOfMoleculesToInsert; i++) {
        sum += fDNAMoleculesWeights[i];
        if (random <= sum) {
            return fDNAMoleculesToInsert[i];
        }
    }
    
    return fDNAMoleculesToInsert[fNumberOfMoleculesToInsert-1];
}

void TsScoreDNADamagePulsed2::CheckForDirectDNABreaks() {
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
    return;
}


void TsScoreDNADamagePulsed2::CheckForIndirectDNABreaks() {
    // Score the indirect strand breaks
    for (int i = 0; i < fNumberOfStrandBreakMolecules; i++) {
        G4int BreakMolID = (fIRT->GetIRTConfiguration()->GetMoleculeIDs())[fStrandBreakMolecules[i]];
        // Retrieves the products of DNA reactions, those that did survive further reactions: OHDeoxyribose, EAQDeoxyribose, HDeoxyribose
        std::vector<TsIRTConfiguration::TsMolecule> SurvivingMolecules = fIRT->GetSurvivingMoleculesWithMolID(BreakMolID);
        G4cout << " --- Retrieved " << SurvivingMolecules.size() << " from molecule type " << BreakMolID << G4endl;
        for (size_t j = 0; j < SurvivingMolecules.size(); j++) {
            std::vector<G4int> thisBreak;
            thisBreak.push_back(SurvivingMolecules[j].volumeID);
            thisBreak.push_back(SurvivingMolecules[j].baseID);
            thisBreak.push_back(SurvivingMolecules[j].strandID);
            thisBreak.push_back(BreakMolID);
            thisBreak.push_back(SurvivingMolecules[j].chemAlgo);
    
            fIndirectStrandBreaks.push_back(thisBreak);
        }
    }
    
    for (size_t i = 0; i < fDirectStrandBreaksInEvent.size(); i++)
        fDirectStrandBreaks.push_back(fDirectStrandBreaksInEvent[i]);
    
    if (fDirectStrandBreaksInEvent.size() != 0)
        fDirectStrandBreaksInEvent.clear();
}


void TsScoreDNADamagePulsed2::Clear() {
    fGValuePerSpeciePerTime.clear();
    fMoleculesPerSpeciePerTime.clear();
    fDeltaGPerReactionPerTime.clear();
    
    UpdateFileNameForUpcomingRun();
}



