// Scorer for TsIRTStrandBreaksPulsed2
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
#include "TsIRTManager.hh"
#include "TsIRTPlasmidSupercoiled.hh"
#include "TsIRTConfiguration.hh"

#include "TsGeometryManager.hh"


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


TsScoreDNADamagePulsed2::TsScoreDNADamagePulsed2(TsParameterManager* pM, TsMaterialManager* mM, TsGeometryManager* gM, TsScoringManager* scM, TsExtensionManager* eM,
                                                     G4String scorerName, G4String quantity, G4String outFileName, G4bool isSubScorer)
: TsVNtupleScorer(pM, mM, gM, scM, eM, scorerName, quantity, outFileName, isSubScorer),
fPm(pM), fEnergyDepositPerEvent(0), fEnergyDepositPerEventEverywhere(0), fName(scorerName), fOldEvent(-1)
{
    SetUnit("");
    
    if (fPm->GetIntegerParameter("Ts/NumberOfThreads") != 1) {
        G4cerr << "TOPAS is exiting due to an error in parameter." << G4endl;
        G4cerr << "Ts/NumberOfThreads has a value different than 1" << G4endl;
        G4cerr << "This example only can be run in single thread mode" << G4endl;
        G4cerr << "Set: i:Ts/NumberOfThreads = 1" << G4endl;
        fPm->AbortSession(1);
    }
    
    fIRT = new TsIRTManager(fPm, fName);
    fUtils = fIRT->GetUtils();
    fStepTimes = fIRT->GetStepTimes();
    for ( size_t u = 0; u < fStepTimes.size(); u++ ) {
        fVEnergyDeposit.push_back(0.0);
        fVEnergyDepositEverywhere.push_back(0.0);
    }
    
    fNtuple->RegisterColumnD(&fGValue, "TotalNumberOfSpeciesAtTime", "");
    fNtuple->RegisterColumnD(&fTime,    "Time", "picosecond");
    fNtuple->RegisterColumnS(&fMoleculeName, "MoleculeName");
    fNtuple->RegisterColumnD(&fEnergy, "EnergyDeposit", "eV");
    
    fPrescribedDose = fPm->GetDoubleParameter(GetFullParmName("PrescribedDose"),"Dose");
    
    fTimeDistribution = fPm->GetStringParameter(GetFullParmName("PulseDistribution"));
    fTimeDistribution.toLower();
    
    fNumberOfPulses = 1;
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
    
    // To filter only physical interactiones within a virtual TsBox.
    fTestIsInside = false;
    if ( fPm->ParameterExists(GetFullParmName("SensitiveVolumeName"))) {
        fTestIsInside = true;
        fSensitiveVolume = fPm->GetStringParameter(GetFullParmName("SensitiveVolumeName"));
    }
    
    if (fPm->ParameterExists(GetFullParmName("DNAMoleculeNames"))) { // From the DNA/chemistry model's name, these are used in IRT
        G4int NbOfDNANames   = fPm->GetVectorLength(GetFullParmName("DNAMoleculesNames"));
        G4int NbOfDNAWeights = fPm->GetVectorLength(GetFullParmName("DNAMoleculesWeights"));
        
        G4String* DNAMoleculesNames   = fPm->GetStringVector(GetFullParmName("DNAMoleculesNames"));
        G4double* DNAMoleculesWeights = fPm->GetUnitlessVector(GetFullParmName("DNAMoleculesWeights"));
        
        if (NbOfDNANames != NbOfDNAWeights) {
            G4cout << "Please provide as many DNA molecules as Weights!!" << G4endl;
            fPm->AbortSession(1);
        }
        
        for (G4int i = 0; i < NbOfDNANames; i++) {
            fDNAMoleculesToInsert.push_back(DNAMoleculesNames[i]);
            fDNAMoleculesWeights.push_back(DNAMoleculesWeights[i]);
        }
    } else {
        fDNAMoleculesToInsert.push_back("Deoxyribose");
        fDNAMoleculesWeights.push_back(1);
    }
    fOutputFile = "TOPAS";
    if ( fPm->ParameterExists(GetFullParmName("OutputFile"))) {
        fOutputFile = fPm->GetStringParameter(GetFullParmName("OutputFile"));
    }
    
    fComponentName = fPm->GetStringParameter(GetFullParmName("Component"));
    
    fTCut = 1.0 * ps;
    
    fNbOfScoredEvents = 0;
    
    fMass = 0.0;
    fDensity = 1.0 * g/(cm*cm*cm);
    fTotalDose = 0.0;
    fDosePerPulse = 0.0;
    fPulseTimeShift = 0.0;
    
    fShiftTime = 0.0 * ps;
    fMinShiftTime = 1.0 * s;
    
    G4String fileName = fPm->GetStringParameter("Sc/RootFileName") + ".bin";
    remove(fileName);
    fTimeOutFile.open(fileName, std::ios::binary | std::ios::app);
    
    GetDNAInformation();
}


TsScoreDNADamagePulsed2::~TsScoreDNADamagePulsed2()
{
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
    std::map<G4String,G4int> MolIDs = fIRT->GetIRTConfiguration()->GetMoleculeIDs();
    for (size_t i = 0; i < fPHSPMolecules.size(); i++) {
        if (fPHSPIsAlive[i]) { // only insert those that are not damage in the direct action
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
            fIRT->AddMolecule(aMol);
        }
    }
}

G4String TsScoreDNADamagePulsed2::SampleDNAMolecule() {
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


G4bool TsScoreDNADamagePulsed2::ProcessHits(G4Step* aStep, G4TouchableHistory*)
{
    if (!fIsActive) {
        fSkippedWhileInactive++;
        return false;
    }
    
    if ( -1 < aStep->GetTrack()->GetTrackID() ) {
        G4double edep = aStep->GetTotalEnergyDeposit() ;
        if ( edep > 0 ) {
            ResolveSolid(aStep);
            edep *= aStep->GetPreStepPoint()->GetWeight();
            
            if ( 0 == fMass )
                fMass = fDensity * fSolid->GetCubicVolume();

            G4double dose = edep / fMass;
            fTotalDose += dose;
            fEnergyDepositPerEvent += edep ;
            fDosePerPulse += dose;
            
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
            
            return true;
        }
    } else {
        fEventID = GetEventID();
        // Sample pulse's time only once at each new G4event/history.
        if ( fEventID != fOldEvent ) { // New history, sample a time within the pulse
            G4bool resample;
            if ( fTimeDistributionType == 1 ) {
                resample = true;
                while ( resample ) {
                    fShiftTime = G4RandGauss::shoot(fTimeMean+fPulseTimeShift, fTimeStdv);
                    if ( fShiftTime > 0 )
                        resample = false;
                }
            } else if ( fTimeDistributionType == 2 ) {
                resample = true;
                while (resample) {
                    fShiftTime = G4RandFlat::shoot(fTimeMean-fTimeFWHM*0.5, fTimeMean+fTimeFWHM*0.5);
                    if ( fShiftTime > 0 )
                        resample = false;
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
                resample = true;
                while (resample) {
                    fShiftTime = G4RandExponential::shoot(fTimeMean + fPulseTimeShift);
                    if ( fShiftTime >= 1.0*ps )
                        resample = false;
                }
            } else {
                fShiftTime = 1.0 * ps;
            }
            
            fOldEvent = fEventID;
            // Force the first event arrival at 1 ps.
            if ((fNbOfScoredEvents == 0 && fNbOfScoredEventsEverywhere == 0))
                fShiftTime = 1.0 * ps;
        }
        
        return true;
    }
    return false;
}


void TsScoreDNADamagePulsed2::UserHookForPreTimeStepAction() {
    if (!G4EventManager::GetEventManager()->GetConstCurrentEvent()->IsAborted()) {
        G4TrackManyList* trackList = G4ITTrackHolder::Instance()->GetMainList();
        G4ManyFastLists<G4Track>::iterator it_begin = trackList->begin();
        G4ManyFastLists<G4Track>::iterator it_end = trackList->end();
        
        for(;it_begin!=it_end;++it_begin){
            if ( fTestIsInside ) {
                const G4String& volumeName = (*it_begin)->GetVolume()->GetName();
                if ( G4StrUtil::contains(volumeName,fSensitiveVolume) ) {
                    fIRT->AddMolecule(*it_begin, fShiftTime, 0, G4ThreeVector());
                }
            }
            else {
                fIRT->AddMolecule(*it_begin, fShiftTime, 0, G4ThreeVector());
            }
        }
        
        G4Scheduler::Instance()->Stop();
    }
}


void TsScoreDNADamagePulsed2::UserHookForEndOfEvent() {
    if ( fEnergyDepositPerEvent > 0 ) {
        fVEnergyDepositPerEvent.push_back(std::make_pair(fShiftTime,fEnergyDepositPerEvent));
        fNbOfScoredEvents++;
    }
    
    fEnergyDepositPerEvent = 0.0;
    
    if(fNumberOfPulses > 1 && fDosePerPulse >= fPrescribedDose/fNumberOfPulses) {
        fDosePerPulse = 0.0;
        fPulseTimeShift += fPulsesTimeDelay;
        G4cout << "-- New Pulse at " << fPulseTimeShift/ps << " ps " << G4endl;
    }
    
    if(fTotalDose > fPrescribedDose )  { // Add everything to IRT
        G4int tBin;
        G4cout << " --- IRT start for event " << GetEventID() << G4endl;
        
        for ( size_t t = 0; t < fVEnergyDepositPerEvent.size(); t++ ) {
            tBin = fUtils->FindBin(fVEnergyDepositPerEvent[t].first, fStepTimes);
            
            for ( int i = tBin; i < (int)fStepTimes.size(); i++)
                fVEnergyDeposit[i] += fVEnergyDepositPerEvent[t].second;
            
            G4double saveTime = fVEnergyDepositPerEvent[t].first/ps;
            G4double saveEdep = fVEnergyDepositPerEvent[t].second/joule;
            fTimeOutFile.write(reinterpret_cast<char*>(&saveTime), sizeof saveTime);
            fTimeOutFile.write(reinterpret_cast<char*>(&saveEdep), sizeof saveEdep);
            
            G4cout << " --- energy deposit at time " << fVEnergyDepositPerEvent[t].first/ps
            << " ps " << fVEnergyDepositPerEvent[t].second/eV << G4endl;
        }
        fTimeOutFile.close();
        
        InsertDNAMolecules();
        fIRT->runIRT();
        
        std::map<G4String, std::map<G4double, G4int>> irt = fIRT->GetGValues();
        
        G4cout << " --- IRT ends." << G4endl;
        for ( auto& nameTimeAndGvalue : irt ) {
            G4String name = nameTimeAndGvalue.first;
            fMoleculeName = name;
            for ( auto& timeAndGvalue : (nameTimeAndGvalue.second) ) {
                G4double time = timeAndGvalue.first;
                G4double gvalue = timeAndGvalue.second;
                fGValue = gvalue;
                fTime = time;
                tBin = fUtils->FindBin(time, fStepTimes);
                fEnergy = fVEnergyDeposit[tBin];
                fNtuple->Fill();
            }
        }
        
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
        
        G4String OutputFileName = fPm->GetStringParameter(GetFullParmName("OutputFile"));
        DeltaGFile.open(fOutputFile + "_DeltaG.phsp", std::ofstream::app);
        DeltaGFile << "# reactionID moleculeA moleculeB prod1 prod2 prod3 Time(ps) DeltaG(/100eV) " << G4endl;
        std::map<G4int, std::map<G4double, G4int>> DeltaG = fIRT->GetDeltaGValues();
        
        for ( auto& indexTimeAndDeltaG : DeltaG) {
            G4int reactionIndex = indexTimeAndDeltaG.first;
            for ( auto& timeAndDeltaG : (indexTimeAndDeltaG.second) ) {
                G4double time   = timeAndDeltaG.first;
                tBin = fUtils->FindBin(time, fStepTimes);
                G4double deltaG = 100 * timeAndDeltaG.second / (fVEnergyDeposit[tBin]/eV);
                G4String ReactA = (fIRT->GetReactants(reactionIndex)).first;
                G4String ReactB = (fIRT->GetReactants(reactionIndex)).second;
                std::vector<G4String> Products = fIRT->GetProducts(reactionIndex);
                DeltaGFile << reactionIndex+1 << "    "  << ReactA << "  "  << ReactB << "  ";
                
                while (Products.size() < 3)
                    Products.push_back("None");
                
                for (size_t prod = 0; prod < Products.size(); prod++)
                    DeltaGFile << Products[prod] << "  ";
                
                DeltaGFile << "  " << time/ps << "     " << deltaG << std::endl;
            }
        }
        DeltaGFile.close();
        DeltaG.clear();
        
        irt.clear();
        fIRT->Clean();
        
        G4cout << " ------ Needed " << fNbOfScoredEvents
        << " to achieve an accumulated dose of " << fTotalDose/gray
        << " Gy" << G4endl;
        
        fNtuple->Write();
        
        Output();
        Clear();
        G4RunManager::GetRunManager()->AbortRun(true);
    }
}


