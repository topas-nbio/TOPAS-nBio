// Scorer for TsIRTGvalue
//
// ********************************************************************
// *                                                                  *
// * This file is part of the TOPAS-nBio extensions to the            *
// *   TOPAS Simulation Toolkit.                                      *
// * The TOPAS-nBio extensions are freely available under the license *
// *   agreement set forth at: https://topas-nbio.readthedocs.io/     *
// *                                                                  *
// ********************************************************************

#include "TsScoreWithIRT.hh"
#include "TsIRT.hh"
#include "TsIRTConfiguration.hh"

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
#include "G4Timer.hh"
#include "G4LogicalVolumeStore.hh"

TsScoreWithIRT::TsScoreWithIRT(TsParameterManager* pM, TsMaterialManager* mM, TsGeometryManager* gM, TsScoringManager* scM, TsExtensionManager* eM,
							   G4String scorerName, G4String quantity, G4String outFileName, G4bool isSubScorer)
: TsVNtupleScorer(pM, mM, gM, scM, eM, scorerName, quantity, outFileName, isSubScorer),
fPm(pM), fEnergyDepositPerEvent(0), fEnergyLossKill(0), fName(scorerName)
{
	SetUnit("");
	
	fNtuple->RegisterColumnD(&fGValue, "GValue: number of molecules per 100 eV of energy deposit", "");
	fNtuple->RegisterColumnD(&fGValueError, "GValue statistical error", "");

	if (fPm->ParameterExists(GetFullParmName("ReportMoleculeYield")) &&
		fPm->GetBooleanParameter(GetFullParmName("ReportMoleculeYield"))){
	    fNtuple->RegisterColumnD(&fMolecules,"Number of molecules","");
	    fNtuple->RegisterColumnD(&fMoleculesError,"Number of molecules error","");		
	}

	fNtuple->RegisterColumnD(&fTime,    "Time", "picosecond");
	fNtuple->RegisterColumnS(&fMoleculeName, "MoleculeName");
	
	fEnergyLossKill  = fPm->GetDoubleParameter(GetFullParmName("KillPrimaryIfEnergyLossExceeds"),"Energy");
	fEnergyLossAbort = fPm->GetDoubleParameter(GetFullParmName("AbortEventIfPrimaryEnergyLossExceeds"),"Energy");
	fMaximumTrackLength = 1*cm;
	if (fPm->ParameterExists(GetFullParmName("KillPrimaryBasedOnTrackLength")))
		fMaximumTrackLength = fPm->GetDoubleParameter(GetFullParmName("KillPrimaryBasedOnTrackLength"), "Length");

	G4double TMiddle = 10 * us;
	if ( fPm->ParameterExists(GetFullParmName("IRTTimeEnd")) )
		TMiddle = fPm->GetDoubleParameter(GetFullParmName("IRTTimeEnd"), "Time");
	
	fNumberOfSteps = 0;
	fMaximumNumberOfSteps = 0;
	fUseMaximumNumberOfSteps = false;
	if ( fPm->ParameterExists(GetFullParmName("MaximumNumberOfSteps")) ) {
		fMaximumNumberOfSteps = fPm->GetIntegerParameter(GetFullParmName("MaximumNumberOfSteps"));
		fUseMaximumNumberOfSteps = true;
	}

	fReportDelta = false;
	if ( fPm->ParameterExists(GetFullParmName("ReportDeltaGValues"))) {
		fReportDelta = fPm->GetBooleanParameter(GetFullParmName("ReportDeltaGValues"));
	}

	fOutputFile = "TOPAS";
	if ( fPm->ParameterExists(GetFullParmName("OutputFile"))) {
		fOutputFile = fPm->GetStringParameter(GetFullParmName("OutputFile"));
	}

	fVolume = 0;

	fIRT       = new TsIRT(fPm, fName);
	
	fNbOfScoredEvents = 0;
	fEnergyLoss = 0.0;
	fTotalTrackLength = 0.0;
	fEnergyDepositPlusEnergyKinetic = 0.0;
	fEnergyDepositPerEvent = 0.0;
	fIRTExecutionTime = 0.0;
	fIRTExecutionTimeStdv = 0.0;


	G4String DeltaGFile = fOutputFile+"_DeltaG.phsp";
	const char* DeltaGFileOutput = DeltaGFile.c_str();
	remove(DeltaGFileOutput);
}


TsScoreWithIRT::~TsScoreWithIRT()
{
	delete fIRT;
}


G4bool TsScoreWithIRT::ProcessHits(G4Step* aStep, G4TouchableHistory*)
{
	if (!fIsActive) {
		fSkippedWhileInactive++;
		return false;
	}
	
	if ( -1 < aStep->GetTrack()->GetTrackID() ) {
		
		if ( aStep->GetTrack()->GetParentID() == 0 ) {
			ResolveSolid(aStep);
			fVolume = fSolid->GetCubicVolume();
			
			if ( !fUseMaximumNumberOfSteps ) {
				fTotalTrackLength += aStep->GetStepLength();
				
				G4double kineticE = aStep->GetPostStepPoint()->GetKineticEnergy();
				G4double eLoss = aStep->GetPreStepPoint()->GetKineticEnergy() - kineticE;
				fEnergyLoss += eLoss;
				
				fLET += aStep->GetTotalEnergyDeposit();
				
				fEnergyDepositPlusEnergyKinetic += aStep->GetTotalEnergyDeposit();
				
				if ( fEnergyLoss > fEnergyLossAbort ) {
					G4RunManager::GetRunManager()->AbortEvent();
					return true;
				}
				if ( fEnergyLoss >= fEnergyLossKill ) {
					aStep->GetTrack()->SetTrackStatus(fStopAndKill);
					return true;
				}
				if ( fTotalTrackLength >= fMaximumTrackLength ) {
					aStep->GetTrack()->SetTrackStatus(fStopAndKill);
					return true;
				}
			} else {
				fNumberOfSteps++;
				if ( fMaximumNumberOfSteps < fNumberOfSteps ) {
					aStep->GetTrack()->SetTrackStatus(fStopAndKill);
				}
			}
		} else if ( aStep->GetTrack()->GetParentID() == 1 && aStep->GetTrack()->GetCurrentStepNumber() == 1 ) {
			G4double ekin = aStep->GetPreStepPoint()->GetKineticEnergy();
			if ( ekin < 100.0 * eV ) {
				fEnergyDepositPlusEnergyKinetic += ekin;
				fLET += ekin;
			}
		}
		
		G4double edep = aStep->GetTotalEnergyDeposit();
		if ( edep > 0 ) {
			fEnergyDepositPerEvent += edep * aStep->GetPreStepPoint()->GetWeight();
			return true;
		}
	} 

	return false;
}


void TsScoreWithIRT::UserHookForEndOfEvent() {
	if(!G4EventManager::GetEventManager()->GetConstCurrentEvent()->IsAborted() && fEnergyDepositPerEvent > 0) {

		fIRTTimer[0].Start();
		fIRT->runIRT();
		fIRTTimer[0].Stop();
			
		std::map<G4String, std::map<G4double, G4int>> irt = fIRT->GetGValues();
		fIRTExecutionTime     += fIRTTimer[0].GetUserElapsed();
		fIRTExecutionTimeStdv += (fIRTTimer[0].GetUserElapsed())*(fIRTTimer[0].GetUserElapsed());

		for ( auto& nameTimeAndGvalue : irt ) {
			G4String name = nameTimeAndGvalue.first;
			for ( auto& timeAndGvalue : (nameTimeAndGvalue.second) ) {
				G4double time = timeAndGvalue.first;
				G4double gvalue = timeAndGvalue.second;

			    fMoleculesPerSpeciePerTime[name][time] += gvalue;
			    fMoleculesPerSpeciePerTime2[name][time] += gvalue*gvalue;

				gvalue *= 100/(fEnergyDepositPerEvent/eV);
				
				fGValuePerSpeciePerTime[name][time] += gvalue;
				fGValuePerSpeciePerTime2[name][time] += gvalue*gvalue;
			}
		}

		if (fReportDelta) {
			std::map<G4int, std::map<G4double, G4int>> DeltaG = fIRT->GetDeltaGValues();
			for ( auto& indexTimeAndDeltaG : DeltaG) {
				G4int reactionIndex = indexTimeAndDeltaG.first;
				for ( auto& timeAndDeltaG : (indexTimeAndDeltaG.second) ) {
					G4double time   = timeAndDeltaG.first;
					G4double deltaG = timeAndDeltaG.second;
					deltaG *= 100/(fEnergyDepositPerEvent/eV);

					fDeltaGPerReactionPerTime[reactionIndex][time] += deltaG;
					fDeltaGPerReactionPerTime2[reactionIndex][time] += deltaG*deltaG;
				}
			}
			DeltaG.clear();
		}
			
		fEnergyDepositPlusEnergyKinetic /= fTotalTrackLength;
		fLET /= fTotalTrackLength;
		
		irt.clear();
		fIRT->Clean();
		fNbOfScoredEvents++;
	}
	
	fEnergyDepositPerEvent = 0.0;
	fEnergyLoss = 0.0;
	fTotalTrackLength = 0.0;
	fNumberOfSteps = 0;
	fEnergyDepositPlusEnergyKinetic = 0.0;
	fLET = 0.0;
	return;
}


void TsScoreWithIRT::UserHookForPreTimeStepAction() {
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
		for(;it_begin!=it_end;++it_begin){
			G4double time = it_begin->GetGlobalTime();
	     	fIRT->AddMolecule(*it_begin, time, 0, G4ThreeVector());
		}

		G4Scheduler::Instance()->Stop();
	}
}


void TsScoreWithIRT::AbsorbResultsFromWorkerScorer(TsVScorer* workerScorer) {
	TsVNtupleScorer::AbsorbResultsFromWorkerScorer(workerScorer);
	
	TsScoreWithIRT* workerGvalueScorer = dynamic_cast<TsScoreWithIRT*>(workerScorer);
	
	fNbOfScoredEvents += workerGvalueScorer->fNbOfScoredEvents;
	fIRTExecutionTime += workerGvalueScorer->fIRTExecutionTime;
	fIRTExecutionTimeStdv += workerGvalueScorer->fIRTExecutionTimeStdv;
	
	std::map<G4String, std::map<G4double, G4double> >::iterator wIter;
	std::map<G4String, std::map<G4double, G4double> >::iterator mIter;
	
	for ( wIter = workerGvalueScorer->fGValuePerSpeciePerTime.begin(); wIter != workerGvalueScorer->fGValuePerSpeciePerTime.end(); wIter++) {
		mIter = fGValuePerSpeciePerTime.find( wIter->first );
		if ( mIter == fGValuePerSpeciePerTime.end() ) {
			fGValuePerSpeciePerTime.insert(std::pair<G4String, std::map<G4double, G4double> > ( wIter->first, wIter->second));
		} else {
			std::map<G4double, G4double>::iterator witer;
			std::map<G4double, G4double>::iterator miter;
			for ( witer = (wIter->second).begin(); witer != (wIter->second).end(); witer++ ) {
				miter = (mIter->second).find( witer->first );
				miter->second += witer->second;
			}
		}
	}
	
	std::map<G4String, std::map<G4double, G4double> >::iterator wIter2;
	std::map<G4String, std::map<G4double, G4double> >::iterator mIter2;
	
	for ( wIter2 = workerGvalueScorer->fGValuePerSpeciePerTime2.begin(); wIter2 != workerGvalueScorer->fGValuePerSpeciePerTime2.end(); wIter2++) {
		mIter2 = fGValuePerSpeciePerTime2.find( wIter2->first );
		if ( mIter2 == fGValuePerSpeciePerTime2.end() ) {
			fGValuePerSpeciePerTime2.insert(std::pair<G4String, std::map<G4double, G4double> > ( wIter2->first, wIter2->second));
		} else {
			std::map<G4double, G4double>::iterator witer2;
			std::map<G4double, G4double>::iterator miter2;
			for ( witer2 = (wIter2->second).begin(); witer2 != (wIter2->second).end(); witer2++ ) {
				miter2 = (mIter2->second).find( witer2->first );
				miter2->second += witer2->second;
			}
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

	// Merge the Number Of Molecules squared Per Time
	for (auto& NameAndElse: workerGvalueScorer->fMoleculesPerSpeciePerTime2) {
		G4String MoleculeName = NameAndElse.first;
		for (auto& TimeAndMolecules:NameAndElse.second) {
			G4double MoleculeTime = TimeAndMolecules.first;
			G4double MoleculesNum = TimeAndMolecules.second;
			fMoleculesPerSpeciePerTime2[MoleculeName][MoleculeTime]+=MoleculesNum;
		}
	}

	// If Report Delta: Merge the DeltaG
	if (fReportDelta) {
		for (auto& indexTimeAndDelta : workerGvalueScorer->fDeltaGPerReactionPerTime){
			G4int index = indexTimeAndDelta.first;
			for (auto& timeAndDelta : indexTimeAndDelta.second) {
				G4double time  = timeAndDelta.first;
				G4double delta = timeAndDelta.second;
				fDeltaGPerReactionPerTime[index][time] += delta;
			}
		}

		for (auto& indexTimeAndDelta : workerGvalueScorer->fDeltaGPerReactionPerTime2){
			G4int index = indexTimeAndDelta.first;
			for (auto& timeAndDelta : indexTimeAndDelta.second) {
				G4double time  = timeAndDelta.first;
				G4double delta = timeAndDelta.second;
				fDeltaGPerReactionPerTime2[index][time] += delta;
			}
		}
	}

	workerGvalueScorer->fGValuePerSpeciePerTime.clear();
	workerGvalueScorer->fGValuePerSpeciePerTime2.clear();
	workerGvalueScorer->fMoleculesPerSpeciePerTime.clear();
	workerGvalueScorer->fMoleculesPerSpeciePerTime2.clear();
	workerGvalueScorer->fNbOfScoredEvents = 0;
	workerGvalueScorer->fEnergyDepositPerEvent = 0.0;
    workerGvalueScorer->fLET = 0.0;
    workerGvalueScorer->fEnergyDepositPlusEnergyKinetic = 0.0;
    workerGvalueScorer->fNumberOfSteps = 0;
    workerGvalueScorer->fEnergyLoss = 0;
    workerGvalueScorer->fTotalTrackLength = 0;
}


void TsScoreWithIRT::Output() {
	if ( fNbOfScoredEvents == 0 )
		return;
	
	std::map<G4String, std::map<G4double, G4double> >::iterator wIter;
	std::map<G4String, std::map<G4double, G4double> >::iterator wIter2;
	std::map<G4double, G4double>::iterator iter;
	std::map<G4double, G4double>::iterator iter2;
	
	for ( wIter = fGValuePerSpeciePerTime.begin(); wIter != fGValuePerSpeciePerTime.end(); wIter++ ) {
		wIter2 = fGValuePerSpeciePerTime2.find(wIter->first);
		
		for ( iter = (wIter->second).begin(); iter != (wIter->second).end(); iter++ ) {
			fTime = iter->first;
			fMoleculeName = wIter->first;

			if (fMoleculeName == "") 
				continue;

			iter2 = (wIter2->second).find( iter->first );
			
			fGValue    = iter->second/fNbOfScoredEvents;
			fMolecules = fMoleculesPerSpeciePerTime[fMoleculeName][fTime]/fNbOfScoredEvents;
			if ( fNbOfScoredEvents > 1 ) {
				fGValueError = sqrt( (1.0/(fNbOfScoredEvents-1)) * ( (iter2->second)/fNbOfScoredEvents - fGValue*fGValue));
				fMoleculesError = sqrt( (1.0/(fNbOfScoredEvents-1)) * ( (fMoleculesPerSpeciePerTime2[fMoleculeName][fTime])/fNbOfScoredEvents - fMolecules*fMolecules));
			} else {
				fGValueError = 1.0;
				fMoleculesError = 1.0;
			}
			fNtuple->Fill();
		}
	}
	
	fNtuple->Write();
	fIRTExecutionTime /= fNbOfScoredEvents;
	if ( fNbOfScoredEvents == 1 )
		fIRTExecutionTimeStdv = 0.0;
	else
		fIRTExecutionTimeStdv = std::sqrt(fNbOfScoredEvents/(fNbOfScoredEvents-1) * 
                                         (fIRTExecutionTimeStdv/fNbOfScoredEvents-fIRTExecutionTime*fIRTExecutionTime));

	G4cout << "######### Averaged IRT-execution-time per history (" << fNbOfScoredEvents<< " histories): " << fIRTExecutionTime << " s +/- " << fIRTExecutionTimeStdv << " s " << G4endl;

	if (fReportDelta) {
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

		for ( wDeltaIter = fDeltaGPerReactionPerTime.begin(); wDeltaIter != fDeltaGPerReactionPerTime.end(); wDeltaIter++ ) {
			wDeltaIter2 = fDeltaGPerReactionPerTime.find(wDeltaIter->first );

			for ( deltaiter = (wDeltaIter->second).begin(); deltaiter != (wDeltaIter->second).end(); deltaiter++) {
				deltaiter2 = (wDeltaIter2->second).find( deltaiter->first );

				DeltaReaction = deltaiter->second/fNbOfScoredEvents;
				if ( fNbOfScoredEvents > 1 ) {
					DeltaError = sqrt( (1.0/(fNbOfScoredEvents-1)) * ( (deltaiter2->second)/fNbOfScoredEvents - DeltaReaction*DeltaReaction));
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
}


void TsScoreWithIRT::Clear() {
	fGValuePerSpeciePerTime.clear();
	fGValuePerSpeciePerTime2.clear();
	fMoleculesPerSpeciePerTime.clear();
	fMoleculesPerSpeciePerTime2.clear();
	fIRTExecutionTime = 0;
	fNbOfScoredEvents = 0;
	if (fReportDelta) {
		fDeltaGPerReactionPerTime.clear();
		fDeltaGPerReactionPerTime2.clear();
	}
	UpdateFileNameForUpcomingRun();
}


