// Scorer for ClusterSizeInRandomCylinders
//
// ********************************************************************
// *                                                                  *
// * This file is part of the TOPAS-nBio extensions to the            *
// *   TOPAS Simulation Toolkit.                                      *
// * The TOPAS-nBio extensions are freely available under the license *
// *   agreement set forth at: https://topas-nbio.readthedocs.io/     *
// *                                                                  *
// ********************************************************************


#include "TsScoreClusterSizeInRandomCylinders.hh"

#include "G4SystemOfUnits.hh"
#include "G4Track.hh"
#include "G4VProcess.hh"
#include "G4ProcessVector.hh"
#include "G4ProcessTable.hh"
#include "G4Electron.hh"
#include "G4ParticleDefinition.hh"
#include <fstream>

#include <map>

TsScoreClusterSizeInRandomCylinders::TsScoreClusterSizeInRandomCylinders(TsParameterManager* pM, TsMaterialManager* mM, TsGeometryManager* gM, TsScoringManager* scM, TsExtensionManager* eM,
                                                                         G4String scorerName, G4String quantity, G4String outFileName, G4bool isSubScorer)
:TsVNtupleScorer(pM, mM, gM, scM, eM, scorerName, quantity, outFileName, isSubScorer), fPm(pM)
{
    SetUnit("");
    
    fNbOfAlgo = 1;
				
    fNtuple->RegisterColumnI(&fClusterSize, "ClusterSize");
    fNtuple->RegisterColumnI(&fClusterFreq, "IonisationClusterFreq");
    fNtuple->RegisterColumnI(&fEClusterFreq, "ExcitationClusterFreq");
    fNtuple->RegisterColumnI(&fTClusterFreq, "TotalClusterFreq");
      
    fElectron = G4Electron::ElectronDefinition();
    
    fRepeatedGeometry = false;
    if (fPm->ParameterExists(GetFullParmName("IsRepeatedGeometry"))){
        fRepeatedGeometry = fPm->GetBooleanParameter(GetFullParmName("IsRepeatedGeometry"));
    }

}


TsScoreClusterSizeInRandomCylinders::~TsScoreClusterSizeInRandomCylinders() {
}


G4bool TsScoreClusterSizeInRandomCylinders::ProcessHits(G4Step* aStep, G4TouchableHistory*)
{
    if (!fIsActive) {
        fSkippedWhileInactive++;
        return false;
    }
    
    const G4String processName = aStep->GetPostStepPoint()->GetProcessDefinedStep()->GetProcessName();
    G4bool scoreIon = false;
    G4bool scoreExc = false;
           
#if GEANT4_VERSION_MAJOR >= 11
    if ( G4StrUtil::contains(processName,"Ionisation") )
#else
    if ( processName.contains("Ionisation") )
#endif
	 scoreIon = true;

#if GEANT4_VERSION_MAJOR >= 11
    if ( G4StrUtil::contains(processName,"Excitation"))
#else
    if ( processName.contains("Excitation"))
#endif
	 scoreExc = true;

    if ( scoreIon || scoreExc ) { 
	G4int index = aStep->GetPreStepPoint()->GetTouchableHandle()->GetCopyNumber();  
        if ( fRepeatedGeometry && index == 0 )
            return false;

	G4int indexParent = fRepeatedGeometry ? aStep->GetPreStepPoint()->GetTouchableHandle()->GetCopyNumber(1):0;  

	if ( scoreIon ) 
	   fIonisations[indexParent][index]++; 
	if ( scoreExc ) 
	   fExcitations[indexParent][index]++;

	fTotal[indexParent][index]++;
    }
    return false;
}


void TsScoreClusterSizeInRandomCylinders::UserHookForEndOfEvent() {
    for ( int i = 0; i < fNbOfAlgo; i++ ) {
        std::map<G4int, G4int>::iterator ionIter;
        for ( ionIter = fIonisations[i].begin(); ionIter != fIonisations[i].end(); ionIter++ ) {
            G4int clusterSize = ionIter->second; 
            if ( fHistogram.find(clusterSize) == fHistogram.end() ) {
	        fHistogram[clusterSize] = 1;
            } else {
                fHistogram[clusterSize]++;
	    }
        }
        fIonisations[i].erase(fIonisations[i].begin(), fIonisations[i].end());

        for ( ionIter = fExcitations[i].begin(); ionIter != fExcitations[i].end(); ionIter++ ) {
            G4int clusterSize = ionIter->second; 
            if ( fHistogramE.find(clusterSize) == fHistogramE.end() ) {
	        fHistogramE[clusterSize] = 1;
            } else {
                fHistogramE[clusterSize]++;
	    }
        }
        fExcitations[i].erase(fExcitations[i].begin(), fExcitations[i].end());
 
        for ( ionIter = fTotal[i].begin(); ionIter != fTotal[i].end(); ionIter++ ) {
            G4int clusterSize = ionIter->second; 
            if ( fHistogramT.find(clusterSize) == fHistogramT.end() ) {
	        fHistogramT[clusterSize] = 1;
            } else {
                fHistogramT[clusterSize]++;
	    }
        }
        fTotal[i].erase(fTotal[i].begin(), fTotal[i].end());
 
    }
    fIonisations.clear();
    fExcitations.clear();
    fTotal.clear();
    fTotalPath = 0.0;    
}


void TsScoreClusterSizeInRandomCylinders::AbsorbResultsFromWorkerScorer(TsVScorer* workerScorer) {
        TsVNtupleScorer::AbsorbResultsFromWorkerScorer(workerScorer);

        TsScoreClusterSizeInRandomCylinders* workerClusterScorer = dynamic_cast<TsScoreClusterSizeInRandomCylinders*>(workerScorer);

        std::map<G4int, G4int>::iterator wIter;
        std::map<G4int, G4int>::iterator mIter;

        for ( wIter = workerClusterScorer->fHistogram.begin(); wIter != workerClusterScorer->fHistogram.end(); wIter++) {
                mIter = fHistogram.find( wIter->first );
                if ( mIter == fHistogram.end() ) {
                        fHistogram.insert(std::pair<G4int, G4int>(wIter->first, wIter->second));
                } else {
			mIter->second += wIter->second;
                }
        }

        workerClusterScorer->fHistogram.clear();

        for ( wIter = workerClusterScorer->fHistogramE.begin(); wIter != workerClusterScorer->fHistogramE.end(); wIter++) {
                mIter = fHistogramE.find( wIter->first );
                if ( mIter == fHistogramE.end() ) {
                        fHistogramE.insert(std::pair<G4int, G4int>(wIter->first, wIter->second));
                } else {
			mIter->second += wIter->second;
                }
        }
        workerClusterScorer->fHistogramE.clear();

        for ( wIter = workerClusterScorer->fHistogramT.begin(); wIter != workerClusterScorer->fHistogramT.end(); wIter++) {
                mIter = fHistogramT.find( wIter->first );
                if ( mIter == fHistogramT.end() ) {
                        fHistogramT.insert(std::pair<G4int, G4int>(wIter->first, wIter->second));
                } else {
			mIter->second += wIter->second;
                }
        }
        workerClusterScorer->fHistogramT.clear();

}


void TsScoreClusterSizeInRandomCylinders::Output() {

        G4int iEnd = fHistogramT.rbegin()->first;

        std::map<G4int,G4int>::iterator itr;

        for ( int i = 0; i <= iEnd; i++ ) {
	    if ( fHistogram.find(i) == fHistogram.end() ) {
                fClusterFreq = 0;
            } else {
                fClusterFreq = fHistogram[i];
            }
            fClusterSize = i;

	    if ( fHistogramE.find(i) == fHistogramE.end() ) {
                fEClusterFreq = 0;
            } else {
                fEClusterFreq = fHistogramE[i];
            }

	    if ( fHistogramT.find(i) == fHistogramT.end() ) {
                fTClusterFreq = 0;
            } else {
                fTClusterFreq = fHistogramT[i];
            }

            fNtuple->Fill();
        }

        fNtuple->Write();
        UpdateFileNameForUpcomingRun();
}

void TsScoreClusterSizeInRandomCylinders::Clear() {
        fHistogram.clear();
        fIonisations.clear();
        fHistogramE.clear();
        fExcitations.clear();
        fHistogramT.clear();
        fTotal.clear();
}

