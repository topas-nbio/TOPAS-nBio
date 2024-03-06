//
// ********************************************************************
// *                                                                  *
// * This file is part of the TOPAS-nBio extensions to the            *
// *   TOPAS Simulation Toolkit.                                      *
// * The TOPAS-nBio extensions are freely available under the license *
// *   agreement set forth at: https://topas-nbio.readthedocs.io/     *
// *                                                                  *
// ********************************************************************

#ifndef TsScoreClusterSizeInRandomCylinders_hh
#define TsScoreClusterSizeInRandomCylinders_hh

#include "TsVNtupleScorer.hh"

#include <map>
#include <fstream>

class TsScoringManager;
class TsParameterManager;
class TsScoreClusterSizeInRandomCylinders : public TsVNtupleScorer
{
public:
    TsScoreClusterSizeInRandomCylinders(TsParameterManager* pM, TsMaterialManager* mM, TsGeometryManager* gM, TsScoringManager* scM, TsExtensionManager* eM,
				G4String scorerName, G4String quantity, G4String outFileName, G4bool isSubScorer);
	
    virtual ~TsScoreClusterSizeInRandomCylinders();

	G4bool ProcessHits(G4Step*,G4TouchableHistory*);
    
	void UserHookForEndOfEvent();
        void AbsorbResultsFromWorkerScorer(TsVScorer* workerScorer);
protected:
	void Output();
	void Clear();

private:
    TsParameterManager* fPm;
    G4double fTotalPath;
    G4bool fRepeatedGeometry;

    G4int fNbOfAlgo;
    G4int fClusterSize;
    G4int fClusterFreq;
    G4int fEClusterFreq;
    G4int fTClusterFreq;

    std::map<G4int, std::map<G4int, G4int>> fIonisations;
    std::map<G4int, std::map<G4int, G4int>> fExcitations;
    std::map<G4int, std::map<G4int, G4int>> fTotal;
    std::map<G4int, G4int> fHistogram; 
    std::map<G4int, G4int> fHistogramE; 
    std::map<G4int, G4int> fHistogramT; 

    G4ParticleDefinition* fElectron;
    std::ofstream fOutputAscii;
   
};
#endif
