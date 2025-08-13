//
// ********************************************************************
// *                                                                  *
// * This file is part of the TOPAS-nBio extensions to the            *
// *   TOPAS Simulation Toolkit.                                      *
// * The TOPAS-nBio extensions are freely available under the license *
// *   agreement set forth at: https://topas-nbio.readthedocs.io/     *
// *                                                                  *
// ********************************************************************

#ifndef TsSBSScoreGValue_hh
#define TsSBSScoreGValue_hh

#include "Randomize.hh"
#include "TsVNtupleScorer.hh"

#include <stdint.h>

class G4MolecularConfiguration;

class TsSBSScoreGValue : public TsVNtupleScorer
{
public:
	TsSBSScoreGValue(TsParameterManager* pM, TsMaterialManager* mM, TsGeometryManager* gM, TsScoringManager* scM, TsExtensionManager* eM,
					 G4String scorerName, G4String quantity, G4String outFileName, G4bool isSubScorer);
	~TsSBSScoreGValue() = default;

	G4bool ProcessHits(G4Step*, G4TouchableHistory*) override;
	void AbsorbResultsFromWorkerScorer(TsVScorer* workerScorer) override;

protected:
	void AccumulateEvent() override;
	void Output() override;
	void Clear() override;

	// Output variables
	G4double fGValue;
	G4double fGValueError;
	G4double fTime;
	G4String fMoleculeName;

	std::map<G4String, std::map<G4double, G4double>> fGValuePerSpeciePerTime;
	std::map<G4String, std::map<G4double, G4double>> fGValuePerSpeciePerTime2;

private:
	TsParameterManager* fPm;

	G4String fMoleculeCounterName;

	G4double fEnergyDepositPerEvent;
	//    G4double* fTimeToRecord;
	std::vector<G4double> fTimesToRecord;

	G4int fNbTimeToRecord;
	G4int fNbOfScoredEvents;
	G4double fEnergyLossKill;
	G4double fEnergyLossAbort;
	G4double fEnergyLoss;
	G4double fMaximumTrackLength;
	G4double fTotalTrackLength;
	G4int fNbOfMoleculesToScavenge;
	G4int* fMoleculeIDToScavenge;
	G4double* fScavengingCapacity;

	std::map<G4double, G4double> fScavengerProducts;
};

#endif
