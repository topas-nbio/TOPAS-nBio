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

#ifndef TsScoreNumberOfMoleculesAtTime_hh
#define TsScoreNumberOfMoleculesAtTime_hh

#include "TsVNtupleScorer.hh"

#include <stdint.h>

class G4MolecularConfiguration;

class TsScoreNumberOfMoleculesAtTime : public TsVNtupleScorer
{
public:
	TsScoreNumberOfMoleculesAtTime(TsParameterManager* pM, TsMaterialManager* mM, TsGeometryManager* gM, TsScoringManager* scM, TsExtensionManager* eM,
								   G4String scorerName, G4String quantity, G4String outFileName, G4bool isSubScorer);
	~TsScoreNumberOfMoleculesAtTime() = default;

	G4bool ProcessHits(G4Step*, G4TouchableHistory*) override;
	void AbsorbResultsFromWorkerScorer(TsVScorer* workerScorer) override;

protected:
	void AccumulateEvent() override;
	void Output() override;
	void Clear() override;

	// Output variables
	G4float fMeanNumberOfMolecules;
	G4int fNumberOfMolecules;
	G4float fTime;
	G4String fMoleculeName;

private:
	TsParameterManager* fPm;

	G4String fMoleculeCounterName;

	std::vector<G4double> fTimesToRecord;
	G4int fNbOfScoredEvents;
	G4double fEnergyLoss;
	G4double fEnergyLossKill;
	G4double fEnergyLossAbort;
	G4double fMaximumTrackLength;
	G4double fTotalTrackLength;

	std::map<G4String, std::map<G4double, G4int>> fMoleculeCountPerIndexPerTime;
};

#endif
