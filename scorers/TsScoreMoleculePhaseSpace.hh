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

#ifndef TsScoreMoleculePhaseSpace_hh
#define TsScoreMoleculePhaseSpace_hh

#include "TsVNtupleScorer.hh"

class TsScoreMoleculePhaseSpace : public TsVNtupleScorer
{
public:
	TsScoreMoleculePhaseSpace(TsParameterManager* pM, TsMaterialManager* mM, TsGeometryManager* gM, TsScoringManager* scM, TsExtensionManager* eM,
						 G4String scorerName, G4String quantity, G4String outFileName, G4bool isSubScorer);
	
	virtual ~TsScoreMoleculePhaseSpace();
	
	G4bool ProcessHits(G4Step*,G4TouchableHistory*);
	
private:
	G4int    fEvt;
	G4String fParticleName;
	G4float  fPosX;
	G4float  fPosY;
	G4float  fPosZ;
	G4float  fTime;

	G4double fTimeCut;

};
#endif

