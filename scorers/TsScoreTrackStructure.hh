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

#ifndef TsScoreTrackStructure_hh
#define TsScoreTrackStructure_hh

#include "TsVNtupleScorer.hh"

class TsScoreTrackStructure : public TsVNtupleScorer
{
public:
	TsScoreTrackStructure(TsParameterManager* pM, TsMaterialManager* mM, TsGeometryManager* gM, TsScoringManager* scM, TsExtensionManager* eM,
							 G4String scorerName, G4String quantity, G4String outFileName, G4bool isSubScorer);

	virtual ~TsScoreTrackStructure();

	G4bool ProcessHits(G4Step*,G4TouchableHistory*);

private:
	G4float fPosX;
	G4float fPosY;
	G4float fPosZ;
	G4float fEdep;
	G4String fParticleName;
	G4int fEvt;
	G4int fTrackID;
	G4int fParentID;
	G4String fVolumeName;

	G4bool fIncludePrimaryPositions;

	G4bool fIncludeEventID;
	G4bool fIncludeTrackID;
	G4bool fIncludeParentID;
	G4bool fIncludeParticleName;
	G4bool fIncludeVolumeName;
};

#endif
