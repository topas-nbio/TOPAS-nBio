//
// ********************************************************************
// *																  *
// * This file is part of the TOPAS-nBio extensions to the			  *
// *   TOPAS Simulation Toolkit.									  *
// * The TOPAS-nBio extensions are freely available under the license *
// *   agreement set forth at: https://topas-nbio.readthedocs.io/	  *
// *																  *
// ********************************************************************
//
// Authors: Alejandro Bertolet, Jan Schuemann

#ifndef TsFociAnalysis_hh
#define TsFociAnalysis_hh

#include "TsVNtupleScorer.hh"

class TsFociAnalysis
{
public:
	TsFociAnalysis();
	~TsFociAnalysis();

	G4int GetNumberOfFoci(std::vector<G4ThreeVector> dsbPositions);
	G4double GetDistance(G4ThreeVector a, G4ThreeVector b);
	G4bool CheckIfAnyIndexIsAvailable(std::vector<G4bool> indexIsAvailable);

	void inline SetFociSize(G4double v)								{ fFociSize = v; }
	void inline Set3DFociImage(G4bool v)							{ fGet3DFociImage = v; }
	void inline Set2DFociImages(G4bool v)							{ fGet2DFociImage = v; }
	void inline SetPlanesFor2DFociImages(std::vector<G4String> v)	{ f2DPlanesForFociImage = v; }
	void inline SetPSFShape(G4String v)								{ fMicroscopePSFShape = v; }
	void inline SetPSFWidth(G4double v)								{ fMicroscopePSFWidth = v; }

private:
	G4double fFociSize;
	G4bool fGet3DFociImage;
	G4bool fGet2DFociImage;
	std::vector<G4String> f2DPlanesForFociImage;
	G4String fMicroscopePSFShape;
	G4double fMicroscopePSFWidth;
};

#endif
