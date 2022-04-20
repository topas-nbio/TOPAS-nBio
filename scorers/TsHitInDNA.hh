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

#ifndef TsHitInDNA_hh
#define TsHitInDNA_hh

#include "G4UnitsTable.hh"

class TsHitInDNA
{
public:
	TsHitInDNA();
	~TsHitInDNA();

	// Get methods
	inline G4int GetEventID()					{ return fEventID; }
	inline G4double GetEdep()					{ return fEdep; }
	inline G4String GetParticleName()			{ return fParticleName; }
	inline G4int GetBasePairID()				{ return fBasePairID; }
	inline G4int GetStrandNumber()				{ return fStrandID; }
	inline G4int GetDNAComponentID()			{ return fComponentID; }
	inline G4ThreeVector GetPosition()			{ return fPos; }

	// Set methods
	inline void SetEventID(G4int id)					{ fEventID = id; }
	inline void SetEdep(G4double edep)					{ fEdep = edep; }
	inline void SetParticleName(G4String name)			{ fParticleName = name; }
	inline void SetBasePairID(G4int bpid)				{ fBasePairID = bpid; }
	inline void SetStrandNumber(G4int stnum)			{ fStrandID = stnum; }
	inline void SetDNAComponentID(G4int cid)			{ fComponentID = cid; }
	inline void SetPosition(G4ThreeVector pos)			{ fPos = pos; }

private:
	G4int fEventID;

	G4double fEdep;
	G4String fParticleName;

	G4int fBasePairID;
	G4int fStrandID;
	G4int fComponentID;
	G4ThreeVector fPos;

};

#endif
