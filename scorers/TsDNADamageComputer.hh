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

#ifndef TsDNADamageComputer_hh
#define TsDNADamageComputer_hh

#include "TsVNtupleScorer.hh"
#include "TsHitInDNA.hh"

class TsDNADamageComputer
{
public:
	TsDNADamageComputer();
	~TsDNADamageComputer();

	// Set methods
	inline void SetScoreDirectDamage(G4bool v)					{ fScoreDirectDamage = v; }
	inline void SetScoreQuasiDirectDamage(G4bool v)				{ fScoreQuasiDirectDamage = v; }
	inline void SetScoreIndirectDamage(G4bool v)				{ fScoreIndirectDamage = v; }
	inline void SetScoreOnBases(G4bool v)						{ fScoreOnBases = v; }
	inline void SetScoreOnBackbones(G4bool v)					{ fScoreOnBackbones = v; }
	inline void SetExcludeShortFragments(G4bool v)				{ fExcludeShortFragments = v; }
	inline void SetLowerThresholdForFragmentDetection(G4int v)	{ fLowerThresholdForFragmentDetection = v; }
	inline void SetUpperThresholdForFragmentDetection(G4int v)	{ fUpperThresholdForFragmentDetection = v; }

	// Get methods
	inline G4int GetSB()										{ return fNumSB; }
	inline G4int GetSBDirect()									{ return fNumSBDirect; }
	inline G4int GetSBQuasiDirect()								{ return fNumSBQuasiDirect; }
	inline G4int GetSBIndirect()								{ return fNumSBIndirect; }
	inline G4int GetSSB()										{ return fNumSSB; }
	inline G4int GetSSBDirect()									{ return fNumSSBDirect; }
	inline G4int GetSSBQuasiDirect()							{ return fNumSSBQuasiDirect; }
	inline G4int GetSSBIndirect()								{ return fNumSSBIndirect; }
	inline G4int GetDSB()										{ return fNumDSB; }
	inline G4int GetDSBDirect()									{ return fNumDSBDirect; }
	inline G4int GetDSBIndirect()								{ return fNumDSBIndirect; }
	inline G4int GetDSBDirectIndirect()							{ return fNumDSBDirectIndirect; }
	inline G4int GetDSBDirectQuasiDirect()						{ return fNumDSBDirectQuasiDirect; }
	inline G4int GetDSBQuasiDirectQuasiDirect()					{ return fNumDSBQuasiDirectQuasiDirect; }
	inline G4int GetDSBIndirectQuasiDirect()					{ return fNumDSBIndirectQuasiDirect; }
	inline G4int GetBD()										{ return fNumBaseDamage; }
	inline G4int GetBDDirect()									{ return fNumBaseDamageDirect; }
	inline G4int GetBDQuasiDirect()								{ return fNumBaseDamageQuasiDirect; }
	inline G4int GetBDIndirect()								{ return fNumBaseDamageIndirect; }

private:
	// Options to define damage


	// Options for whether counting types of damage
	G4bool fScoreDirectDamage;
	G4bool fScoreQuasiDirectDamage;
	G4bool fScoreIndirectDamage;

	G4bool fScoreOnBases;
	G4bool fScoreOnBackbones;

	G4bool fExcludeShortFragments;
	G4int fLowerThresholdForFragmentDetection;
	G4int fUpperThresholdForFragmentDetection;

	// Quantification of damage
	G4int fNumSB;
	G4int fNumSBDirect;
	G4int fNumSBQuasiDirect;
	G4int fNumSBIndirect;
	G4int fNumSSB;
	G4int fNumSSBDirect;
	G4int fNumSSBQuasiDirect;
	G4int fNumSSBIndirect;
	G4int fNumDSB;
	G4int fNumDSBDirect;
	G4int fNumDSBIndirect;
	G4int fNumDSBDirectIndirect;
	G4int fNumDSBDirectQuasiDirect;
	G4int fNumDSBQuasiDirectQuasiDirect;
	G4int fNumDSBIndirectQuasiDirect;
	G4int fNumBaseDamage;
	G4int fNumBaseDamageDirect;
	G4int fNumBaseDamageQuasiDirect;
	G4int fNumBaseDamageIndirect;


};

#endif
