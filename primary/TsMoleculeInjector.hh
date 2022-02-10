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
// Authors: Alejandro Bertolet, Giulia Tamborino, Jan Schuemann

#ifndef TsMoleculeInjector_hh
#define TsMoleculeInjector_hh

#include "TsVGenerator.hh"
#include "G4MoleculeGun.hh"

class TsMoleculeInjector : public TsVGenerator
{
public:
	TsMoleculeInjector(TsParameterManager* pM, TsGeometryManager* gM, TsGeneratorManager* pgM, G4String sourceName);
	virtual ~TsMoleculeInjector();

	void GeneratePrimaries(G4Event* anEvent);

protected:
	void ResolveParameters();

private:
	//std::unique_ptr<G4MoleculeGun> fMoleculeGun;
	G4MoleculeGun* fMoleculeGun;
	G4bool MoleculeExists(G4String name);

	std::map<G4String, G4String> fExistingMolecules;

	G4String fSolidName;
	G4String fMolecule;
	G4bool fInjectOnSurface;
	G4int fNumberOfMolecules;
	G4double fPartialPressuremmHg;
	G4double fTemperature;
};

#endif
