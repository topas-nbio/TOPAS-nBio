// Particle Generator for MoleculeInjector
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

#include "TsMoleculeInjector.hh"

#include "TsParameterManager.hh"

#include "G4VSolid.hh"
#include "G4LogicalVolume.hh"
#include "G4LogicalVolumeStore.hh"
#include "G4VPhysicalVolume.hh"
#include "G4VoxelLimits.hh"
#include "G4AffineTransform.hh"

#include "G4SystemOfUnits.hh"
#include "G4PhysicalConstants.hh"

#include "G4DNAChemistryManager.hh"
#include "G4DNAMolecularReactionTable.hh"
#include "G4MoleculeTable.hh"
#include "G4Molecule.hh"
#include "G4MoleculeGun.hh"

#include "G4Electron_aq.hh"
#include "G4Hydrogen.hh"
#include "G4H2.hh"
#include "G4H3O.hh"
#include "G4H2O2.hh"
#include "TsO2.hh"
#include "G4OH.hh"

#include "Randomize.hh"

TsMoleculeInjector::TsMoleculeInjector(TsParameterManager* pM, TsGeometryManager* gM, TsGeneratorManager* pgM, G4String sourceName)
					  : TsVGenerator(pM, gM, pgM, sourceName)
{
	fExistingMolecules["solvatedelectron"] = "e_aq";
    fExistingMolecules["hydrogen"] = "H";
    fExistingMolecules["hydronium"] = "H3Op";
    fExistingMolecules["dyhydrogen"] = "H2";
    fExistingMolecules["hydroxide"] = "OHm";
    fExistingMolecules["hydrogenperoxide"] = "H2O2";
    fExistingMolecules["oxygen"] = "O2";
    fExistingMolecules["hydroxyl"] = "OH";
    ResolveParameters();
}

TsMoleculeInjector::~TsMoleculeInjector() {}

void TsMoleculeInjector::ResolveParameters()
{
	TsVGenerator::ResolveParameters();

	fSolidName = fPm->GetStringParameter(GetFullParmName("Component"));
	if (!fPm->ParameterExists(GetFullParmName("Molecule")))
	{
		G4cerr << "TOPAS is exiting due to an error in the MoleculeInjector." << G4endl;
		G4cerr << "Molecule has to be specified." << G4endl;
		exit(1);
	}
	fMolecule = fPm->GetStringParameter(GetFullParmName("Molecule"));
	fMolecule.toLower();
	if (!MoleculeExists(fMolecule))
	{
		G4cerr << "TOPAS is exiting due to an error in MoleculeInjector." << G4endl;
		G4cerr << "Molecule " << fMolecule << " was not found in the database." << G4endl;
		exit(1);
	}

	fInjectOnSurface = false;
	if (fPm->ParameterExists(GetFullParmName("InsufflateOnSurface")))
		fInjectOnSurface = fPm->GetBooleanParameter(GetFullParmName("InsufflateOnSurface"));

	fPartialPressuremmHg = 0;
	if (fPm->ParameterExists(GetFullParmName("PartialPressure_mmHg")))
		fPartialPressuremmHg = fPm->GetUnitlessParameter(GetFullParmName("PartialPressure_mmHg"));

	fTemperature = 20.0;
	if (fPm->ParameterExists(GetFullParmName("Temperature_Celsius")))
		fTemperature = fPm->GetUnitlessParameter(GetFullParmName("Temperature_Celsius"));

	fNumberOfMolecules = 0;
	if (fPm->ParameterExists(GetFullParmName("NumberOfMolecules")))
		fNumberOfMolecules = fPm->GetIntegerParameter(GetFullParmName("NumberOfMolecules"));


	if (fPm->ParameterExists(GetFullParmName("PartialPressure_mmHg")) && fPm->ParameterExists(GetFullParmName("NumberOfMolecules")))
	{
		G4cerr << "TOPAS is exiting due to an error in the molecule Injector." << G4endl;
		G4cerr << "PartialPressure and NumberOfMolecules cannot be specified at the same time." << G4endl;
		exit(1);
	}
	if (!fPm->ParameterExists(GetFullParmName("PartialPressure_mmHg")) && !fPm->ParameterExists(GetFullParmName("NumberOfMolecules")))
		{
			G4cerr << "TOPAS is exiting due to an error in the molecule Injector." << G4endl;
			G4cerr << "Either PartialPressure or NumberOfMolecules have to be specified." << G4endl;
			exit(1);
		}
}

void TsMoleculeInjector::GeneratePrimaries(G4Event* anEvent)
{
	fMoleculeGun = (G4MoleculeGun*)new G4MoleculeGun();
	G4DNAChemistryManager::Instance()->SetGun(fMoleculeGun);
	G4VSolid* solid = G4LogicalVolumeStore::GetInstance()->GetVolume(fSolidName)->GetSolid();

	fExistingMolecules["solvatedelectron"] = "e_aq";
    fExistingMolecules["hydrogen"] = "H";
    fExistingMolecules["hydronium"] = "H3Op";
    fExistingMolecules["dyhydrogen"] = "H2";
    fExistingMolecules["hydroxide"] = "OHm";
    fExistingMolecules["hydrogenperoxide"] = "H2O2";

    if (strstr(fMolecule, "solvatedelectron"))
    	G4MoleculeTable::Instance()->CreateConfiguration("e_aq", G4Electron_aq::Definition());
    else if (strstr(fMolecule, "hydrogen"))
    	G4MoleculeTable::Instance()->CreateConfiguration("H", G4Hydrogen::Definition());
    else if (strstr(fMolecule, "hydronium"))
    	G4MoleculeTable::Instance()->CreateConfiguration("H3Op", G4H3O::Definition());
    else if (strstr(fMolecule, "dyhydrogen"))
    	G4MoleculeTable::Instance()->CreateConfiguration("H2", G4H2::Definition());
    else if (strstr(fMolecule, "hydroxide"))
    {
    	G4MolecularConfiguration* OHm =	G4MoleculeTable::Instance()->CreateConfiguration("OHm", G4OH::Definition(), -1, 5.0e-9 * (m2 / s));
    	OHm->SetMass(17.0079 * g / Avogadro * c_squared);
    }
    else if (strstr(fMolecule, "hydrogenperoxide"))
   		G4MoleculeTable::Instance()->CreateConfiguration("H2O2", G4H2O2::Definition());
    else if (strstr(fMolecule, "oxygen"))
		G4MoleculeTable::Instance()->CreateConfiguration("O2", TsO2::Definition());
	else if (strstr(fMolecule, "hydroxyl"))
		G4MoleculeTable::Instance()->CreateConfiguration("OH", G4OH::Definition());

	if (fPartialPressuremmHg > 0)
	{
		G4double vol = solid->GetCubicVolume() / 1e9; // from mm3 to m3
		G4double R = 8.31446261815324;
		G4double tempKelvin = fTemperature + STP_Temperature;
		G4double pPa = fPartialPressuremmHg * 133.322;
		G4double nMol = (vol * pPa) / (R * tempKelvin);
		fNumberOfMolecules = (G4int)floor(nMol * Avogadro);
		G4cout << "vol: " << vol<< " - temp: " << tempKelvin << " - nMol: " << nMol << " - Number of molecules: " << fNumberOfMolecules << G4endl;
	}
	for (G4int i = 0; i < fNumberOfMolecules; i++)
	{
		G4ThreeVector position;
		if (fInjectOnSurface)
			position = solid->GetPointOnSurface();
		else
		{
			// Getting limits of the solid for molecules to be insufflated
			G4VoxelLimits voxelLimits;
			G4AffineTransform affineTrans;
			G4double xmin, xmax, ymin, ymax, zmin, zmax;

			solid->CalculateExtent(kXAxis, voxelLimits, affineTrans, xmin, xmax);
			solid->CalculateExtent(kYAxis, voxelLimits, affineTrans, ymin, ymax);
			solid->CalculateExtent(kZAxis, voxelLimits, affineTrans, zmin, zmax);

			G4double x, y, z;
			while (1)
			{
				x = G4RandFlat::shoot(xmin, xmax);
				y = G4RandFlat::shoot(ymin, ymax);
				z = G4RandFlat::shoot(zmin, zmax);
				if (solid->Inside(G4ThreeVector(x, y, z)) == kInside) break;
			}
			position = G4ThreeVector(x, y, z);
		}

		G4double iniTime = 5*picosecond;
		fMoleculeGun->AddMolecule(fExistingMolecules[fMolecule], position, iniTime);
	}
}

G4bool TsMoleculeInjector::MoleculeExists(G4String name)
{
	if (fExistingMolecules.find(name) == fExistingMolecules.end()) return false;
	else return true;
}
