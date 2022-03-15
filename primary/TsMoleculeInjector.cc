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

#include "G4PhysicalConstants.hh"

#include "G4DNAChemistryManager.hh"
#include "G4DNAMolecularReactionTable.hh"
#include "G4MoleculeTable.hh"
#include "G4Molecule.hh"
#include "G4MoleculeGun.hh"
#include "G4VMoleculeCounter.hh"

#include "G4Electron_aq.hh"
#include "G4Hydrogen.hh"
#include "G4H2.hh"
#include "G4H3O.hh"
#include "G4H2O2.hh"
#include "G4H2O.hh"
#include "TsO2.hh"


#include "Randomize.hh"

TsMoleculeInjector::TsMoleculeInjector(TsParameterManager* pM, TsGeometryManager* gM, TsGeneratorManager* pgM, G4String sourceName)
					  : TsVGenerator(pM, gM, pgM, sourceName), fDifferentTimesForEvents(false)
{
	fExistingMolecules["solvatedelectron"] = "e_aq";
	fExistingMolecules["e_aq^-1"] = "e_aq";
	fExistingMolecules["hydrogen"] = "H";
	fExistingMolecules["H^0"] = "H";
	fExistingMolecules["hydronium"] = "H3Op";
	fExistingMolecules["H3O^1"] = "H3Op";
	fExistingMolecules["dyhydrogen"] = "H2";
	fExistingMolecules["H_2^0"] = "H2";
	fExistingMolecules["hydroxide"] = "OHm";
	fExistingMolecules["OH^0"] = "OHm";
	fExistingMolecules["hydrogenperoxide"] = "H2O2";
	fExistingMolecules["H2O2^0"] = "H2O2";
	fExistingMolecules["water"] = "H2O";
	fExistingMolecules["H2O^0"] = "H2O";
	fExistingMolecules["oxygen"] = "O2";
	fExistingMolecules["O2^0"] = "O2";

    ResolveParameters();
}

TsMoleculeInjector::~TsMoleculeInjector() {}

void TsMoleculeInjector::ResolveParameters()
{
	TsVGenerator::ResolveParameters();

	fSolidName = fPm->GetStringParameter(GetFullParmName("Component"));
	fUseMoleculePhaseSpace = false;
	if (fPm->ParameterExists(GetFullParmName("MoleculePhaseSpaceFileName")))
	{
		fUseMoleculePhaseSpace = true;
		fFileName = fPm->GetStringParameter(GetFullParmName("MoleculePhaseSpaceFileName"));
		if (fPm->ParameterExists(GetFullParmName("TimesForEachEvent")))
		{
			fDifferentTimesForEvents = true;
			fNEvents = fPm->GetVectorLength(GetFullParmName("TimesForEachEvent"));
			G4double* times = fPm->GetDoubleVector(GetFullParmName("TimesForEachEvent"), "Time");
			for (G4int i=0; i < fNEvents; i++)
				fTimesForEvents.push_back(times[i]);
		}
	}
	else
	{
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

}

void TsMoleculeInjector::GeneratePrimaries(G4Event* anEvent)
{

	fMoleculeGun = (G4MoleculeGun*)new G4MoleculeGun();
	G4DNAChemistryManager::Instance()->SetGun(fMoleculeGun);
	G4VSolid* solid = G4LogicalVolumeStore::GetInstance()->GetVolume(fSolidName)->GetSolid();

    G4MoleculeTable::Instance()->CreateConfiguration("e_aq", G4Electron_aq::Definition());
    G4MoleculeTable::Instance()->CreateConfiguration("H", G4Hydrogen::Definition());
    G4MoleculeTable::Instance()->CreateConfiguration("H3Op", G4H3O::Definition());
    G4MoleculeTable::Instance()->CreateConfiguration("H2", G4H2::Definition());
    //G4MoleculeTable::Instance()->CreateConfiguration("OHm", G4OH::Definition(), -1, 5.0e-9 * (m2 / s));;
    //fOHm->SetMass(17.0079 * g / Avogadro * c_squared);
    G4MoleculeTable::Instance()->CreateConfiguration("H2O2", G4H2O2::Definition());
    G4MoleculeTable::Instance()->CreateConfiguration("OH", G4OH::Definition());
    G4MoleculeTable::Instance()->CreateConfiguration("H2O", G4H2O::Definition());
    G4MoleculeTable::Instance()->CreateConfiguration("O2", TsO2::Definition());

	if (fUseMoleculePhaseSpace)
		ShootMoleculePhaseSpace();
	else
	{
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

			G4double iniTime = 1*picosecond;
			fMoleculeGun->AddMolecule(fExistingMolecules[fMolecule], position, iniTime);
		}
	}
}

G4bool TsMoleculeInjector::ShootMoleculePhaseSpace()
{
	G4String dataFileSpec = fFileName + ".phsp";
	fDataFile.open(dataFileSpec);
	if (!fDataFile)
	{
		G4cerr << "Error opening phase space data file: " << dataFileSpec << G4endl;
		fPm->AbortSession(1);
	}
	G4int firstEventID = -1;
	while (fDataFile.tellg() != -1)
	{
		G4String asciiLine = "";
		getline(fDataFile, asciiLine);
		if (!fDataFile.good()) return true;
		std::istringstream input(asciiLine);
		input >> fMoleculeID >> fPosX >> fPosY >> fPosZ >> fGlobalTime >> fEvtID;
		if (firstEventID < 0) firstEventID = fEvtID;
		G4double initialTime = 0;
		if (fDifferentTimesForEvents)
		{
			if (fEvtID-firstEventID < fNEvents)
				initialTime = fTimesForEvents[fEvtID-firstEventID];
		}
		G4ThreeVector position;
		position = G4ThreeVector(fPosX*um, fPosY*um, fPosZ*um);
		if (fExistingMolecules[fMoleculeID] != "H2O")
			fMoleculeGun->AddMolecule(fExistingMolecules[fMoleculeID], position, initialTime+fGlobalTime*picosecond);
	}

	fDataFile.close();
	return true;
}


G4bool TsMoleculeInjector::MoleculeExists(G4String name)
{
	if (fExistingMolecules.find(name) == fExistingMolecules.end()) return false;
	else return true;
}
