// Scorer for MoleculePhaseSpace
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

#include "TsScoreMoleculePhaseSpace.hh"
#include "G4ProcessManager.hh"
#include "G4ParticleDefinition.hh"
#include "G4VProcess.hh"
#include "G4RunManager.hh"
#include "G4Run.hh"

#include "G4SystemOfUnits.hh"
#include "G4SteppingManager.hh"
#include "G4VTouchable.hh"
#include "G4VPhysicalVolume.hh"

#include "G4Molecule.hh"

TsScoreMoleculePhaseSpace::TsScoreMoleculePhaseSpace(TsParameterManager* pM, TsMaterialManager* mM, TsGeometryManager* gM, TsScoringManager* scM, TsExtensionManager* eM,
										   G4String scorerName, G4String quantity, G4String outFileName, G4bool isSubScorer)
:TsVNtupleScorer(pM, mM, gM, scM, eM, scorerName, quantity, outFileName, isSubScorer)
{
	SetUnit("");
	
	fNtuple->RegisterColumnS(&fParticleName, "Particle name");
	fNtuple->RegisterColumnF(&fPosX, "Position X", "um");
	fNtuple->RegisterColumnF(&fPosY, "Position Y", "um");
	fNtuple->RegisterColumnF(&fPosZ, "Position Z", "um");
	fNtuple->RegisterColumnF(&fTime, "Global time", "ps");
	fNtuple->RegisterColumnI(&fEvt, "Event ID");
	
	fTimeCut = 1.0 * us;
	if ( fPm->ParameterExists(GetFullParmName("TimeCut") ) )
		fTimeCut = fPm->GetDoubleParameter(GetFullParmName("TimeCut"),"Time");
}


TsScoreMoleculePhaseSpace::~TsScoreMoleculePhaseSpace() {;}


G4bool TsScoreMoleculePhaseSpace::ProcessHits(G4Step* aStep,G4TouchableHistory*)
{
	if (!fIsActive) {
		fSkippedWhileInactive++;
		return false;
	}
	
	G4Track* aTrack = aStep->GetTrack();
	
	if (aTrack->GetTrackID() < 0)
	{
		fTime = aStep->GetPreStepPoint()->GetGlobalTime();
		fEvt = GetEventID();
		G4ThreeVector pos = aStep->GetPreStepPoint()->GetPosition();
		fPosX = pos.x();
		fPosY = pos.y();
		fPosZ = pos.z();
		fParticleName = GetMolecule(aTrack)->GetName();
		fNtuple->Fill();
		if (fTime >= fTimeCut)
		{
			aStep->GetTrack()->SetTrackStatus(fStopAndKill);
			return true;
		}
		return true;
	}
	return false;
}

