// Scorer for TrackStructure
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

#include "TsScoreTrackStructure.hh"

TsScoreTrackStructure::TsScoreTrackStructure(TsParameterManager* pM, TsMaterialManager* mM, TsGeometryManager* gM, TsScoringManager* scM, TsExtensionManager* eM,
											 G4String scorerName, G4String quantity, G4String outFileName, G4bool isSubScorer)
					  :TsVNtupleScorer(pM, mM, gM, scM, eM, scorerName, quantity, outFileName, isSubScorer),
fIncludeEventID(true),  fIncludeTrackID(true), fIncludeParentID(false), fIncludeParticleName(false), fIncludeVolumeName(false), fIncludePrimaryPositions(true)
{
	SetUnit("");

	fNtuple->RegisterColumnF(&fPosX, "Position X", "um");
	fNtuple->RegisterColumnF(&fPosY, "Position Y", "um");
	fNtuple->RegisterColumnF(&fPosZ, "Position Z", "um");
	fNtuple->RegisterColumnF(&fEdep, "Energy deposited", "keV");

	if (fPm->ParameterExists(GetFullParmName("IncludePrimaryPositions")))
		fIncludePrimaryPositions = fPm->GetBooleanParameter(GetFullParmName("IncludePrimaryPositions"));

	if (fPm->ParameterExists(GetFullParmName("IncludeEventID")))
		fIncludeEventID = fPm->GetBooleanParameter(GetFullParmName("IncludeEventID"));

	if (fPm->ParameterExists(GetFullParmName("IncludeTrackID")))
		fIncludeTrackID = fPm->GetBooleanParameter(GetFullParmName("IncludeTrackID"));

	if (fPm->ParameterExists(GetFullParmName("IncludeParentID")))
		fIncludeParentID = fPm->GetBooleanParameter(GetFullParmName("IncludeParentID"));

	if (fPm->ParameterExists(GetFullParmName("IncludeParticleName")))
		fIncludeParticleName = fPm->GetBooleanParameter(GetFullParmName("IncludeParticleName"));

	if (fPm->ParameterExists(GetFullParmName("IncludeVolumeName")))
		fIncludeVolumeName = fPm->GetBooleanParameter(GetFullParmName("IncludeVolumeName"));

	if (fIncludeEventID)
		fNtuple->RegisterColumnI(&fEvt, "EventID");

	if (fIncludeTrackID)
		fNtuple->RegisterColumnI(&fTrackID, "Track ID");

	if (fIncludeParticleName)
		fNtuple->RegisterColumnS(&fParticleName, "Particle name");

	if (fIncludeVolumeName)
		fNtuple->RegisterColumnS(&fVolumeName, "Volume name");

	if (fIncludeParentID)
		fNtuple->RegisterColumnI(&fParentID, "Parent ID");
}

TsScoreTrackStructure::~TsScoreTrackStructure() { }

G4bool TsScoreTrackStructure::ProcessHits(G4Step* aStep, G4TouchableHistory*)
{
	if (!fIsActive)
	{
		fSkippedWhileInactive++;
		return false;
	}

	fEdep = aStep->GetTotalEnergyDeposit();
	G4Track* aTrack = aStep->GetTrack();

	if (fEdep > 0 || (fIncludePrimaryPositions && aTrack->GetTrackID() == 0))
	{

		fEvt = GetEventID();
		G4ThreeVector pos = aStep->GetPreStepPoint()->GetPosition();
		fPosX = pos.x();
		fPosY = pos.y();
		fPosZ = pos.z();


		fTrackID = aTrack->GetTrackID();
		fParentID = aTrack->GetParentID();
		fParticleName = aTrack->GetParticleDefinition()->GetParticleName();

		G4TouchableHistory* touchable = (G4TouchableHistory*)(aStep->GetPreStepPoint()->GetTouchable());
		fVolumeName = touchable->GetVolume()->GetName();

		fNtuple->Fill();

		return true;
	}
	else
		return false;
}

