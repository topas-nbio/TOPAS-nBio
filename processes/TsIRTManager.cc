// Extra Class for TsEmDNAChemistry
#include "TsIRTManager.hh"
#include "TsIRT.hh"
#include "TsIRTContinuous.hh"
#include "TsParameterManager.hh"

TsIRTManager::TsIRTManager(TsParameterManager* pM, G4String parmName)
: fPm(pM), fName(parmName) {
	G4String chemistryList = fPm->GetStringParameter("Ch/ChemistryName");

	G4String parName = "Ch/" + chemistryList + "/IRTProcedure";
	
	G4String IRTType = "pure";
	fIRTType = 0;
	if ( fPm->ParameterExists(parName) )
		IRTType = fPm->GetStringParameter(parName);

	G4StrUtil::to_lower(IRTType);
	if (IRTType == "pure")
		fIRTProcedure = new TsIRT(fPm,fName);
	else if (IRTType == "continuous"){                 // This is Wook-Geun Implementation
		fIRTProcedure = new TsIRTContinuous(fPm,fName);
		fIRTType      = 1;
	}
	else {
		G4cout << "- Error in TOPAS-nBio IRT Manager!" << G4endl;
		G4cout << "   There is no " << IRTType << " IRT Type" << G4endl;
		G4cout << "   Available options are: pure, hybrid, continuous" << G4endl;
		G4cout << "   Please refer to the IRT section of the user documentation" << G4endl;
		exit(1);
	}
}

TsIRTManager::~TsIRTManager() {
	delete fIRTProcedure;
}

void TsIRTManager::runIRT(G4double startTime, G4double finalTime, G4double transTime, G4bool isContinuation) {
	if (fIRTType == 0)
		fIRTProcedure->runIRT();

	else
		fIRTProcedure->runIRT(startTime,finalTime,transTime,isContinuation);
}

std::pair<G4String, G4String> TsIRTManager::GetReactants(G4int reactIndex) {
	return fIRTProcedure->GetReactants(reactIndex);
}

std::vector<G4String> TsIRTManager::GetProducts(G4int reactIndex) {
	return fIRTProcedure->GetProducts(reactIndex);
}

// temporal function
void TsIRTManager::AddMolecule(G4int molID, G4ThreeVector position, G4double time,
                               G4int trackID, G4bool isDNA, G4int volumeID, G4int baseID, G4int strandID) {
	fIRTProcedure->AddMolecule(molID, position, time, trackID, isDNA, volumeID, baseID, strandID);
}

void TsIRTManager::AddMolecule(TsIRTConfiguration::TsMolecule aMol) {
	fIRTProcedure->AddMolecule(aMol);
}

void TsIRTManager::AddMolecule(G4Track* aTrack, G4double time, G4int moleculeID, G4ThreeVector offset) {
	fIRTProcedure->AddMolecule(aTrack, time, moleculeID, offset);
}

void TsIRTManager::AddMolecule(const G4Track* aTrack, G4double time, G4int moleculeID, G4ThreeVector offset) {
	fIRTProcedure->AddMolecule(aTrack, time, moleculeID, offset);
}

void TsIRTManager::AddMolecule(G4Track* aTrack, G4double time, G4int moleculeID, G4ThreeVector offset, G4bool isDNA) {
	fIRTProcedure->AddMolecule(aTrack, time, moleculeID, offset, isDNA);
}

void TsIRTManager::Clean() {
	fIRTProcedure->Clean();
}

void TsIRTManager::SetContainersForNextPulse() {
	fIRTProcedure->SetContainersForNextPulse();
}

TsIRTUtils* TsIRTManager::GetUtils() {
	return fIRTProcedure->GetUtils();
}


TsVIRTProcedure* TsIRTManager::GetIRTProcedure() {
	return fIRTProcedure;
}


std::vector<G4double> TsIRTManager::GetStepTimes() {
	return fIRTProcedure->GetStepTimes();
}

TsIRTConfiguration::TsMolecule TsIRTManager::ConstructMolecule(G4Track* aTrack, G4double time, G4int moleculeID, G4ThreeVector offset) {
	return fIRTProcedure->ConstructMolecule(aTrack, time, moleculeID, offset);
}

TsIRTConfiguration* TsIRTManager::GetIRTConfiguration() {
	return fIRTProcedure->GetIRTConfiguration();
}

std::map<G4String, std::map<G4double, G4int>> TsIRTManager::GetGValues() {
	return fIRTProcedure->GetGValues();
}

std::map<G4int, std::map<G4double, G4int>>    TsIRTManager::GetDeltaGValues() {
	return fIRTProcedure->GetDeltaGValues();
}

std::vector<TsIRTConfiguration::TsMolecule> TsIRTManager::GetSurvivingMoleculesWithMolID(G4int molID) {
	return fIRTProcedure->GetSurvivingMoleculesWithMolID(molID);
}

std::map<G4int, std::pair<G4int,G4int>>       TsIRTManager::GetReactedDNA() {
	return fIRTProcedure->GetReactedDNA();
}
