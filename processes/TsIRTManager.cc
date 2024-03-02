// Extra Class for TsEmDNAChemistry
#include "TsIRTManager.hh"
#include "TsIRT.hh"
#include "TsHybridIRT.hh"
#include "TsParameterManager.hh"

TsIRTManager::TsIRTManager(TsParameterManager* pM, G4String parmName)
: fPm(pM), fName(parmName) {
	G4String chemistryList = fPm->GetStringParameter("Ch/ChemistryName");

	G4String parName = "Ch/" + chemistryList + "/IRTProcedure";
	G4String IRTType = "pure";
	if ( fPm->ParameterExists(parName) )
		IRTType = fPm->GetStringParameter(parName);

	IRTType.toLower();
	if (IRTType == "pure")
		fIRTProcedure = new TsIRT(fPm,fName);
	else if (IRTType == "hybrid")
		fIRTProcedure = new TsHybridIRT(fPm,fName);
	//else if (IRTType == "continuous")
	//	fIRTProcedure = new TsIRTContinuous(fPm,fName);
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

void TsIRTManager::runIRT() {
	fIRTProcedure->runIRT();
}

std::pair<G4String, G4String> TsIRTManager::GetReactants(G4int reactIndex) {
	return fIRTProcedure->GetReactants(reactIndex);
}

std::vector<G4String> TsIRTManager::GetProducts(G4int reactIndex) {
	return fIRTProcedure->GetProducts(reactIndex);
}

void TsIRTManager::AddMolecule(G4Track* aTrack, G4double time, G4int moleculeID, G4ThreeVector offset) {
	fIRTProcedure->AddMolecule(aTrack, time, moleculeID, offset);
}

void TsIRTManager::AddMolecule(G4Track* aTrack, G4double time, G4int moleculeID, G4ThreeVector offset, G4bool isDNA) {
	fIRTProcedure->AddMolecule(aTrack, time, moleculeID, offset, isDNA);
}

void TsIRTManager::Clean() {
	fIRTProcedure->Clean();
}

std::map<G4String, std::map<G4double, G4int>> TsIRTManager::GetGValues() {
	return fIRTProcedure->GetGValues();
}
std::map<G4int, std::map<G4double, G4int>>    TsIRTManager::GetDeltaGValues() {
	return fIRTProcedure->GetDeltaGValues();
}