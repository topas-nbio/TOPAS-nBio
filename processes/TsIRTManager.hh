#ifndef TsIRTManager_hh
#define TsIRTManager_hh

#include "globals.hh"
#include <map>
#include "G4Track.hh"
#include "G4ThreeVector.hh"

class TsVIRTProcedure;
class TsParameterManager;

class TsIRTManager {
public:
	TsIRTManager(TsParameterManager*, G4String);
	~TsIRTManager();

	void runIRT();
	void AddMolecule(G4Track*, G4double, G4int, G4ThreeVector);
	void AddMolecule(G4Track*, G4double, G4int, G4ThreeVector, G4bool);
	void Clean();

	std::pair<G4String, G4String> GetReactants(G4int);
	std::vector<G4String> GetProducts(G4int);

	std::map<G4String, std::map<G4double, G4int>> GetGValues();
	std::map<G4int, std::map<G4double, G4int>>    GetDeltaGValues();

private:
	G4String fName;
	TsVIRTProcedure* fIRTProcedure;
	TsParameterManager* fPm;
};

#endif