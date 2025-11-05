#ifndef TsIRTManager_hh
#define TsIRTManager_hh

#include "globals.hh"
#include <map>
#include "G4Track.hh"
#include "G4ThreeVector.hh"
#include "TsIRTUtils.hh"
#include "TsIRTConfiguration.hh"
#include "TsIRTConfiguration.hh"

class TsVIRTProcedure;
class TsParameterManager;

class TsIRTManager {
public:
	TsIRTManager(TsParameterManager*, G4String);
	~TsIRTManager();

	void runIRT(G4double startTime = -1, G4double finalTime = -1, G4double transTime=-1, G4bool isContinuation=false);
	void AddMolecule(TsIRTConfiguration::TsMolecule);
	void AddMolecule(G4Track*, G4double, G4int, G4ThreeVector);
	void AddMolecule(const G4Track*, G4double, G4int, G4ThreeVector);
	void AddMolecule(G4Track*, G4double, G4int, G4ThreeVector, G4bool);
        // temporal function
	void AddMolecule(G4int, G4ThreeVector, G4double,
                               G4int, G4bool, G4int, G4int, G4int);
	void Clean();
	void SetContainersForNextPulse();

	TsIRTUtils* GetUtils();
	TsVIRTProcedure* GetIRTProcedure();
	std::vector<G4double> GetStepTimes();
	TsIRTConfiguration::TsMolecule ConstructMolecule(G4Track*, G4double, G4int, G4ThreeVector);
	TsIRTConfiguration* GetIRTConfiguration();

	std::pair<G4String, G4String> GetReactants(G4int);
	std::vector<G4String> GetProducts(G4int);

	std::map<G4String, std::map<G4double, G4int>> GetGValues();
	std::map<G4int, std::map<G4double, G4int>>    GetDeltaGValues();
	std::map<G4int, std::pair<G4int,G4int>>       GetReactedDNA();

	std::vector<TsIRTConfiguration::TsMolecule> GetSurvivingMoleculesWithMolID(G4int);

private:
	G4String fName;
	TsVIRTProcedure* fIRTProcedure;
	TsParameterManager* fPm;

	G4int fIRTType;
};

#endif
