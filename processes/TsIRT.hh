#ifndef TsIRT_hh
#define TsIRT_hh

#include "TsVIRTProcedure.hh"
#include "TsIRTConfiguration.hh"
#include "globals.hh"
#include "G4ThreeVector.hh"
#include "G4Step.hh"
#include "G4VSolid.hh"

#include <vector>
#include <map>
#include <unordered_map>

class TsParameterManager;
class TsIRTUtils;

class TsIRT: public TsVIRTProcedure {
public:
	TsIRT(TsParameterManager*, G4String);
	~TsIRT();
	
	void runIRT(G4double startTime = -1, G4double finalTime = -1, G4double transTime=-1, G4bool isContinuation=false) override;
	void Clean() override;

	void SetContainersForNextPulse() override;

protected:
	void CleanIRTVariables() override;
	void ConductReactions() override;

	void contactReactions(G4int i,std::unordered_map<G4int, G4bool> used) override;

	void RemoveMolecule(G4int Index) override;

protected:
	TsParameterManager* fPm;
	G4String fName;

	G4int fVerbosity;

	G4bool fHighTimeScavenger;
	G4bool fScorersInitialized;
	
};
#endif

