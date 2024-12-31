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
	void AddMolecule(G4Track*, G4double, G4int, G4ThreeVector);
	void AddMolecule(G4int, G4ThreeVector, G4double,
                               G4int, G4bool, G4int, G4int, G4int);

	TsIRTConfiguration::TsMolecule ConstructMolecule(G4Track*, G4double, G4int, G4ThreeVector);

	void SetContainersForNextPulse() override;

protected:
	void CleanIRTVariables() override;
	void ConductReactions() override;

	void Sampling();
	void sampleReactions(G4int i);

	void RemoveMolecule(G4int Index) override;
	G4bool MoleculeExists(G4int Index);

protected:

//	std::vector<TsIRTConfiguration::TsMolecule> fChemicalSpecies;
	std::unordered_map<G4int,TsIRTConfiguration::TsMolecule> fChemicalSpecies;
//	std::unordered_map<G4int,std::unordered_map<G4int,std::unordered_map<G4int,std::vector<G4int>>>> fSpaceBinned;
	std::unordered_map<G4int,std::unordered_map<G4int,std::unordered_map<G4int,std::unordered_map<G4int,G4bool>>>> fSpaceBinned;
	std::map<G4int, G4String> fMoleculesName;

	TsIRTConfiguration* fReactionConf;
	TsParameterManager* fPm;
	G4String fName;

	G4int fChemVerbosity;

	G4bool fHighTimeScavenger;
	
};
#endif

