#ifndef TsHybridIRT_hh
#define TsHybridIRT_hh

#include "TsIRTConfiguration.hh"
#include "TsVIRTProcedure.hh"

#include "globals.hh"
#include "G4ThreeVector.hh"
#include "G4Step.hh"
#include "G4VSolid.hh"

#include <vector>
#include <map>
#include <unordered_map>

class TsParameterManager;
class TsIRTUtils;
class TsHybridIRT : public TsVIRTProcedure {
public:
	TsHybridIRT(TsParameterManager*, G4String);
	~TsHybridIRT();

	void runIRT(G4double startTime = -1, G4double finalTime = -1, G4double transTime=-1, G4bool isContinuation=false) override;
	void Clean() override;

	//Prepare the remaining species for next Pulse
	void SetContainersForNextPulse() override;

protected:
	void contactReactions(G4int i,std::unordered_map<G4int, G4bool> j) override;
	void sampleReactions(G4int i) override;
	void ConductReactions() override;
	void CleanIRTVariables() override;
	void RemoveMolecule(G4int Index) override;
	G4bool IsNextIRTStillPossible();

	//Homogeneous Chemistry Methods
	// Gillespie Direct Method
	void UpdateGillespieDirect(G4double);
	G4bool DoGillespieDirect(G4double,G4bool);
	std::pair<G4int,G4int> ChooseReactantsBasedOnDistance(std::vector<G4int>,std::vector<G4int>,G4int,G4int,G4int,G4double);
	G4int ChooseReactantRandomly(std::vector<G4int>,G4int,G4double);

	// Adaptative Tau Leaping Method; Yang Cao (2006) 
	void UpdateGillespieTauLeaping(G4double);
	G4bool DoGillespieTauLeaping(G4double);
	std::vector<G4int> GetAllSpeciesOfAKindAtATime(G4double,G4int);

	void GetGillespieDirectReactants(G4double currentTime);
	void GetGillespieDirectTime();

protected:
	TsParameterManager* fPm;
	G4String fName;

	G4int fVerbosity;

	G4bool fHighTimeScavenger;
	G4bool fScorersInitialized;

	// Transition Time Variables
	G4double fTransCut;

	// Homogeneous Chemistry containers
	G4bool fUseVariableBackground;
	G4double fCubicVolume;
	std::unordered_map<G4int,G4double> fPropensityAtThisTime;
	std::unordered_map<G4int,std::vector<G4int>> fConcentrationsAtThisTime;
	std::unordered_map<G4int,G4int> fHomogeneousConcentrations;
	std::unordered_map<G4int,G4int> fInitialBackgroundConcentrations;
	std::unordered_map<G4int,G4int> fCurrentBackgroundConcentrations;

	// Gillespie Direct Method
	G4int    fHomogeneousReactIndex;
	G4double fHomogeneousTimeStep;
	G4double fTotalPropensityAtThisTime;
	G4bool fGillespieFinalTau;

	G4int fGillespieIndexA;
	G4int fGillespieIndexB;

	// Adaptative Tau Leaping Method
	std::unordered_map<G4int,G4double> fMuForMolecule;
	std::unordered_map<G4int,G4double> fSigmaForMolecule;
	std::unordered_map<G4int,G4double> fHORForMolecule;
	std::unordered_map<G4int,G4double> fGFactorForMolecule;
	std::unordered_map<G4int,G4double> fMinForMolecule;
	std::unordered_map<G4int,G4double> fMaxForMolecule;
	std::unordered_map<G4int,G4double> fLimitsForMolecule;
	std::unordered_map<G4int,G4int>    fChangesForMolecule;
	G4int fTotalChanges;

};
#endif

