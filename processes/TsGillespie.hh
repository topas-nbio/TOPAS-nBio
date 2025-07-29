#ifndef TsGillespie_hh
#define TsGillespie_hh

#include "globals.hh"
#include "TsVIRTProcedure.hh"
#include "TsIRTConfiguration.hh"
#include "G4VSolid.hh"

class TsParameterManager;
class TsGillespie{
public:
	TsGillespie(TsParameterManager*, G4String, TsVIRTProcedure*, G4double, G4double, G4int, G4double, G4double, G4double);
	~TsGillespie();

	void Initialize();
	void InitializeContainers();
	void RecoverReactionData();
	void Clean();
	void IncludeSpeciesAtTime();
	std::pair<G4int, G4double> Propensity();
	G4double TimeIncrement(G4double);
	G4double GetNextEscapeYieldTime();
	G4int GetMoleculeIndex(G4int);
	void AddMolecule(G4int,G4int,G4int);
	void AddMolecule(G4int);
	void UpdateContainers();
	void DoReaction(G4double, G4int);
	void RunAlt(G4double iniT = -1, G4double finT = -1);
	void DiffuseAllMolecules();
	void SampleNewMolecules();
	G4int GetSpeciesIndex() {return fSpeciesIndex;}
	void AddEscapeYields(std::map<G4int,G4int>,G4double);
	void AddDeltaYields(std::map<G4int,G4int>,G4double);

	// Debug Functions
	void PrintConcentrations();
	void PrintMoleculesAtTime();
	void PrintReactionInfo(G4int);
	void PrintIndividualPropensities();

	std::vector<G4double> GetStepTimes() {return fTimeSteps;}
	std::map<G4String, std::map<G4double, G4int>> GetGValues() {return fMoleculesAtTime;}
	std::map<G4int,std::map<G4double,G4int>> GetDeltaGValues() {return fDeltaGPerReactionPerTime;}
	std::unordered_map<G4int, TsIRTConfiguration::TsMolecule> GetChemicalSpecies() {return fMoleculesFromIRT;}

	void SetMolPerTime(G4int, G4double);

	void SetGValues(std::map<G4String, std::map<G4double, G4int>> gval) {fMoleculesAtTime = gval;}
	void SetDeltaGValues(std::map<G4int,std::map<G4double,G4int>> dgval) {fDeltaGPerReactionPerTime = dgval;}

	void SetCointanersStatus(G4bool isInit) {fContainersInitialized = isInit;}
	void UpdateConcentrations();
	void SaveState();

private:
	TsParameterManager* fPm;
	TsVIRTProcedure*    fIRT;

	G4String fName;

	G4int fSpeciesIndex;
	G4int fNoneIndex;
	G4long fTotalNumberOfMolecules;
	G4int fLastReactionIndex;
	G4int fNumberOfSavedStates;
	G4int fNumberOfStatedToEquilibrium;
	G4int fTotalNumberOfMoleculesSum1;
	G4int fTotalNumberOfMoleculesSum2;

	G4double fPercentageDifference;
	G4double fCubicVolume;
	G4double fInitialTime;
	G4double fFinalTime;
	G4double fInitialTimeInSeconds;
	G4double fFinalTimeInSeconds;
	G4double fCurrentTime;

	G4double fpVolume;

	std::map<G4int, G4double> fMolPerTime;

	std::vector<G4double>    fTimeSteps;
	std::map<G4int,G4int>    fConcentrations;
	std::map<G4int,G4int>    fConcentrationsSum1;
	std::map<G4int,G4int>    fConcentrationsSum2;
	std::map<G4int,G4double> fTotalMolecules;
	std::map<G4int,G4int>    fDeltas;
	std::map<G4int,G4String> fMoleculeNames;
	std::map<G4int,std::map<G4double, G4int>>    fDeltaGPerReactionPerTime;
	std::map<G4String,std::map<G4double, G4int>> fMoleculesAtTime;
	std::map<G4int,TsIRTConfiguration::TsMolecularReaction>   fReactions;
	std::unordered_map<G4int, TsIRTConfiguration::TsMolecule> fMoleculesFromIRT;

	// Alternative Gillespie Yield Load
	std::vector<G4double> fPulseTimes;
	std::vector<std::map<G4int,G4int>> fEscapeYieldsAtTime;
	std::map<G4double,std::map<G4int,G4int>> fDeltaYieldsAtTime;
	std::vector<G4double> prop;
	std::vector<G4int> equilibrium;

	TsIRTUtils* fUtils;
	G4bool fContainersInitialized;
	G4bool fPrintConcentrations;

	std::map<G4int, std::vector<G4int>> fMoleculeCanReactWith;
	G4double fPropensity;

	G4double fDoserate;
	G4double fInsertPropensity;
	G4double fIrradiationTime;
	std::vector<G4int> fReverse;
};

#endif
