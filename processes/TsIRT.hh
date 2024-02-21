#ifndef TsIRT_hh
#define TsIRT_hh

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
class TsIRT {
public:
	TsIRT(TsParameterManager*, G4String);
	~TsIRT();
	
	void runIRT(G4double startTime = -1, G4double finalTime = -1);
	
	void AddMolecule(TsIRTConfiguration::TsMolecule);
	void AddMolecule(G4Step*, G4double, G4int, G4ThreeVector);
	void AddMolecule(G4Track*, G4double, G4int, G4ThreeVector);
	void AddMolecule(const G4Track*, G4double, G4int, G4ThreeVector);
	void AddMolecule(G4Track*, G4double, G4int, G4ThreeVector, G4bool);
	TsIRTConfiguration::TsMolecule ConstructMolecule(G4Track*, G4double, G4int, G4ThreeVector);
	
	void Clean();
	
	G4bool Inside(G4ThreeVector);
	
	inline std::vector<G4double> GetStepTimes() {return fStepTimes; };
	inline TsIRTUtils* GetUtils() { return fUtils; };
	
	void SetGValues(std::map<G4String, std::map<G4double, G4int>>);
	void SetDeltaGValues(std::map<G4int, std::map<G4double, G4int>>);
	std::map<G4String, std::map<G4double, G4int>> GetGValues();
	std::map<G4int, std::map<G4double, G4int>>    GetDeltaGValues();
	std::map<G4String, std::map<G4double, G4int>> GetGValuesInVolume();
	std::map<G4int, std::pair<G4int,G4int>>       GetReactedDNA();
	
	std::pair<G4String, G4String> GetReactants(G4int);
	
	std::vector<G4String> GetProducts(G4int);

	std::vector<G4int> GetSurvivingMolecules();

	inline TsIRTConfiguration* GetIRTConfiguration() {return fReactionConf;}
	inline std::unordered_map<G4int, TsIRTConfiguration::TsMolecule> GetChemicalSpecies() {return fChemicalSpecies;}
	G4int GetSpeciesIndex() {return fSpeciesIndex;}
	void SetCointanersStatus(G4bool isInit) {fScorersInitialized = isInit;}
	std::map<G4int,G4int> GetConcentrations() {return fConcentrations;}
	void PrintConcentrations();
	std::map<G4int,G4int> GetEscapeYields();
	std::map<G4int,G4int> GetDeltaYields();

private:
	void initializeScorers();
	void FindBinIndexes(G4ThreeVector thisPos, G4double rcutOff);
	void contactReactions(G4int i,std::unordered_map<G4int, G4bool> j);
	void sampleReactions(G4int i);
	G4String GetFullParmName(G4String name);

	void VoxelizeAndSortSpace();
	void TestForContactReactions();
	void SampleIndependantReactionTimes();
	void SortIndependantReactionTimes();
	void ConductReactions();
	void CleanIRTVariables();
	void UpdateGValues();
	void CleanReactedMolecules();
	void AddToIRT(G4double Time, G4int Reaction, G4int molA, G4int molB, G4double OrigTime, G4ThreeVector aPos, G4ThreeVector bPos, G4bool isBack);
	void AddIRTinAscendantOrder(G4double Time);
	void RemoveFirstIRTElement();
	void RemoveMolecule(G4int Index);
	G4int CountSurvivingMolecules();
	G4bool MoleculeExists(G4int Index);

private:
	
	TsParameterManager* fPm;
	TsIRTConfiguration* fReactionConf;
	TsIRTUtils* fUtils;
	
	G4String fName;
	
	std::unordered_map<G4int,TsIRTConfiguration::TsMolecule> fChemicalSpecies;
	
	std::vector<G4double> fStepTimes;
	std::map<G4int, G4String> fMolecules;
	std::map<G4String, G4int> fMoleculesIDs;
	std::map<G4String, std::map<G4double, G4int>> fGValues;
	std::map<G4int, std::map<G4double, G4int>> fDeltaGValues;
	std::map<G4String, std::map<G4double, G4int>> fGValuesInVolume;
	std::map<G4int, std::pair<G4int,G4int>> fReactedDNA;
	
	G4bool fUseSpinScaled;
	G4bool fHighTimeScavenger;
	G4bool fReportDelta;
	G4double fRCutOff;
	G4double fBinWidth;
	G4double fXMin;
	G4double fXMax;
	G4double fYMin;
	G4double fYMax;
	G4double fZMin;
	G4double fZMax;
	G4int fNx;
	G4int fNy;
	G4int fNz;
	G4int fVerbosity;
	G4double fDx;
	G4double fDy;
	G4double fDz;
	
	G4int fxiniIndex;
	G4int fxendIndex;
	G4int fyiniIndex;
	G4int fyendIndex;
	G4int fziniIndex;
	G4int fzendIndex;
	
	G4bool fTestForContactReactions;
	G4bool fSortByTime;
	G4bool fTestIsInside;
	G4bool fScorersInitialized;
	G4bool fReactedByContact;
	G4bool fSampleIRTatStart;
	
	std::map<G4int,G4int>    fConcentrations;
	std::map<G4int, std::map<G4int, G4int>> fTheGvalue;
	std::map<G4int, std::map<G4int, G4int>> fTheGvalueInVolume;
	std::unordered_map<G4int,std::unordered_map<G4int,std::unordered_map<G4int,std::unordered_map<G4int,G4bool>>>> fSpaceBinned;
	std::unordered_map<G4int,G4bool> fUsed;
	std::vector<G4int> fContactProducts;

	G4int fChunk;
	G4int fIndex;
	G4int fSpeciesIndex;
	G4double fMinTime;
	G4double fTimeCut;
	G4double fCurrentTimeScale;
	std::vector<std::pair<G4double, G4int > > fIRTValues;
	std::unordered_map<G4int,G4int> fIRTIndex;
	std::unordered_map<G4int,G4int> fIRTMolA;
	std::unordered_map<G4int,G4int> fIRTMolB;
	std::unordered_map<G4int,G4bool> fIRTIsBackground;
	std::unordered_map<G4int,G4double> fIRTOrigTime;
	std::unordered_map<G4int,std::pair<G4ThreeVector, G4ThreeVector> > fIRTPositions;
	
};
#endif

