#ifndef TsIRTConfiguration_hh
#define TsIRTConfiguration_hh

#include "globals.hh"
#include "G4ThreeVector.hh"
#include "Randomize.hh"

#include <vector>
#include <map>
#include <unordered_map>

class TsParameterManager;
class TsIRTUtils;

class TsIRTConfiguration {
public:
	TsIRTConfiguration(G4String, TsParameterManager*);
	
	~TsIRTConfiguration();
	
	void AddMolecule(G4String name, G4int moleculeID, G4double diffusionCoefficient, G4double charge, G4double radius);
	
	void AddMolecule(G4String name, G4double diffusionCoefficient, G4double charge, G4double radius);
	
	void AddMolecule(G4String name);
	
	void AdjustReactionRateForPH(G4String);

	void AdjustReactionAndDiffusionRateForTemperature();
	
	G4bool MoleculeExists(G4String name);
	
	G4double GetMoleculeRadius(G4int);
	G4int GetMoleculeCharge(G4int);
	
	G4int GetReactionIndex(G4int pdgA, G4int pdgB);
	
	G4double GetOnsagerRadius(G4int molA, G4int molB);
	
	void ResolveReactionParameters(G4int molA, G4int molB, G4double kobs, G4int reactionType);
	
	void ResolveReactionRateCoefficients();
	
	void CalculateContactProbabilities();
	
	void ResolveRemainerReactionParameters();
	
	void InsertReaction(G4String A, G4String B, G4String p1, G4String p2, G4String p3,
						G4double kobs, G4int reactionType);

	void InsertReaction(G4String A, G4String B, std::vector<G4String> p,
						G4double kobs, G4int reactionType);
	
	void InsertReaction(G4int molA, G4int molB, std::vector<G4int> products,
						G4double kobs, G4int reactionType);
	
	void InsertBackgroundReaction(G4String A, G4String B, G4String p1, G4String p2, G4String p3,
								  G4double kobs, G4double concentration, G4bool sampleExponential);

	void InsertBackgroundReaction(G4String A, G4String B, std::vector<G4String> p,
								  G4double kobs, G4double concentration, G4bool sampleExponential);
	
	void InsertBackgroundReaction(G4String A, G4String B, G4String p1, G4String p2, G4String p3,
								  G4double scavengingCapacity, G4bool sampleExponential);

	void InsertBackgroundReaction(G4String A, G4String B, std::vector<G4String> p,
								  G4double scavengingCapacity, G4bool sampleExponential);
	
	void QuitIfMoleculeNotFound(G4String mol);
	
	void Quit(const G4String& name, G4String message);
	
	void SetTimeLimits(G4double, G4double);
	
	G4int GetNumberOfReactions();
	
	std::pair<G4String, G4String> GetReactants(G4int);
	
	std::vector<G4String> GetProducts(G4int);
	
	// PH adjust functions
	
	std::vector<G4double> GetH2SO4ComponentsConcentrationP(G4double);
	
	std::vector<G4double> GetH2SO4ComponentsConcentrationPH(G4double);
	
	G4double IonicRate(G4double, G4double, G4int, G4int);
	
	G4double GetIonicStrength(std::vector<G4double>);

	inline std::map<G4String, G4int> GetMoleculeIDs() {return fMoleculesID;};
	
	inline std::map<G4int, G4String> GetMoleculeNames() { return fMoleculesName;};
	
private:
	TsParameterManager* fPm;
	TsIRTUtils* fUtils;
	G4String fName;
	
	struct TsMoleculeDefinition {
		G4double diffusionCoefficient;
		G4double charge;
		G4double radius;
	};
	
public:
	struct TsMolecularReaction {
		G4int    index;
		
		G4int reactorA;
		G4int reactorB;
		std::vector<G4int> products;
		
		G4double kobs;
		G4double kdif;
		G4double kact;
		G4double reactionRadius;
		G4double effectiveReactionRadius;
		G4double effectiveTildeReactionRadius;
		G4double probabilityOfReaction;
		G4double alpha;
		G4int    reactionType;
		G4bool   sampleExponential;
		G4double OnsagerRadius;
		G4double concentration;
		G4double scavengingCapacity;
	};

	struct TsMolecule {
		G4int id = -1;
		G4int trackID = -1;
		G4int spin = -1;

		G4double time = 0;
		G4ThreeVector position = G4ThreeVector();

		G4bool reacted = false;
		G4bool isDNA = false;
		G4bool isNew = true;
		//std::vector<G4int> tested;
	};
	
private:
	std::map<G4int, TsMolecule> fMolecules;
	
	std::map<G4int, TsMoleculeDefinition> fMoleculesDefinition;
	std::map<G4String, G4int> fMoleculesID;
	std::map<G4int, G4String> fMoleculesName;
	
	std::map<G4int, TsMolecularReaction > fReactions;
	std::map<G4int, std::vector<std::pair<G4int,G4int>>> fMoleculeCanReactWith;
	
	std::map<G4String, G4String> fExistingMolecules;
	
	G4double fUpperTime;
	G4double fLowerTime;
	G4int fReactionID;
	G4int fTotalBinaryReaction;
	G4int fLastMoleculeID;
	
	G4double fTemperature;
	G4bool fScaleForTemperature;
	G4bool fKick;
	
	G4String fpHSolvent;
	G4double fpHSolventConcentration;
	G4double fpHValue;
	
	G4bool fQualityAssurance;
	
	G4bool fAllTotallyDiffusionControlled;
	
public:
	
	G4double SampleExponentialTime(G4int pdgA, G4int pdgB, G4int indexOfReaction);
	G4double GetIndependentReactionTime(TsMolecule molA, TsMolecule molB, G4int indexOfReaction);
	G4double SampleIRTTotallyDiffusionControlled(TsMolecule molA, TsMolecule molB, G4int indexOfReaction);
	G4double SampleIRTPartiallyDiffusionControlled(TsMolecule molA, TsMolecule molB, G4int indexOfReaction);
	std::pair<G4int, G4double> SampleIRTFirstOrderAndBackgroundReactions(TsMolecule molA );
	std::vector<std::pair<G4int, G4double>> SampleAllIRTFirstOrderAndBackgroundReactions(TsMolecule molA );
	
	G4int ContactFirstOrderAndBackgroundReactions(TsMolecule molA );
	
	G4double CalculateProbabilityPartiallyDiffusionControlled(TsMolecule molA, TsMolecule molB, G4int indexOfReaction, G4double t);
	G4double CalculateProbabilityOfScavenger(TsMolecule molA, G4int indexOfReaction, G4double t);
	
	G4double brents_fun(TsMolecule molA, TsMolecule molB, G4int indexOfReaction, G4double offset);
	G4double brents_fun_scav(TsMolecule molA, G4int indexOfReaction, G4double offset);
	G4double SolveTime(TsMolecule molA, G4int indexOfReaction, G4double offset);
	
	void ResampleReactantsPosition(TsMolecule& molA, TsMolecule& molB, G4int index, G4double time);
	std::vector<G4ThreeVector> GetPositionOfProducts(TsMolecule molA, TsMolecule molB, G4int index);
	G4double GetRCutOff(G4double tCutOff);
	
	G4bool MakeReaction(std::unordered_map<G4int,TsMolecule> &initialSpecies, G4int& speciesIndex,
						std::unordered_map<G4int, std::unordered_map<G4int, std::unordered_map<G4int, std::unordered_map<G4int,G4bool>>>> &spaceBinned,
						G4int NX, G4int NY, G4int NZ, G4double XMin, G4double XMax, G4double YMin, G4double YMax, G4double ZMin, G4double ZMax,
						std::map<G4int, std::map<G4int, G4int>> &theGvalue, std::vector<G4double> timeSteps,
						G4int iM, G4int indexOfReaction, G4double irt, std::unordered_map<G4int,G4bool> &used, std::vector<G4int>& prods);
	
	G4bool MakeReaction(std::unordered_map<G4int,TsMolecule> &initialSpecies, G4int& speciesIndex,
						std::unordered_map<G4int, std::unordered_map<G4int, std::unordered_map<G4int, std::unordered_map<G4int,G4bool>>>> &spaceBinned,
						G4int NX, G4int NY, G4int NZ, G4double XMin, G4double XMax, G4double YMin, G4double YMax, G4double ZMin, G4double ZMax,
						std::map<G4int, std::map<G4int, G4int>> &theGvalue, std::vector<G4double> timeSteps,
						G4int iM, G4int jM, G4int indexOfReaction, G4double irt, G4double probabilityOfReaction, std::unordered_map<G4int,G4bool> &used, std::vector<G4int>& prods);
	
	G4bool MakeReaction(std::unordered_map<G4int,TsMolecule> &initialSpecies, G4int& speciesIndex,
						std::unordered_map<G4int, std::unordered_map<G4int, std::unordered_map<G4int, std::unordered_map<G4int,G4bool>>>> &spaceBinned,
						G4int NX, G4int NY, G4int NZ, G4double XMin, G4double XMax, G4double YMin, G4double YMax, G4double ZMin, G4double ZMax,
						std::map<G4int, std::map<G4int, G4int>> &theGvalue,
						std::map<G4int, std::map<G4int, G4int>> &theGvalueInVolume, std::vector<G4double> timeSteps,
						G4int iM, G4int jM, G4int indexOfReaction, G4double irt, G4double probabilityOfReaction, std::unordered_map<G4int,G4bool> &used, std::vector<G4int>& prods);

	G4bool Inside(G4ThreeVector p);
	
	void ScoreGvalue(std::vector<TsMolecule> &initialSpecies,
					 std::map<G4int, std::map<G4int, G4int>> &theGvalueInVolume,
					 std::vector<G4double> timeSteps,
					 G4int iM, G4int jM, G4int indexOfReaction, G4double irt);
	
	TsMolecularReaction GetReaction(G4int index);
	std::map<G4int, G4String> GetMoleculeName() {return fMoleculesName;}
	std::map<G4int, TsIRTConfiguration::TsMolecularReaction> GetReactions() {return fReactions;}
	
	void Diffuse(TsMolecule& mol, G4double dt);
	
	
	void TestSampling(G4int indexOfReaction, G4int nHistories);
	
	void PrintMoleculesInformation();
	void PrintReactionsInformation();
	
	G4int GetLastMoleculeID() { return fLastMoleculeID; };
	G4int GetLastReactionID() { return fReactionID-1; };
	
	G4double IonicRate(G4double, TsMolecularReaction);
	
};
#endif

