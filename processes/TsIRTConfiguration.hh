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

	void AddMolecule(G4String name, G4double diffusionCoefficient, G4double charge, G4double radius, G4int moleculeID=0);
	G4bool MoleculeExists(G4String name);
	void ResolveReactionParameters();
	void InsertReaction(G4String, G4String, std::vector<G4int>, G4double, G4int);
	void InsertBackgroundReaction(G4String, G4String, std::vector<G4int>,
								  G4double, G4double, G4bool);
	void SetTimeLimits(G4double, G4double);

	

	void QuitIfMoleculeNotFound(G4String mol);
	void Quit(const G4String& name, G4String message);


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
		G4double probabilityOfReaction;
		G4double diffusionCoefficient;
		G4int    reactionType;
		G4bool   sampleExponential;
		G4bool   positionSensitive=false;
		G4double OnsagerRadius;
		G4double concentration;
		G4double scavengingCapacity;

		// Multiple Channel Reaction
		G4bool isMultiChannel = false;
		std::vector<std::vector<G4int>> productsChannel;
		std::vector<G4double> weightsChannels;
	};

	struct TsMolecule {
		G4int id       = -1;
		G4int spin     = -1;

		G4double time  = 0;
		G4ThreeVector position = G4ThreeVector();

		G4bool isDNA   = false;
		G4bool reacted = false;
		G4bool isNew   = true;

		G4int trackID  = -1;
		G4int parentID = 0;

		G4int volumeID = -1;
		G4int baseID   = -1;
		G4int strandID = -1;

		G4int chemAlgo = -1; // To differentiate IRT from direct Gillespie
	};
	
private:
	std::map<G4int, TsMolecule> fMolecules;
	std::map<G4int, TsMoleculeDefinition> fMoleculesDefinition;
	std::map<G4String, G4int> fMoleculesID;
	std::map<G4int, G4String> fMoleculesName;
	std::map<G4String, G4String> fExistingMolecules;
	
	std::map<G4int, TsMolecularReaction > fReactions;
	std::map<G4int, std::vector<std::pair<G4int,G4int>>> fMoleculeCanReactWith;
	
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
	
	G4bool fAllTotallyDiffusionControlled;
	
public:
	
	G4double GetIndependentReactionTime(TsMolecule molA, TsMolecule molB, G4int indexOfReaction);
	std::pair<G4int, G4double> SampleIRTFirstOrderAndBackgroundReactions(TsMolecule molA );
	G4int ContactFirstOrderAndBackgroundReactions(TsMolecule molA );
	
	std::vector<G4ThreeVector> ResampleReactantsPosition(TsMolecule& molA, TsMolecule& molB, G4int index, G4double time);
	std::vector<G4ThreeVector> GetBackgroundPositionOfProducts(TsMolecule molA, G4int index);

	G4bool Inside(G4ThreeVector p);
	
	void ScoreGvalue(std::vector<TsMolecule> &initialSpecies,
					 std::map<G4int, std::map<G4int, G4int>> &theGvalueInVolume,
					 std::vector<G4double> timeSteps,
					 G4int iM, G4int jM, G4int indexOfReaction, G4double irt);
	
	// Get functions
	G4double GetMoleculeRadius(G4int);
	G4int GetMoleculeCharge(G4int);
	G4int GetReactionIndex(G4int pdgA, G4int pdgB);
	G4double GetOnsagerRadius(G4int molA, G4int molB);
	G4double GetRCutOff(G4double tCutOff);
	std::pair<G4String, G4String> GetReactants(G4int);
	std::vector<G4String> GetProducts(G4int);
	std::vector<G4int> GetReactionProducts(G4int index);
	inline G4int GetNumberOfReactions() { return (G4int)fReactions.size();}
	inline G4String GetMoleculeNameFromMoleculeID(G4int molID) {return fMoleculesName[molID];}
	inline std::map<G4int, TsIRTConfiguration::TsMolecularReaction> GetReactions() {return fReactions;}
	inline TsMolecularReaction GetReaction(G4int index) {return fReactions[index];}
	inline G4int GetLastMoleculeID() { return fLastMoleculeID; }
	inline G4int GetLastReactionID() { return fReactionID-1; }

	std::map<G4int, std::vector<std::pair<G4int,G4int>>> GetReactability() {return fMoleculeCanReactWith;};
	void Diffuse(TsMolecule& mol, G4double dt);
	
	void PrintMoleculesInformation();
	void PrintReactionsInformation();
	
	// PH temperature adjust functions
	void AdjustReactionRateForPH(G4String);
	void AdjustReactionAndDiffusionRateForTemperature();
	std::vector<G4double> GetH2SO4ComponentsConcentrationP(G4double);
	std::vector<G4double> GetH2SO4ComponentsConcentrationPH(G4double);
	G4double GetIonicStrength(std::vector<G4double>);
	inline std::map<G4String, G4int> GetMoleculeIDs() {return fMoleculesID;};
	inline std::map<G4int, G4String> GetMoleculeNames() { return fMoleculesName;};
	G4double IonicRate(G4double, TsMolecularReaction);

};
#endif

