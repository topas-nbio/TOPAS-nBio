//
// ********************************************************************
// *                                                                  *
// * This file is part of the TOPAS-nBio extensions to the            *
// *   TOPAS Simulation Toolkit.                                      *
// * The TOPAS-nBio extensions are freely available under the license *
// *   agreement set forth at: https://topas-nbio.readthedocs.io/     *
// *                                                                  *
// ********************************************************************
//

#ifndef TsScoreDNADamagePulsed_hh
#define TsScoreDNADamagePulsed_hh

#include "TsVNtupleScorer.hh"
#include "TsIRTConfiguration.hh"
#include "Randomize.hh"

#include <stdint.h>

class TsIRTManager;
class TsIRTUtils;

class TsScoreDNADamagePulsed : public TsVNtupleScorer
{
public:
	TsScoreDNADamagePulsed(TsParameterManager* pM, TsMaterialManager* mM, TsGeometryManager* gM, TsScoringManager* scM, TsExtensionManager* eM,
							  G4String scorerName, G4String quantity, G4String outFileName, G4bool isSubScorer);
	~TsScoreDNADamagePulsed();
	
	G4bool ProcessHits(G4Step*, G4TouchableHistory*);
	
	void UserHookForEndOfEvent();
    
    void UserHookForPreTimeStepAction();
    void AbsorbResultsFromWorkerScorer(TsVScorer*);
    void SampleShiftTime();
    void GetDNAInformation();
    void InsertDNAMolecules();
    void RunAndSaveInfo();
    void RunChemistry();

    void FilterMoleculesInsideDNAVolumes();
    void CheckForDirectDNABreaks();
    void CheckForIndirectDNABreaks();

    G4String SampleDNAMolecule();
		
protected:
	void Output();
	
	// Output variables
	G4double fGValue;
	G4double fGValueError;
	G4double fTime;
	G4String fMoleculeName;
	G4double fEnergy;
	G4double fEnergyError;
	G4double fNumberOfMolecules;
	G4double fNumberOfMoleculesError;

	G4double fTotalEnergyDeposit;
	
private:
	G4bool Inside(G4ThreeVector);
	
	TsParameterManager* fPm;
	TsIRTManager* fIRT;
	TsIRTUtils* fUtils;
	
	std::vector<std::pair<G4double,G4double>> fPulseInformation;
	std::vector<std::pair<G4double,G4double>> fVEnergyDepositPerEvent;
	std::vector<std::pair<G4double,G4double>> fVEnergyDepositPerEventEverywhere;

	std::map<G4String, std::map<G4double,G4double>> fGValuePerSpeciePerTime;
	std::map<G4String, std::map<G4double,G4double>> fGValuePerSpeciePerTime2;

	std::map<G4String, std::map<G4double,G4double>> fSpeciePerTime;
	std::map<G4String, std::map<G4double,G4double>> fSpeciePerTime2;

	std::map<G4int, std::map<G4double, G4double>> fDeltaGValues;
	std::map<G4int, std::map<G4double, G4double>> fDeltaGValues2;

	std::vector<G4double> fRandomTimes;
	
	G4double fEnergyDepositPerEvent;
	G4double fEnergyDepositPerEventEverywhere;

	G4double fMass;
	G4double fTotalEDep;
	G4double fPrescribedEDep;
	G4double fTotalDose;
	G4double fPrescribedDose;
	G4int fNbOfScoredEvents;
	G4int fNbOfScoredEventsEverywhere;
	G4int fNbOfIRTRuns;
	G4int fNumberOfRepetitions;
	G4int fNumberOfThreads;
	G4double fRadius;

	G4double fTCut;
	G4String fName;
	G4String fComponent;
	
	G4double fTimeMean;
	G4double fTimeStdv;
	G4double fTimeFWHM;
	G4String fTimeDistribution;
	G4int fTimeDistributionType;
    
    G4int fNumberOfPulses;
    G4double fPulsesTimeDelay;
    G4double fDosePerPulse;
	G4double fEDepPerPulse;
    G4double fPulseTimeShift;
	
	std::vector<G4double> fStepTimes;
	std::vector<G4double> fVEnergyDeposit;
	std::vector<G4double> fAccumulatedEnergy;
	std::vector<G4double> fAccumulatedEnergy2;
	std::vector<G4double> fPrimaryTimes;

	G4int fCurrentPulse;
	G4int fEventID;
	G4int fOldEvent;
	G4double fShiftTime;
	G4double fMinShiftTime;
	G4double fTimeSimulationForPulse;
	
	G4double* fTimeValues;
	G4double* fTimeWeights;
	G4double fVolume;
	std::vector<G4double> fTimeTops;
	G4int fNbOfTimes;
	std::ofstream fTimeOutFile;
	G4bool fLowLimitTo1ps;
	G4bool fSavePulseInformation;
	G4bool fTimeSampled;
	G4bool fReportDelta;
	G4bool fPulseRecycle;

	G4bool fTestIsInside;
    G4String fSensitiveVolume;
    G4String fRootFileName;

    G4bool fSimulationFinished;
    G4String fOutputFile;

    // Molecules
    std::map<G4int,std::pair<G4double,G4double>> fPulseTimeLimits;
    std::map<G4int, std::vector<TsIRTConfiguration::TsMolecule>> fMolecules; 

    std::map<G4String, std::map<G4double, G4int>> fGValuesInRun;
    std::map<G4int, std::map<G4double, G4int>>    fDeltaGInRun;

    // Molecules to insert
    std::vector<G4String> fDNAMoleculesToInsert;
    std::vector<G4double> fDNAMoleculesWeights;

    // DNA Information from Geometry
    std::vector<G4bool>   fPHSPIsAlive;
    std::vector<G4String> fPHSPMolecules;
    std::vector<G4double> fPHSPTime;
    std::vector<G4ThreeVector> fPHSPPosition;
    std::vector<std::vector<G4int>> fDNADetails;

    // DNA Strand Break quantities
    G4bool fDNAHasBeenInserted;

    // Direct DNA damage containers
    std::map<G4int, std::map<G4int, std::map<G4int, G4double>>> fDNAEnergyDeposit;
    std::map<G4int,std::vector<std::vector<G4int>>> fDirectStrandBreaks;
    std::vector<std::vector<G4int>> fDirectStrandBreaksInEvent;

    // Indirect DNA damage containers
    //TsMoleculeGun fMoleculeGun;
    std::vector<G4String> fStrandBreakMolecules;
    std::map<G4int,std::vector<std::vector<G4int>>> fIndirectStrandBreaks;
};

#endif


