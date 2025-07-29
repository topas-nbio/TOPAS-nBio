//
//
// ********************************************************************
// *                                                                  *
// * This file is part of the TOPAS-nBio extensions to the            *
// *   TOPAS Simulation Toolkit.                                      *
// * The TOPAS-nBio extensions are freely available under the license *
// *   agreement set forth at: https://topas-nbio.readthedocs.io/     *
// *                                                                  *
// ********************************************************************

#ifndef TsScoreWithIRTAndGillespieContinuous_hh
#define TsScoreWithIRTAndGillespieContinuous_hh

#include "TsVNtupleScorer.hh"
#include "TsIRTConfiguration.hh"
#include "Randomize.hh"
#include "TsGillespie.hh"

#include <stdint.h>

class TsVIRTProcedure;
class TsIRTContinuous;
class TsIRTUtils;

class TsScoreWithIRTAndGillespieContinuous : public TsVNtupleScorer
{
public:
	TsScoreWithIRTAndGillespieContinuous(TsParameterManager* pM, TsMaterialManager* mM, TsGeometryManager* gM, TsScoringManager* scM, TsExtensionManager* eM,
							  G4String scorerName, G4String quantity, G4String outFileName, G4bool isSubScorer);
	~TsScoreWithIRTAndGillespieContinuous();
	
	virtual G4bool ProcessHits(G4Step*, G4TouchableHistory*);
	
	void UserHookForEndOfEvent();
    
    virtual void UserHookForPreTimeStepAction();
    void AbsorbResultsFromWorkerScorer(TsVScorer*);
    void SampleShiftTime();

    void RunAndSaveInfo();
    void RunChemistry();
		
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
	TsIRTContinuous* fIRT;
	TsIRTUtils* fUtils;
	TsGillespie* fGillespie;

	std::map<G4String, std::map<G4double,G4double>> fGValuePerSpeciePerTime;
	std::map<G4String, std::map<G4double,G4double>> fGValuePerSpeciePerTime2;

	std::map<G4int, std::map<G4double, G4double>> fDeltaGValues;
	std::map<G4int, std::map<G4double, G4double>> fDeltaGValues2;
	

	std::map<G4int, std::map<G4double, G4double>> fDeltaGValues_nh;
	std::map<G4int, std::map<G4double, G4double>> fDeltaGValues2_nh;

	G4double fEnergyDepositPerEvent;

	G4double fTotalDose;
	G4double fPrescribedDose;
	G4int fNbOfScoredEvents;
	G4int fNbOfScoredEventsEverywhere;
	G4int fNbOfIRTRuns;
	G4int fNumberOfRepetitions;
	G4int fNumberOfThreads;

	G4double fTCut;
	G4String fName;
	
	G4double fTimeMean;
	G4double fTimeStdv;
	G4double fTimeFWHM;
	G4String fTimeDistribution;
	G4int fTimeDistributionType;
    
    G4int fNumberOfPulses;
    G4double fPulsesTimeDelay;
    G4double fDosePerPulse;
    G4double fPulseTimeShift;
	
	std::vector<G4double> fStepTimes;

	G4int fCurrentPulse;
	G4int fEventID;
	G4int fOldEvent;
	G4double fShiftTime;
	G4double fMinShiftTime;
	
	G4double* fTimeValues;
	G4double* fTimeWeights;
	G4int fNbOfTimes;
	std::ofstream fTimeOutFile;
	G4bool fLowLimitTo1ps;
	G4bool fSavePulseInformation;
	G4bool fReportDelta;
	G4bool fPulseRecycle;

	G4bool fTestIsInside;
    G4String fSensitiveVolume;
    G4String fRootFileName;

    G4bool fSimulationFinished;
    G4String fOutputFile;

    // Molecules
    std::map<G4int, std::vector<TsIRTConfiguration::TsMolecule>> fMolecules;

    std::map<G4String, std::map<G4double, G4int>> fGValuesInRun;
    std::map<G4int, std::map<G4double, G4int>>    fDeltaGInRun;

    std::map<G4int, std::map<G4double, G4int>>    fDeltaGInRun_nh;

    std::vector<std::pair<G4double, G4double>> fEnergyDep;

    G4double fIrradiationTime;
    G4double fGillespieWorld;
    G4double fSideLength;
};

#endif


