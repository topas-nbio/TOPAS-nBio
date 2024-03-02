//
// ********************************************************************
// *                                                                  *
// * This  code  implementation is the  intellectual property  of the *
// * TOPAS collaboration.                                             *
// * Use or redistribution of this code is not permitted without the  *
// * explicit approval of the TOPAS collaboration.                    *
// * Contact: Joseph Perl, perl@slac.stanford.edu                     *
// *                                                                  *
// ********************************************************************
//

#ifndef TsScoreWithIRTCummulativeInTOPAS2_hh
#define TsScoreWithIRTCummulativeInTOPAS2_hh

#include "TsVNtupleScorer.hh"
#include "TsIRTConfiguration.hh"
#include "Randomize.hh"

#include <stdint.h>

class TsHybridIRT;
class TsIRTUtils;

class TsScoreWithIRTCummulativeInTOPAS2 : public TsVNtupleScorer
{
public:
	TsScoreWithIRTCummulativeInTOPAS2(TsParameterManager* pM, TsMaterialManager* mM, TsGeometryManager* gM, TsScoringManager* scM, TsExtensionManager* eM,
							  G4String scorerName, G4String quantity, G4String outFileName, G4bool isSubScorer);
	~TsScoreWithIRTCummulativeInTOPAS2();
	
	virtual G4bool ProcessHits(G4Step*, G4TouchableHistory*);
	
	void UserHookForEndOfEvent();
    
    virtual void UserHookForPreTimeStepAction();
    void AbsorbResultsFromWorkerScorer(TsVScorer*);
    void SampleShiftTime();

    void RunAndSaveInfo();
    void RunChemistry();
		
protected:
	void Output();
	void OutputForThread(G4int);
	
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
	TsHybridIRT* fIRT;
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
	std::vector<G4double> fVEnergyDeposit;
	std::vector<G4double> fAccumulatedEnergy;
	std::vector<G4double> fAccumulatedEnergy2;
	std::vector<G4double> fPrimaryTimes;

	G4int fCurrentPulse;
	G4int fEventID;
	G4int fOldEvent;
	G4double fShiftTime;
	G4double fMinShiftTime;
	
	G4double* fTimeValues;
	G4double* fTimeWeights;
	G4double fVolume;
	std::vector<G4double> fTimeTops;
	G4int fNbOfTimes;
	std::ofstream fTimeOutFile;
	G4bool fLowLimitTo1ps;
	G4bool fSavePulseInformation;
	G4bool fTimeSampled;
	G4bool fReportEachHistory;
	G4bool fReportDelta;

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
};

#endif


