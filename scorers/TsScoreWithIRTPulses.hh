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

#ifndef TsScoreWithIRTPulses_hh
#define TsScoreWithIRTPulses_hh

#include "TsVNtupleScorer.hh"
#include "TsIRTConfiguration.hh"

#include "G4Timer.hh"
#include "Randomize.hh"

#include <stdint.h>

class TsIRTManager;

class TsScoreWithIRTPulses : public TsVNtupleScorer
{
public:
    TsScoreWithIRTPulses(TsParameterManager* pM, TsMaterialManager* mM, TsGeometryManager* gM, TsScoringManager* scM, TsExtensionManager* eM,
                      G4String scorerName, G4String quantity, G4String outFileName, G4bool isSubScorer);
    ~TsScoreWithIRTPulses();
    
    virtual G4bool ProcessHits(G4Step*, G4TouchableHistory*);
    void AbsorbResultsFromWorkerScorer(TsVScorer* workerScorer);
 
    void UserHookForEndOfEvent();
	
	virtual void UserHookForPreTimeStepAction();
	
protected:
    
    void Output();
    void Clear();
    void SampleTimeShift();
 
    // Output variables
    G4double fGValue;
    G4double fTime;
    G4String fMoleculeName;
    G4double fMolecules;
    G4double fEnergyDepositGvalue;
	
    std::map<G4String, std::map<G4double, G4double> > fGValuePerSpeciePerTime;
    std::map<G4String, std::map<G4double, G4double> > fMoleculesPerSpeciePerTime;
    std::map<G4int, std::map<G4double, G4double> > fDeltaGPerReactionPerTime;
   
private:
    TsParameterManager* fPm;
	TsIRTManager* fIRT;
	TsIRTUtils* fUtils;

	std::vector<TsIRTConfiguration::TsMolecule> fSpecies;

    G4double fEnergyDepositPerEvent;

    // Pulse definition
    G4double fTotalEnergyDeposit;
    G4double fPrescribedDose;
    G4double fDosePerPulse;
    G4double fDensity;
    G4double fMass;
    
    G4double fShiftTime;
    G4int fEventID;
    G4int fOldEvent;
    G4double fTimeMean;
    G4double fTimeFWHM;
    G4double fTimeStdv;
    G4int fNumberOfPulses;
    G4int fPulseCount;
    G4double fPulsesTimePeriod; // = 1/frequency
    G4double fPulseTimeShift; // Added to the molecules
    	
	std::vector<G4double> fVEnergyDeposits;
	std::vector<G4double> fVStepTimes;
    std::vector<std::pair<G4double,G4double>> fVEnergyDepositPerSampledTime;

    std::ofstream fTimeOutFile;
};

#endif

