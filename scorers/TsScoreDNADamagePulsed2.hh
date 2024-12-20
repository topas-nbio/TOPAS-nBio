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

#ifndef TsScoreDNADamagePulsed2_hh
#define TsScoreDNADamagePulsed2_hh

#include "TsVNtupleScorer.hh"
#include "TsIRTConfiguration.hh"
#include "Randomize.hh"

#include <stdint.h>

class TsIRTManager;
class TsIRTUtils;

class TsScoreDNADamagePulsed2 : public TsVNtupleScorer
{
public:
	TsScoreDNADamagePulsed2(TsParameterManager* pM, TsMaterialManager* mM, TsGeometryManager* gM, TsScoringManager* scM, TsExtensionManager* eM,
							  G4String scorerName, G4String quantity, G4String outFileName, G4bool isSubScorer);
	~TsScoreDNADamagePulsed2();
	
	virtual G4bool ProcessHits(G4Step*, G4TouchableHistory*);
	
	void UserHookForEndOfEvent();
    
    virtual void UserHookForPreTimeStepAction();
    
    void InsertDNAMolecules();
    G4String SampleDNAMolecule();
    void GetDNAInformation();

		
protected:
	
	// Output variables
	G4double fGValue;
	G4double fTime;
	G4String fMoleculeName;
	G4double fEnergy;
	
private:
	G4bool Inside(G4ThreeVector);
	
	TsParameterManager* fPm;
	TsIRTManager* fIRT;
	TsIRTUtils* fUtils;
	
	std::vector<std::pair<G4double,G4double>> fVEnergyDepositPerEvent;
	std::vector<std::pair<G4double,G4double>> fVEnergyDepositPerEventEverywhere;
	std::vector<G4double> fRandomTimes;
	
	G4double fEnergyDepositPerEvent;
	G4double fEnergyDepositPerEventEverywhere;

    G4double fMass;
    G4double fDensity;
	G4double fTotalDose;
	G4double fPrescribedDose;
	G4int fNbOfScoredEvents;
	G4int fNbOfScoredEventsEverywhere;

	G4double fTCut;
	G4String fName;
    G4String fOutputFile;
    G4String fComponentName;
	
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
	std::vector<G4double> fVEnergyDepositEverywhere;

	G4int fEventID;
	G4int fOldEvent;
	G4double fShiftTime;
	G4double fMinShiftTime;
	
	G4double* fTimeValues;
	G4double* fTimeWeights;
	std::vector<G4double> fTimeTops;
	G4int fNbOfTimes;
	std::ofstream fTimeOutFile;
	G4bool fLowLimitTo1ps;
	
	G4bool fTestIsInside;
    G4String fSensitiveVolume;
    
    // DNA Information from Geometry
    std::vector<G4bool>   fPHSPIsAlive;
    std::vector<G4String> fPHSPMolecules;
    std::vector<G4double> fPHSPTime;
    std::vector<G4ThreeVector> fPHSPPosition;
    std::vector<std::vector<G4int>> fDNADetails;
    // Molecules to insert
    std::vector<G4String> fDNAMoleculesToInsert;
    std::vector<G4double> fDNAMoleculesWeights;
    
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


