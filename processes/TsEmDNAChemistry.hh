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

#ifndef TsEmDNAChemistry_hh
#define TsEmDNAChemistry_hh 1

#include "G4VPhysicsConstructor.hh"
#include "G4VUserChemistryList.hh"
#include "globals.hh"
#include <map>

class TsParameterManager;

class G4DNAMolecularReactionTable;

class TsEmDNAChemistry
: public G4VUserChemistryList,
  public G4VPhysicsConstructor
{
    
public:
    
    explicit TsEmDNAChemistry();
    TsEmDNAChemistry(TsParameterManager* pM);
    
    virtual ~TsEmDNAChemistry();
    
    virtual void ConstructParticle()
    {
        ConstructMolecule();
    }
    
    virtual void ConstructMolecule();
    virtual void ConstructProcess();
    
    virtual void ConstructDissociationChannels();
    virtual void ConstructReactionTable(G4DNAMolecularReactionTable* reactionTable);
    virtual void ConstructTimeStepModel(G4DNAMolecularReactionTable* reactionTable);
    
    G4bool IsWaterMolecularConfigurationActivated();
    
private:
    void DefineParameters();
    G4bool MoleculeExists(G4String);
    G4String GetFullParmName(G4String);
    void Quit(G4String, G4String);
    
private:
    TsParameterManager* fPm;
    std::map<G4String, G4String> fExistingMolecules;
    G4String fName;
    
    std::map<G4String, G4double> fDiffusionCoefficients;
    std::vector< std::vector<G4String> > fReactionSpecies;
    std::vector< std::vector<G4String> > fReactionProducts;
    std::vector<G4double> fReactionRates;
    
    G4bool fSetWaterConfiguration;

    // Ionisation Branching ratio
    G4double fIonizationStates;

    // Fifth Layer Excitation
    G4double fA1B1Relaxation;
    G4double fA1B1DissociativeDecay;

    // Fourth Layer Excitation
    G4double fB1A1Relaxation;
    G4double fB1A1DissociativeDecay;
    G4double fB1A1AutoIonization;
    G4double fA1B1DissociativeDecay2;
    G4double fB1A1DissociativeDecay2;

    // Third, Second and First Layer Excitation
    G4double fRydDiffAutoIonization;
    G4double fRydDiffRelaxation;

    // Dissociative Attachment
    G4double fDissociativeAttachment;

    // Electron Hole Recombination
    G4double fH2OvibDissociationDecay1;
    G4double fH2OvibDissociationDecay2;
    G4double fH2OvibDissociationDecay3;
    G4double fH2OvibDissociationDecay4;

    std::vector<G4String> fScavengedMolecules;
    std::vector<G4double> fScavengedCapacities;
    std::vector<G4bool> fScavengerHasProducts;
    std::vector<std::vector<G4String>> fScavengerProducts;
	std::vector<G4String> fRemoveInMaterialTheseMolecules;
	G4String fRemoveInMaterial;
};

#endif
