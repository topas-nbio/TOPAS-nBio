//
// ********************************************************************
// * License and Disclaimer                                           *
// *                                                                  *
// * The  Geant4 software  is  copyright of the Copyright Holders  of *
// * the Geant4 Collaboration.  It is provided  under  the terms  and *
// * conditions of the Geant4 Software License,  included in the file *
// * LICENSE and available at  http://cern.ch/geant4/license .  These *
// * include a list of copyright holders.                             *
// *                                                                  *
// * Neither the authors of this software system, nor their employing *
// * institutes,nor the agencies providing financial support for this *
// * work  make  any representation or  warranty, express or implied, *
// * regarding  this  software system or assume any liability for its *
// * use.  Please see the license in the file  LICENSE  and URL above *
// * for the full disclaimer and the limitation of liability.         *
// *                                                                  *
// * This  code  implementation is the result of  the  scientific and *
// * technical work of the GEANT4 collaboration.                      *
// * By using,  copying,  modifying or  distributing the software (or *
// * any work based  on the software)  you  agree  to acknowledge its *
// * use  in  resulting  scientific  publications,  and indicate your *
// * acceptance of all terms of the Geant4 Software license.          *
// ********************************************************************
//
//

#ifndef TsDNAKRIonMillerGreenExcitationModel_h
#define TsDNAKRIonMillerGreenExcitationModel_h 1

#include "G4VEmModel.hh"
#include "G4ParticleChangeForGamma.hh"
#include "G4ProductionCutsTable.hh"

#include "G4Proton.hh"
#include "G4DNAGenericIonsManager.hh"
#include "G4DNAWaterExcitationStructure.hh"
#include "Randomize.hh"
#include "G4NistManager.hh"

#include <deque>

class TsDNAKRIonMillerGreenExcitationModel : public G4VEmModel
{

public:

  TsDNAKRIonMillerGreenExcitationModel(const G4ParticleDefinition* p = 0, const G4String& nam = "DNAKRIonMillerGreenExcitationModel");
  virtual ~TsDNAKRIonMillerGreenExcitationModel();
  virtual void Initialise(const G4ParticleDefinition*, const G4DataVector&);
  virtual G4double CrossSectionPerVolume(const G4Material*,const G4ParticleDefinition*,G4double,G4double,G4double);
  virtual G4double GetPartialCrossSection(const G4Material*,G4int,const G4ParticleDefinition*,G4double);
  virtual void SampleSecondaries(std::vector<G4DynamicParticle*>*,const G4MaterialCutsCouple*,const G4DynamicParticle*,G4double,G4double);
  inline void SelectStationary(G4bool input); 

protected:
  G4ParticleChangeForGamma* fParticleChangeForGamma;

private:
  G4double IonEffectiveCharge(const G4ParticleDefinition*, G4double, G4double);
  G4double IonSlaterEffectiveCharge(G4int, G4int);
  
  G4bool statCode;

  // Water density table
  const std::vector<G4double>* fpMolWaterDensity;

  std::map<G4String,G4double,std::less<G4String> > lowEnergyLimit;
  std::map<G4String,G4double,std::less<G4String> > highEnergyLimit;

  G4bool isInitialised;
  G4int verboseLevel;
  
  G4double fAtomicMassNumber;
  G4double fAtomicNumber;
  G4double fQeff;

  G4double fLowEnergyLimit;
  G4double fHighEnergyLimit;

  std::vector<std::vector<G4double>> EFFECTIVECHART;

  // Cross section

  // Partial cross section
  G4double PartialCrossSection(G4double energy,G4int level, const G4ParticleDefinition* particle);
  G4double Sum(G4double energy, const G4ParticleDefinition* particle);
  G4int RandomSelect(G4double energy, const G4ParticleDefinition* particle);

  G4int nLevels;

  G4DNAWaterExcitationStructure waterExcitation;
   
  TsDNAKRIonMillerGreenExcitationModel & operator=(const  TsDNAKRIonMillerGreenExcitationModel &right);
  TsDNAKRIonMillerGreenExcitationModel(const  TsDNAKRIonMillerGreenExcitationModel&);

};

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

inline void TsDNAKRIonMillerGreenExcitationModel::SelectStationary (G4bool input)
{ 
    statCode = input; 
}		 

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

#endif
