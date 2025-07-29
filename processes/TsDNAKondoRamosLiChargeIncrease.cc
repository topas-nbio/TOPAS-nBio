// Extra Class for TsEmDNAPhysics
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

#include "TsDNAKondoRamosLiChargeIncrease.hh"
#include "G4PhysicalConstants.hh"
#include "G4SystemOfUnits.hh"
#include "G4DNAMolecularMaterial.hh"
#include "G4DNAGenericIonsManager.hh"

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

using namespace std;

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

TsDNAKondoRamosLiChargeIncrease::TsDNAKondoRamosLiChargeIncrease(const G4ParticleDefinition*,
                                                                       const G4String& nam) :
G4VEmModel(nam), isInitialised(false)
{
  fpMolWaterDensity = 0;

  verboseLevel = 1;
  // Verbosity scale:
  // 0 = nothing
  // 1 = warning for energy non-conservation
  // 2 = details of energy budget
  // 3 = calculation of cross sections, file openings, sampling of atoms
  // 4 = entering in methods

  if (verboseLevel > 0)
  {
    G4cout << "Dingfelder charge increase model is constructed " << G4endl;
  }
  fParticleChangeForGamma = 0;

  // Selection of stationary mode

  statCode = false;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

TsDNAKondoRamosLiChargeIncrease::~TsDNAKondoRamosLiChargeIncrease()
{}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

void TsDNAKondoRamosLiChargeIncrease::Initialise(const G4ParticleDefinition* particle,
                                                    const G4DataVector& /*cuts*/)
{
  if (verboseLevel > 3)
  {
    G4cout << "Calling TsDNAKondoRamosLiChargeIncrease::Initialise()"
        << G4endl;
  }

  // Energy limits
  G4DNAGenericIonsManager *instance;
  instance = G4DNAGenericIonsManager::Instance();
  lithium3pDef = instance->GetIon("lithium+++");
  lithium2pDef = instance->GetIon("lithium++");
  lithium1pDef = instance->GetIon("lithium+");
  lithium0pDef = instance->GetIon("lithium");

  G4String lithium2p;
  G4String lithium1p;
  G4String lithium0p;

  // Limits
  lithium2p = lithium2pDef->GetParticleName();
  lowEnergyLimit[lithium2p]  = 7 * keV;
  highEnergyLimit[lithium2p] = 70 * MeV;

  lithium1p = lithium1pDef->GetParticleName();
  lowEnergyLimit[lithium1p]  = 7 * keV;
  highEnergyLimit[lithium1p] = 70 * MeV;

  lithium0p = lithium0pDef->GetParticleName();
  lowEnergyLimit[lithium0p]  = 7 * keV;
  highEnergyLimit[lithium0p] = 70 * MeV;

  if (particle==lithium2pDef) {
    SetLowEnergyLimit(lowEnergyLimit[lithium2p]);
    SetHighEnergyLimit(highEnergyLimit[lithium2p]);   
  }

  if (particle==lithium1pDef) {
    SetLowEnergyLimit(lowEnergyLimit[lithium1p]);
    SetHighEnergyLimit(highEnergyLimit[lithium1p]);   
  }

  if (particle==lithium0pDef) {
    SetLowEnergyLimit(lowEnergyLimit[lithium0p]);
    SetHighEnergyLimit(highEnergyLimit[lithium0p]);   
  }

  // Final state
  //                      |   Energy      | Li0+->Li1+ | Li1+->Li2+ | Li2+->Li3+ |
  CrossSections.push_back({7.000000E+03*eV,4.837164E-22,4.394210E-22,2.963844E-22});
  CrossSections.push_back({1.947792E+04*eV,3.900489E-22,2.186432E-22,1.226901E-22});
  CrossSections.push_back({5.419846E+04*eV,7.207329E-22,2.361835E-22,1.226873E-22});
  CrossSections.push_back({1.508104E+05*eV,3.651596E-21,1.477138E-22,1.077393E-22});
  CrossSections.push_back({4.196390E+05*eV,9.008062E-21,2.839547E-22,2.915472E-22});
  CrossSections.push_back({1.167670E+06*eV,1.006124E-20,1.405591E-21,1.860977E-21});
  CrossSections.push_back({3.249112E+06*eV,7.392020E-21,2.569215E-21,3.688974E-21});
  CrossSections.push_back({9.040848E+06*eV,4.144544E-21,3.090820E-21,2.579072E-21});
  CrossSections.push_back({2.515670E+07*eV,2.294121E-21,1.885519E-21,1.315700E-21});
  CrossSections.push_back({7.000000E+07*eV,1.169267E-21,1.146429E-21,6.486156E-22});

  IonIonizationEnergy[lithium2pDef] = 122.45*eV;
  IonIonizationEnergy[lithium1pDef] = 75.637*eV;
  IonIonizationEnergy[lithium0pDef] =   5.39*eV;

  if( verboseLevel>0 )
  {
    G4cout << "Kondo Ramos charge increase model is initialized " << G4endl
    << "Energy range: "
    << LowEnergyLimit() / keV << " keV - "
    << HighEnergyLimit() / MeV << " MeV for "
    << particle->GetParticleName()
    << G4endl;
  }

  // Initialize water density pointer
  fpMolWaterDensity = G4DNAMolecularMaterial::Instance()->GetNumMolPerVolTableFor(G4Material::GetMaterial("G4_WATER"));

  if (isInitialised)
  { return;}
  fParticleChangeForGamma = GetParticleChangeForGamma();
  isInitialised = true;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

G4double TsDNAKondoRamosLiChargeIncrease::CrossSectionPerVolume(const G4Material* material,
                                                                   const G4ParticleDefinition* particleDefinition,
                                                                   G4double k,
                                                                   G4double,
                                                                   G4double)
{
  if (verboseLevel > 3){
    G4cout
        << "Calling CrossSectionPerVolume() of TsDNAKondoRamosChargeIncreaseModel"
        << G4endl;
  }

  // Check if this is a supported Particle
  if (particleDefinition != lithium2pDef &&
      particleDefinition != lithium1pDef &&
      particleDefinition != lithium0pDef)
    return 0;

  // Calculate total cross section for model
  G4double crossSection = 0.;
  G4double waterDensity = (*fpMolWaterDensity)[material->GetIndex()];
  const G4String& particleName = particleDefinition->GetParticleName();

  crossSection = Sum(k,particleDefinition);

  if (verboseLevel > 2)
  {
    G4cout << "_______________________________________" << G4endl;
    G4cout << "TsDNAKondoRamosLiChargeIncrease" << G4endl;
    G4cout << "Kinetic energy(eV)=" << k/eV << "particle :" << particleName << G4endl;
    G4cout << "Cross section per water molecule (cm^2)=" << crossSection/cm/cm << G4endl;
    G4cout << "Cross section per water molecule (cm^-1)=" << crossSection*
    waterDensity/(1./cm) << G4endl;
  }
  return crossSection*waterDensity;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

void TsDNAKondoRamosLiChargeIncrease::SampleSecondaries(std::vector<G4DynamicParticle*>* fvect,
                                                           const G4MaterialCutsCouple*,
                                                           const G4DynamicParticle* aDynamicParticle,
                                                           G4double,
                                                           G4double) {
  if (verboseLevel > 3) {
    G4cout
        << "Calling SampleSecondaries() of TsDNAKondoRamosLiChargeIncrease"
        << G4endl;
  }
  
  if (!statCode) fParticleChangeForGamma->ProposeLocalEnergyDeposit(0.);
  
  G4ParticleDefinition* definition = aDynamicParticle->GetDefinition();
  G4double particleMass = definition->GetPDGMass();

  G4double inK = aDynamicParticle->GetKineticEnergy();


  G4int finalStateIndex = 1;
  G4int n = 1;

  G4double outK = 0.;
  
  if (!statCode) outK = inK - IncomingParticleBindingEnergyConstant(definition,finalStateIndex);

  else outK = inK;
  
  if (statCode) 
    fParticleChangeForGamma-> ProposeLocalEnergyDeposit(IncomingParticleBindingEnergyConstant(definition,finalStateIndex));

  fParticleChangeForGamma->ProposeTrackStatus(fStopAndKill);

  G4double electronK = inK*electron_mass_c2/(particleMass);

  if (outK<0) {
    G4Exception("TsDNAKondoRamosLiChargeIncrease::SampleSecondaries","em0004",
        FatalException,"Final kinetic energy is negative.");
  }

  G4DynamicParticle* dp = new G4DynamicParticle(OutgoingParticleDefinition(definition,finalStateIndex),
      aDynamicParticle->GetMomentumDirection(),outK);

  fvect->push_back(dp);

  n = n - 1;

  while (n>0) {
    n--;
    fvect->push_back(new G4DynamicParticle
        (G4Electron::Electron(), aDynamicParticle->GetMomentumDirection(), electronK) );
  }
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

G4int TsDNAKondoRamosLiChargeIncrease::NumberOfFinalStates(G4ParticleDefinition* particleDefinition,
                                                              G4int finalStateIndex)

{
  if (particleDefinition == lithium2pDef ||
      particleDefinition == lithium1pDef ||
      particleDefinition == lithium0pDef )
    return 1;
  return 0;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

G4ParticleDefinition* TsDNAKondoRamosLiChargeIncrease::OutgoingParticleDefinition(G4ParticleDefinition* particleDefinition,
                                                                                     G4int finalStateIndex) {
  if (particleDefinition == lithium0pDef){
    return lithium1pDef;
  }

  else if (particleDefinition == lithium1pDef){
    return lithium2pDef;
  }

  else if (particleDefinition == lithium2pDef){
    return lithium3pDef;
  }

  return 0;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

G4double TsDNAKondoRamosLiChargeIncrease::IncomingParticleBindingEnergyConstant(G4ParticleDefinition* particleDefinition,
                                                                                   G4int)
{
  // Li3+ -> None :   0.00 eV
  // Li2+ -> Li3+ : 122.45 eV
  // Li1+ -> Li2+ : 75.637 eV
  // Li0+ -> Li1+ :   5.39 eV
  if (particleDefinition == lithium0pDef ||
      particleDefinition == lithium1pDef ||
      particleDefinition == lithium2pDef)
    return IonIonizationEnergy[particleDefinition];

  return 0;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

G4double TsDNAKondoRamosLiChargeIncrease::PartialCrossSection(G4double k,
                                                                 G4int,
                                                                 const G4ParticleDefinition* particleDefinition)
{
  size_t index = 0;
  if (particleDefinition == lithium0pDef)
    index = 1;
  else if (particleDefinition == lithium1pDef)
    index = 2;
  else if (particleDefinition == lithium2pDef)
    index = 3;

  G4double crossSection = CrossSections[0][index];

  // Linear interpolation for energies in between
  for (size_t i = 0; i < CrossSections.size() - 1; i++) {
    G4double Energy1 = CrossSections[i+0][0];
    G4double Energy2 = CrossSections[i+1][0];
    G4double Cross1  = CrossSections[i+0][index];
    G4double Cross2  = CrossSections[i+1][index];

    if ((Energy1 <= k) and (Energy2 >= k)) {
      G4double slope = (Cross2 - Cross1)/(Energy2 - Energy1);
      G4double inter = Cross1 - (slope*Energy1);
      crossSection = slope*k + inter;
      break;
    }
  }
  return crossSection * m * m;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

G4int TsDNAKondoRamosLiChargeIncrease::RandomSelect(G4double,
                                                       const G4ParticleDefinition*)
{
  return 1;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

G4double TsDNAKondoRamosLiChargeIncrease::Sum(G4double k,
                                                 const G4ParticleDefinition* particleDefinition)
{
  G4double totalCrossSection = 0.;
  totalCrossSection += PartialCrossSection(k, 0, particleDefinition);

  return totalCrossSection;
}
