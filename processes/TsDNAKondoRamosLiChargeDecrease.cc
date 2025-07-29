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

#include "TsDNAKondoRamosLiChargeDecrease.hh"
#include "G4PhysicalConstants.hh"
#include "G4SystemOfUnits.hh"
#include "G4DNAMolecularMaterial.hh"
#include "G4DNAChemistryManager.hh"
#include "G4DNAGenericIonsManager.hh"

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

using namespace std;

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

TsDNAKondoRamosLiChargeDecrease::TsDNAKondoRamosLiChargeDecrease(const G4ParticleDefinition*,const G4String& nam):
G4VEmModel(nam), isInitialised(false)
{
  verboseLevel = 0;
  // Verbosity scale:
  // 0 = nothing
  // 1 = warning for energy non-conservation
  // 2 = details of energy budget
  // 3 = calculation of cross sections, file openings, sampling of atoms
  // 4 = entering in methods

  if (verboseLevel > 0)
  {
    G4cout << "Kondo-Ramos charge decrease model is constructed " << G4endl;
  }
  // Selection of stationary mode

  statCode = false;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

void TsDNAKondoRamosLiChargeDecrease::Initialise(const G4ParticleDefinition* particle,const G4DataVector& /*cuts*/)
{

  if (verboseLevel > 3)
  {
    G4cout << "Calling TsDNAKondoRamosLiChargeDecrease::Initialise()"
        << G4endl;
  }

  // Energy limits
  G4DNAGenericIonsManager *instance;
  instance = G4DNAGenericIonsManager::Instance();
  lithium3pDef = instance->GetIon("lithium+++");
  lithium2pDef = instance->GetIon("lithium++");
  lithium1pDef = instance->GetIon("lithium+");
  lithium0pDef = instance->GetIon("lithium");

  G4String lithium3p;
  G4String lithium2p;
  G4String lithium1p;

  // Limits
  lithium3p = lithium3pDef->GetParticleName();
  lowEnergyLimit[lithium3p]  = 7 * keV;
  highEnergyLimit[lithium3p] = 70 * MeV;

  lithium2p = lithium2pDef->GetParticleName();
  lowEnergyLimit[lithium2p]  = 7 * keV;
  highEnergyLimit[lithium2p] = 70 * MeV;

  lithium1p = lithium1pDef->GetParticleName();
  lowEnergyLimit[lithium1p]  = 7 * keV;
  highEnergyLimit[lithium1p] = 70 * MeV;


  if (particle==lithium3pDef) {
    SetLowEnergyLimit(lowEnergyLimit[lithium3p]);
    SetHighEnergyLimit(highEnergyLimit[lithium3p]);   
  }

  if (particle==lithium2pDef) {
    SetLowEnergyLimit(lowEnergyLimit[lithium2p]);
    SetHighEnergyLimit(highEnergyLimit[lithium2p]);   
  }

  if (particle==lithium1pDef) {
    SetLowEnergyLimit(lowEnergyLimit[lithium1p]);
    SetHighEnergyLimit(highEnergyLimit[lithium1p]);   
  }

  // Final state
  // Li1+ -> Li0+
  CrossSections[lithium1pDef].push_back({7.000000E+03*eV,1.388372E-20,5.281229E-20,7.772002E-20,7.455183E-20,7.299424E-20});
  CrossSections[lithium1pDef].push_back({1.947792E+04*eV,7.455694E-21,3.042259E-20,4.997593E-20,5.304524E-20,5.442689E-20});
  CrossSections[lithium1pDef].push_back({5.419846E+04*eV,1.759910E-21,1.276839E-20,2.952734E-20,3.571974E-20,4.136280E-20});
  CrossSections[lithium1pDef].push_back({1.508104E+05*eV,2.333072E-22,3.555758E-21,1.542107E-20,1.872622E-20,2.230467E-20});
  CrossSections[lithium1pDef].push_back({4.196390E+05*eV,1.924410E-24,1.356233E-21,4.971195E-21,5.492790E-21,6.351585E-21});
  CrossSections[lithium1pDef].push_back({1.167670E+06*eV,4.834599E-27,3.785287E-22,5.987575E-22,6.857671E-22,7.872640E-22});
  CrossSections[lithium1pDef].push_back({3.249112E+06*eV,2.419467E-26,2.080827E-23,3.474311E-23,3.881663E-23,5.434617E-23});
  CrossSections[lithium1pDef].push_back({9.040848E+06*eV,4.411336E-26,7.626291E-25,6.389405E-25,3.940552E-25,3.079267E-25});
  CrossSections[lithium1pDef].push_back({2.515670E+07*eV,0.000000E+00,0.000000E+00,0.000000E+00,0.000000E+00,0.000000E+00});
  CrossSections[lithium1pDef].push_back({7.000000E+07*eV,0.000000E+00,0.000000E+00,0.000000E+00,0.000000E+00,0.000000E+00});

  // Li2+ -> Li1+
  CrossSections[lithium2pDef].push_back({7.000000E+03*eV,2.183450E-20,7.571585E-20,1.097461E-19,9.557622E-20,9.822393E-20});
  CrossSections[lithium2pDef].push_back({1.947792E+04*eV,1.440076E-20,5.209166E-20,8.319713E-20,7.660796E-20,7.731910E-20});
  CrossSections[lithium2pDef].push_back({5.419846E+04*eV,5.612017E-21,2.170462E-20,4.477613E-20,4.700860E-20,5.194373E-20});
  CrossSections[lithium2pDef].push_back({1.508104E+05*eV,4.147418E-22,4.776992E-21,1.654478E-20,2.139252E-20,2.644916E-20});
  CrossSections[lithium2pDef].push_back({4.196390E+05*eV,2.389398E-23,1.240803E-21,5.100385E-21,5.821666E-21,6.828634E-21});
  CrossSections[lithium2pDef].push_back({1.167670E+06*eV,6.800994E-26,3.673736E-22,5.634587E-22,6.213708E-22,6.897739E-22});
  CrossSections[lithium2pDef].push_back({3.249112E+06*eV,7.515617E-26,2.504243E-23,3.860696E-23,4.099779E-23,4.415964E-23});
  CrossSections[lithium2pDef].push_back({9.040848E+06*eV,6.670820E-26,1.578399E-24,7.253453E-25,5.778991E-25,5.964252E-25});
  CrossSections[lithium2pDef].push_back({2.515670E+07*eV,0.000000E+00,0.000000E+00,0.000000E+00,0.000000E+00,0.000000E+00});
  CrossSections[lithium2pDef].push_back({7.000000E+07*eV,0.000000E+00,0.000000E+00,0.000000E+00,0.000000E+00,0.000000E+00});

  // Li3+ -> Li2+
  CrossSections[lithium3pDef].push_back({7.000000E+03*eV,3.078846E-20,1.127884E-19,1.597307E-19,1.436304E-19,1.441565E-19});
  CrossSections[lithium3pDef].push_back({1.947792E+04*eV,2.303027E-20,8.657864E-20,1.287623E-19,1.215686E-19,1.221344E-19});
  CrossSections[lithium3pDef].push_back({5.419846E+04*eV,9.314896E-21,4.483006E-20,8.068468E-20,8.030346E-20,8.921234E-20});
  CrossSections[lithium3pDef].push_back({1.508104E+05*eV,1.385545E-21,9.178010E-21,2.864419E-20,3.625581E-20,4.513544E-20});
  CrossSections[lithium3pDef].push_back({4.196390E+05*eV,3.077101E-23,1.958089E-21,8.385574E-21,1.013619E-20,1.260920E-20});
  CrossSections[lithium3pDef].push_back({1.167670E+06*eV,4.373020E-24,4.628501E-22,8.705857E-22,1.083426E-21,1.150265E-21});
  CrossSections[lithium3pDef].push_back({3.249112E+06*eV,2.358760E-26,2.046628E-23,3.564790E-23,4.711811E-23,5.519495E-23});
  CrossSections[lithium3pDef].push_back({9.040848E+06*eV,7.894388E-26,6.535859E-25,3.967109E-25,5.149671E-25,2.505713E-25});
  CrossSections[lithium3pDef].push_back({2.515670E+07*eV,0.000000E+00,0.000000E+00,0.000000E+00,0.000000E+00,0.000000E+00});
  CrossSections[lithium3pDef].push_back({7.000000E+07*eV,0.000000E+00,0.000000E+00,0.000000E+00,0.000000E+00,0.000000E+00});

  WaterBindingEnergy[1] = 12.62*eV;
  WaterBindingEnergy[2] = 14.75*eV;
  WaterBindingEnergy[3] = 18.51*eV;
  WaterBindingEnergy[4] =  32.4*eV;
  WaterBindingEnergy[5] = 539.7*eV;

  IonIonizationEnergy[lithium3pDef] =   0.00*eV;
  IonIonizationEnergy[lithium2pDef] = 122.45*eV;
  IonIonizationEnergy[lithium1pDef] = 75.637*eV;
  IonIonizationEnergy[lithium0pDef] =   5.39*eV;

  if( verboseLevel>0 )
  {
    G4cout << "Kondo Ramos charge decrease model is initialized " << G4endl
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

G4double TsDNAKondoRamosLiChargeDecrease::CrossSectionPerVolume(const G4Material* material,
                                                                   const G4ParticleDefinition* particleDefinition,
                                                                   G4double k,
                                                                   G4double,
                                                                   G4double)
{
  if (verboseLevel > 3)
  {
    G4cout
        << "Calling CrossSectionPerVolume() of TsDNAKondoRamosLiChargeDecrease"
        << G4endl;
  }

  // Check if this is a supported Particle
  if (particleDefinition != lithium3pDef &&
      particleDefinition != lithium2pDef &&
      particleDefinition != lithium1pDef)
    return 0;

  // Calculate total cross section for model
  G4double crossSection = 0.;
  G4double waterDensity = (*fpMolWaterDensity)[material->GetIndex()];
  const G4String& particleName = particleDefinition->GetParticleName();

  crossSection = Sum(k,particleDefinition);

  if (verboseLevel > 2)
  {
    G4cout << "_______________________________________" << G4endl;
    G4cout << "TsDNAKondoRamosLiChargeDecrease" << G4endl;
    G4cout << "Kinetic energy(eV)=" << k/eV << "particle :" << particleName << G4endl;
    G4cout << "Cross section per water molecule (cm^2)=" << crossSection/cm/cm << G4endl;
    G4cout << "Cross section per water molecule (cm^-1)=" << crossSection*
    waterDensity/(1./cm) << G4endl;
  }

  return crossSection*waterDensity;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

void TsDNAKondoRamosLiChargeDecrease::SampleSecondaries(std::vector<G4DynamicParticle*>* fvect,
                                                           const G4MaterialCutsCouple* /*couple*/,
                                                           const G4DynamicParticle* aDynamicParticle,
                                                           G4double,
                                                           G4double)
{
  if (verboseLevel > 3)
  {
    G4cout
        << "Calling SampleSecondaries() of TsDNAKondoRamosLiChargeDecrease"
        << G4endl;
  }

  G4double inK = aDynamicParticle->GetKineticEnergy();
  G4ParticleDefinition* definition = aDynamicParticle->GetDefinition();

  G4double particleMass = definition->GetPDGMass();
  G4int finalStateIndex = RandomSelect(inK,definition);

  G4int n = 1; // No multiple electron capture
  G4double waterBindingEnergy = WaterBindingEnergyConstant(definition, finalStateIndex);
  G4double outgoingParticleBindingEnergy = OutgoingParticleBindingEnergyConstant(definition, finalStateIndex);

  G4double outK = 0.;
  
  if (!statCode) outK = inK - n*(inK*electron_mass_c2/particleMass) 
                             - waterBindingEnergy + outgoingParticleBindingEnergy;
  else outK = inK;
  
  if (outK<0) {
    G4Exception("TsDNAKondoRamosLiChargeDecrease::SampleSecondaries","em0004",
        FatalException,"Final kinetic energy is negative.");
  }

  fParticleChangeForGamma->ProposeTrackStatus(fStopAndKill);

  if (!statCode) {
    fParticleChangeForGamma->ProposeLocalEnergyDeposit(waterBindingEnergy);
  }
  
  else {
     fParticleChangeForGamma->ProposeLocalEnergyDeposit(n*(inK*electron_mass_c2/particleMass) 
     + waterBindingEnergy - outgoingParticleBindingEnergy);
  }

  G4DynamicParticle* dp = new G4DynamicParticle (OutgoingParticleDefinition(definition, finalStateIndex),
      aDynamicParticle->GetMomentumDirection(),
      outK);
  fvect->push_back(dp);

  const G4Track * theIncomingTrack = fParticleChangeForGamma->GetCurrentTrack();
      G4DNAChemistryManager::Instance()->CreateWaterMolecule(eIonizedMolecule,
          1,
          theIncomingTrack);

}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

G4int TsDNAKondoRamosLiChargeDecrease::NumberOfFinalStates(G4ParticleDefinition* particleDefinition,
                                                              G4int finalStateIndex)

{
  if (particleDefinition == lithium3pDef ||
      particleDefinition == lithium2pDef ||
      particleDefinition == lithium1pDef )
    return 1;
  return 0;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

G4ParticleDefinition* TsDNAKondoRamosLiChargeDecrease::OutgoingParticleDefinition(G4ParticleDefinition* particleDefinition,
                                                                                     G4int finalStateIndex)
{
  if (particleDefinition == lithium3pDef)
    return lithium2pDef;

  else if (particleDefinition == lithium2pDef)
    return lithium1pDef;

  else if (particleDefinition == lithium1pDef)
    return lithium0pDef;

  return 0;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

G4double TsDNAKondoRamosLiChargeDecrease::WaterBindingEnergyConstant(G4ParticleDefinition* particleDefinition,
                                                                        G4int finalStateIndex)
{
  return WaterBindingEnergy[finalStateIndex];
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

G4double TsDNAKondoRamosLiChargeDecrease::OutgoingParticleBindingEnergyConstant(G4ParticleDefinition* particleDefinition,
                                                                                   G4int finalStateIndex)
{
  // Li3+ -> None :   0.00 eV
  // Li2+ -> Li3+ : 122.45 eV
  // Li1+ -> Li2+ : 75.637 eV
  // Li0+ -> Li1+ :   5.39 eV
  return IonIonizationEnergy[particleDefinition];

  return 0;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

G4double TsDNAKondoRamosLiChargeDecrease::PartialCrossSection(G4double k,
                                                                 G4int index,
                                                                 const G4ParticleDefinition* particleDefinition)
{
  G4double crossSection = CrossSections[particleDefinition][0][index+1];

  // Linear interpolation for energies in between
  for (size_t i = 0; i < CrossSections[particleDefinition].size() - 1; i++) {
    G4double Energy1 = CrossSections[particleDefinition][i+0][0];
    G4double Energy2 = CrossSections[particleDefinition][i+1][0];
    G4double Cross1  = CrossSections[particleDefinition][i+0][index+1];
    G4double Cross2  = CrossSections[particleDefinition][i+1][index+1];

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

G4int TsDNAKondoRamosLiChargeDecrease::RandomSelect(G4double k,
                                                       const G4ParticleDefinition* particleDefinition)
{
  const G4int n = 5;
  G4double* values(new G4double[n]);
  G4double value(0);
  G4int i = n;

  while (i > 0)
  {
    i--;
    values[i] = PartialCrossSection(k, i, particleDefinition);
    value += values[i];
  }

  value *= G4UniformRand();

  i = n;
  while (i > 0)
  {
    i--;

    if (values[i] > value)
      break;

    value -= values[i];
  }

  delete[] values;

  return i;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

G4double TsDNAKondoRamosLiChargeDecrease::Sum(G4double k,const G4ParticleDefinition* particleDefinition)
{
  G4double totalCrossSection = 0.;

  for (G4int i = 0; i < 5; i++)
  {
    totalCrossSection += PartialCrossSection(k, i, particleDefinition);
  }

  return totalCrossSection;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......