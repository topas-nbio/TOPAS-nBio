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

#include "TsDNAKRIonBornExcitationModel.hh"
#include "G4SystemOfUnits.hh"
#include "G4DNAChemistryManager.hh"
#include "G4DNAMolecularMaterial.hh"
#include <map>

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

using namespace std;

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

TsDNAKRIonBornExcitationModel::TsDNAKRIonBornExcitationModel(const G4ParticleDefinition*,
                                                     const G4String& nam) :
G4VEmModel(nam), isInitialised(false), fTableData(0) {
  fpMolWaterDensity = 0;
  fHighEnergy = 0;
  fLowEnergy = 0;
  fParticleDefinition = 0;

  verboseLevel = 1;
  // Verbosity scale:
  // 0 = nothing
  // 1 = warning for energy non-conservation
  // 2 = details of energy budget
  // 3 = calculation of cross sections, file openings, sampling of atoms
  // 4 = entering in methods

  if (verboseLevel > 0) {
    G4cout << "KR Born excitation model for ions is constructed " << G4endl;
  }
  fParticleChangeForGamma = 0;

  // Selection of stationary mode

  statCode = false;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

TsDNAKRIonBornExcitationModel::~TsDNAKRIonBornExcitationModel() {
  // Cross section
  if (fTableData)
    delete fTableData;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

void TsDNAKRIonBornExcitationModel::Initialise(const G4ParticleDefinition* particle,const G4DataVector&) {

  if (verboseLevel > 3) {
    G4cout << "Calling TsDNAKRIonBornExcitationModel::Initialise()" << G4endl;
  }

  fParticleDefinition = particle;
  fAtomicMassNumber   = particle->GetAtomicMass();
  fAtomicNumber       = particle->GetLeptonNumber();

  G4double Observers  = fAtomicMassNumber - particle->GetPDGCharge()/eplus;
  fQeff               = IonSlaterEffectiveCharge(fAtomicMassNumber,Observers);

  fTableFile = "dna/sigma_excitation_p_born";
  fLowEnergy  = fAtomicMassNumber * 500. * keV;
  fHighEnergy = fAtomicMassNumber * 100. * MeV;

  SetLowEnergyLimit(fLowEnergy);
  SetHighEnergyLimit(fHighEnergy);

  G4double scaleFactor = (1.e-22 / 3.343) * m*m;
  fTableData = new G4DNACrossSectionDataSet(new G4LogLogInterpolation, fAtomicMassNumber*eV,scaleFactor );
  fTableData->LoadData(fTableFile);

  if( verboseLevel>0 )
  {
    G4cout << "KR Born excitation model is initialized " << G4endl
    << "Energy range: "
    << LowEnergyLimit() / eV << " eV - "
    << HighEnergyLimit() / keV << " keV for "
    << particle->GetParticleName()
    << G4endl;
  }

  // Initialize water density pointer
  fpMolWaterDensity = G4DNAMolecularMaterial::Instance()->GetNumMolPerVolTableFor(G4Material::GetMaterial("G4_WATER"));

  if (isInitialised)
  { return;}
  fParticleChangeForGamma = GetParticleChangeForGamma();
  isInitialised = true;

  EFFECTIVECHART.push_back({1.000, 0.0000, 0.000, 0.0000});
  EFFECTIVECHART.push_back({2.625, 1.2996, 1.770, 1.1402});
  EFFECTIVECHART.push_back({2.164, 0.9764, 1.750, 0.6821});
  EFFECTIVECHART.push_back({1.300, 0.6465, 1.880, 0.5547});
  EFFECTIVECHART.push_back({1.031, 0.4924, 2.000, 0.4939});
  EFFECTIVECHART.push_back({1.065, 0.4800, 2.130, 0.4434});
  EFFECTIVECHART.push_back({1.179, 0.4677, 2.270, 0.4143});
  EFFECTIVECHART.push_back({1.360, 0.4613, 2.410, 0.3925});
  EFFECTIVECHART.push_back({1.508, 0.4602, 2.590, 0.3755});
  EFFECTIVECHART.push_back({1.792, 0.4515, 2.710, 0.3671});
  EFFECTIVECHART.push_back({1.712, 0.3923, 2.850, 0.3469});
  EFFECTIVECHART.push_back({1.492, 0.3452, 3.010, 0.3269});
  EFFECTIVECHART.push_back({1.170, 0.3191, 3.170, 0.3087});
  EFFECTIVECHART.push_back({1.012, 0.2933, 3.260, 0.2958});
  EFFECTIVECHART.push_back({0.954, 0.2659, 3.330, 0.2857});
  EFFECTIVECHART.push_back({0.926, 0.2478, 3.392, 0.2739});
  EFFECTIVECHART.push_back({0.933, 0.2368, 3.447, 0.2633});
  EFFECTIVECHART.push_back({0.957, 0.2165, 3.500, 0.2560});
  EFFECTIVECHART.push_back({0.964, 0.2151, 3.516, 0.2509});
  EFFECTIVECHART.push_back({0.941, 0.2248, 3.570, 0.2404});
  EFFECTIVECHART.push_back({0.950, 0.2324, 3.627, 0.2328});
  EFFECTIVECHART.push_back({0.998, 0.2345, 3.667, 0.2238});
  EFFECTIVECHART.push_back({1.061, 0.2243, 3.709, 0.2171});
  EFFECTIVECHART.push_back({1.138, 0.2291, 3.745, 0.2187});
  EFFECTIVECHART.push_back({1.207, 0.2408, 3.803, 0.2090});
  EFFECTIVECHART.push_back({1.308, 0.2391, 3.840, 0.2088});
  EFFECTIVECHART.push_back({1.397, 0.2462, 3.891, 0.2048});
  EFFECTIVECHART.push_back({1.455, 0.2397, 3.973, 0.1925});
  EFFECTIVECHART.push_back({1.520, 0.2246, 4.000, 0.1985});
  EFFECTIVECHART.push_back({1.538, 0.2106, 4.050, 0.1878});
  EFFECTIVECHART.push_back({1.541, 0.1988, 4.110, 0.2001});
  EFFECTIVECHART.push_back({1.512, 0.1914, 4.182, 0.1897});
  EFFECTIVECHART.push_back({1.492, 0.1990, 4.230, 0.1782});
  EFFECTIVECHART.push_back({1.460, 0.1857, 4.290, 0.1772});
  EFFECTIVECHART.push_back({1.407, 0.1897, 4.369, 0.1686});
  EFFECTIVECHART.push_back({1.351, 0.1872, 4.418, 0.1611});
  EFFECTIVECHART.push_back({1.286, 0.1686, 4.494, 0.1619});
  EFFECTIVECHART.push_back({1.129, 0.1784, 4.618, 0.1509});
  EFFECTIVECHART.push_back({1.139, 0.1702, 4.680, 0.1485});
  EFFECTIVECHART.push_back({1.136, 0.1694, 4.749, 0.1412});
  EFFECTIVECHART.push_back({1.197, 0.1601, 4.769, 0.1435});
  EFFECTIVECHART.push_back({1.246, 0.1587, 4.829, 0.1397});
  EFFECTIVECHART.push_back({1.205, 0.1358, 4.904, 0.1414});
  EFFECTIVECHART.push_back({1.130, 0.1395, 4.990, 0.1324});
  EFFECTIVECHART.push_back({1.050, 0.1354, 5.050, 0.1314});
  EFFECTIVECHART.push_back({1.044, 0.1107, 5.101, 0.1316});
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

G4double TsDNAKRIonBornExcitationModel::CrossSectionPerVolume(const G4Material* material,
                                                          const G4ParticleDefinition* particleDefinition,
                                                          G4double ekin,
                                                          G4double,
                                                          G4double)
{
  if (verboseLevel > 3) {
    G4cout << "Calling CrossSectionPerVolume() of TsDNAKRIonBornExcitationModel"
        << G4endl;
  }

  if(particleDefinition != fParticleDefinition) return 0;

  // Calculate total cross section for model

  G4double sigma=0;

  G4double waterDensity = (*fpMolWaterDensity)[material->GetIndex()];

  if (ekin >= fLowEnergy && ekin <= fHighEnergy) {
    G4int level = RandomSelect(ekin);
    G4double excitationEnergy = waterStructure.ExcitationEnergy(level);
    G4double Zeff = IonEffectiveCharge(particleDefinition,ekin,excitationEnergy);
    sigma = fTableData->FindValue(ekin);
    sigma *= Zeff*Zeff;
  }

  if (verboseLevel > 2) {
    G4cout << "__________________________________" << G4endl;
    G4cout << "TsDNAKRIonBornExcitationModel - XS INFO START" << G4endl;
    G4cout << "Kinetic energy(eV)=" << ekin/eV << " particle : " << particleDefinition->GetParticleName() << G4endl;
    G4cout << "Cross section per water molecule (cm^2)=" << sigma/cm/cm << G4endl;
    G4cout << "Cross section per water molecule (cm^-1)=" << sigma*waterDensity/(1./cm) << G4endl;
    G4cout << "TsDNAKRIonBornExcitationModel - XS INFO END" << G4endl;
  }

  return sigma*waterDensity;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

void TsDNAKRIonBornExcitationModel::SampleSecondaries(std::vector<G4DynamicParticle*>* /*fvect*/,
                                                  const G4MaterialCutsCouple* /*couple*/,
                                                  const G4DynamicParticle* aDynamicParticle,
                                                  G4double,
                                                  G4double)
{

  if (verboseLevel > 3)
  {
    G4cout << "Calling SampleSecondaries() of TsDNAKRIonBornExcitationModel"
        << G4endl;
  }

  G4double k = aDynamicParticle->GetKineticEnergy();

  G4int level = RandomSelect(k);
  G4double excitationEnergy = waterStructure.ExcitationEnergy(level);
  G4double newEnergy = k - excitationEnergy;

  if (newEnergy > 0)
  {
    fParticleChangeForGamma->ProposeMomentumDirection(aDynamicParticle->GetMomentumDirection());
    
    if (!statCode) fParticleChangeForGamma->SetProposedKineticEnergy(newEnergy);
    else fParticleChangeForGamma->SetProposedKineticEnergy(k);
    
    fParticleChangeForGamma->ProposeLocalEnergyDeposit(excitationEnergy);
  }

  const G4Track * theIncomingTrack = fParticleChangeForGamma->GetCurrentTrack();
  G4DNAChemistryManager::Instance()->CreateWaterMolecule(eExcitedMolecule,
      level,
      theIncomingTrack);
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

G4double TsDNAKRIonBornExcitationModel::GetPartialCrossSection(const G4Material*,
                                                           G4int level,
                                                           const G4ParticleDefinition* particle,
                                                           G4double kineticEnergy)
{
  if (fParticleDefinition != particle)
  {
    G4Exception("TsDNAKRIonBornExcitationModel::GetPartialCrossSection",
                "bornParticleType",
                FatalException,
                "Model initialized for another particle type.");
  }

  return fTableData->GetComponent(level)->FindValue(kineticEnergy);
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

G4int TsDNAKRIonBornExcitationModel::RandomSelect(G4double k)
{
  G4int level = 0;

  G4double* valuesBuffer = new G4double[fTableData->NumberOfComponents()];
  const size_t n(fTableData->NumberOfComponents());
  size_t i(n);
  G4double value = 0.;

  while (i > 0)
  {
    i--;
    valuesBuffer[i] = fTableData->GetComponent(i)->FindValue(k);
    value += valuesBuffer[i];
  }

  value *= G4UniformRand();
  i = n;

  while (i > 0)
  {
    i--;

    if (valuesBuffer[i] > value)
    {
      delete[] valuesBuffer;
      return i;
    }
    value -= valuesBuffer[i];
  }

  if (valuesBuffer)
    delete[] valuesBuffer;

  return level;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

G4double TsDNAKRIonBornExcitationModel::IonEffectiveCharge(const G4ParticleDefinition* particle, G4double k, G4double energyTransferred) {
    G4int Z = fAtomicNumber;
    G4int Q = particle->GetPDGCharge() / eplus;
    G4int N = Z-Q;
    G4double Mass = particle->GetPDGMass();

    if (N == 0) return Z;

    G4double tElectron = 0.511*k / Mass;
    G4double H = 2.*13.60569172 * eV;
    G4double R  = std::sqrt( 2. * tElectron / H ) / ( energyTransferred / H )  *  fQeff;

    // Lithium value 0.05, this value might be different for different Ions
    R *= 0.05;

    G4double ZE = EFFECTIVECHART[N][0]-(EFFECTIVECHART[N][1]*(Z-Q-1));
    G4double NU = EFFECTIVECHART[N][2]-(EFFECTIVECHART[N][3]*(Z-Q-1));
    G4double Screening = 0;
    if (R < 100)
        Screening = 1 / ((NU / ZE)*(exp(R * ZE) - 1) + 1);

    //G4cout << Z << "  " << N << "  " << Screening << "  " << R << "  " << energyTransferred/eV << "  " << Mass << "  " << ZE << "  " << NU << "  " << Z - (N * (1 - Screening)) << G4endl;
    return Z - (N * (1 - Screening));
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

bool BornCompareAtomicShell(std::vector<G4int> orb1,std::vector<G4int> orb2) {
    if (orb1[0] < orb2[0])
        return true;
    else if ((orb1[0] == orb2[0]) and orb1[1] < orb2[1])
        return true;

    return false;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

G4double TsDNAKRIonBornExcitationModel::IonSlaterEffectiveCharge(G4int AtomicNumber, G4int NumberOfElectrons) {
    if (NumberOfElectrons <= 1) {
        return AtomicNumber;
    }

    std::vector<G4int> Orbitals  = {2,6,10,14,18,22,26};
    std::vector<G4int> Labels    = {1,2,3,4,5,6,7};
    std::vector<std::vector<G4int>> OrbitalLb;

    // Fill Orbital
    G4int i = 0;
    G4int j = 0;
    G4int k = 0;
    G4int l = 1;
    G4int nShell    = 1;
    G4int Electrons = 0;
    while (Electrons < NumberOfElectrons) {
        Electrons += Orbitals[i];
        
        if (Electrons <= NumberOfElectrons)
            OrbitalLb.push_back({l,Labels[i],Orbitals[i]});
        else 
            OrbitalLb.push_back({l,Labels[i],Orbitals[i]-(Electrons-NumberOfElectrons)});

        if (l > nShell)
            nShell = l;

        i -= 1;
        l += 1;
        if (i < 0) {
            i = j;

            if (k == 0){
                k += 1;
            }

            else {
                k = 0;
                j -= 1;
            }

            j += 1;
            l = j+1;
        }
    }

    //for (i = 0; i < OrbitalLb.size(); i++)
    //    G4cout << OrbitalLb[i][0] << "  " << OrbitalLb[i][1] << " : " << OrbitalLb[i][2] << G4endl;
    //G4cout << G4endl;

    // Sort the orbitals and remove excess electrons
    std::sort(OrbitalLb.begin(),OrbitalLb.end(),BornCompareAtomicShell);

    //for (i = 0; i < OrbitalLb.size(); i++)
    //    G4cout << OrbitalLb[i][0] << "  " << OrbitalLb[i][1] << " : " << OrbitalLb[i][2] << G4endl;
    //G4cout << G4endl;
    //G4cout << G4endl;

    G4double S = 0;
    size_t i_inv = OrbitalLb.size();
    G4int Shell = OrbitalLb[i_inv-1][0];
    for (i = 0; i < OrbitalLb.size(); i++) {
        G4double Factor = 0.35;
        Electrons       = OrbitalLb[i_inv-i-1][2];

        if (i == 0) {
            Electrons = Electrons - 1;

            if (OrbitalLb[i_inv-1][0] == 1 and OrbitalLb[i_inv-1][1] == 1)
                Factor    = 0.30;
        }

        else if (Shell-OrbitalLb[i_inv-i-1][0] == 1) {
            if (OrbitalLb[i_inv-1][1] == 1 or OrbitalLb[i_inv-1][1] == 2)
                Factor = 0.85;
            else
                Factor = 1.0;
        }

        else if (i > 1)
            Factor = 1.0;

        S += Factor*Electrons;
    }

    //G4cout << AtomicNumber << "  " << S << "  " << Shell << "  " << AtomicNumber-S << G4endl;
    return (AtomicNumber-S)/Shell;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......