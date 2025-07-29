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

#include "TsDNAKRIonMillerGreenExcitationModel.hh"
#include "G4SystemOfUnits.hh"
#include "G4DNAChemistryManager.hh"
#include "G4DNAMolecularMaterial.hh"
#include "G4Exp.hh"

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

using namespace std;

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

TsDNAKRIonMillerGreenExcitationModel::TsDNAKRIonMillerGreenExcitationModel(const G4ParticleDefinition*,
                                                                 const G4String& nam)
:G4VEmModel(nam),isInitialised(false)
{
  fpMolWaterDensity = 0;

  nLevels=0;

  fLowEnergyLimit  = 0 * eV;
  fHighEnergyLimit = 0 * keV;

  verboseLevel= 1;
  // Verbosity scale:
  // 0 = nothing
  // 1 = warning for energy non-conservation
  // 2 = details of energy budget
  // 3 = calculation of cross sections, file openings, sampling of atoms
  // 4 = entering in methods

  if( verboseLevel>0 ) {
    G4cout << "KR Miller & Green excitation model is constructed " << G4endl;
  }
  fParticleChangeForGamma = 0;

  // Selection of stationary mode

  statCode = false;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

TsDNAKRIonMillerGreenExcitationModel::~TsDNAKRIonMillerGreenExcitationModel()
{}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

void TsDNAKRIonMillerGreenExcitationModel::Initialise(const G4ParticleDefinition* particle,
                                                      const G4DataVector& /*cuts*/)
{

  if (verboseLevel > 3)
    G4cout << "Calling TsDNAKRIonMillerGreenExcitationModel::Initialise()" << G4endl;

  // Energy limits

  G4DNAGenericIonsManager *instance;
  instance = G4DNAGenericIonsManager::Instance();

  fAtomicMassNumber   = particle->GetAtomicMass();
  fAtomicNumber       = particle->GetLeptonNumber();

  G4double Observers  = fAtomicMassNumber - particle->GetPDGCharge()/eplus;
  fQeff               = IonSlaterEffectiveCharge(fAtomicMassNumber,Observers);

  G4String proton;

  // LIMITS AND CONSTANTS
  fLowEnergyLimit  = fAtomicMassNumber * 10. * eV;
  fHighEnergyLimit = fAtomicMassNumber * 500. * keV;

  SetLowEnergyLimit(fLowEnergyLimit);
  SetHighEnergyLimit(fHighEnergyLimit);

  nLevels = waterExcitation.NumberOfLevels();

  //
  if( verboseLevel>0 )
  {
    G4cout << "KR Miller & Green excitation model is initialized " << G4endl
           << "Energy range: "
           << LowEnergyLimit() / eV << " eV - "
           << HighEnergyLimit() / keV << " keV for "
           << particle->GetParticleName()
           << G4endl;
  }

  // Initialize water density pointer
  fpMolWaterDensity = G4DNAMolecularMaterial::Instance()->GetNumMolPerVolTableFor(G4Material::GetMaterial("G4_WATER"));

  if (isInitialised) { return; }
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

G4double TsDNAKRIonMillerGreenExcitationModel::CrossSectionPerVolume(const G4Material* material,
                                                                const G4ParticleDefinition* particleDefinition,
                                                                G4double k,
                                                                G4double,
                                                                G4double)
{
  if (verboseLevel > 3)
    G4cout << "Calling CrossSectionPerVolume() of TsDNAKRIonMillerGreenExcitationModel" << G4endl;

  // Calculate total cross section for model

  G4DNAGenericIonsManager *instance;
  instance = G4DNAGenericIonsManager::Instance();

  if (particleDefinition != G4Proton::ProtonDefinition()   &&
      particleDefinition != instance->GetIon("lithium+")   &&
      particleDefinition != instance->GetIon("lithium++")  &&
      particleDefinition != instance->GetIon("lithium+++") &&
      particleDefinition != instance->GetIon("lithium"))
     return 0;

  G4double lowLim  = fLowEnergyLimit;
  G4double highLim = fHighEnergyLimit;
  G4double crossSection = 0.;

  G4double waterDensity = (*fpMolWaterDensity)[material->GetIndex()];

  const G4String& particleName = particleDefinition->GetParticleName();

  if (k >= lowLim && k <= highLim) {
    crossSection = Sum(k,particleDefinition);    
  }

  if (verboseLevel > 2)
  {
    G4cout << "__________________________________" << G4endl;
    G4cout << "TsDNAKRIonMillerGreenExcitationModel - XS INFO START" << G4endl;
    G4cout << "Kinetic energy(eV)=" << k/eV << " particle : " << particleDefinition->GetParticleName() << G4endl;
    G4cout << "Cross section per water molecule (cm^2)=" << crossSection/cm/cm << G4endl;
    G4cout << "Cross section per water molecule (cm^-1)=" << crossSection*waterDensity/(1./cm) << G4endl;
    // G4cout << " - Cross section per water molecule (cm^-1)=" << sigma*material->GetAtomicNumDensityVector()[1]/(1./cm) << G4endl;
    G4cout << "TsDNAKRIonMillerGreenExcitationModel - XS INFO END" << G4endl;
  }

    return crossSection*waterDensity;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

void TsDNAKRIonMillerGreenExcitationModel::SampleSecondaries(std::vector<G4DynamicParticle*>* /*fvect*/,
                                                        const G4MaterialCutsCouple* /*couple*/,
                                                        const G4DynamicParticle* aDynamicParticle,
                                                        G4double,
                                                        G4double)
{

  if (verboseLevel > 3)
    G4cout << "Calling SampleSecondaries() of TsDNAKRIonMillerGreenExcitationModel" << G4endl;

  G4double particleEnergy0 = aDynamicParticle->GetKineticEnergy();

  G4int level = RandomSelect(particleEnergy0,aDynamicParticle->GetDefinition());

  // Dingfelder's excitation levels
  const G4double excitation[]={ 8.17*eV, 10.13*eV, 11.31*eV, 12.91*eV, 14.50*eV};
  G4double excitationEnergy = excitation[level];

  G4double newEnergy = 0.;
  
  if (!statCode) newEnergy = particleEnergy0 - excitationEnergy;

  else newEnergy = particleEnergy0;
  
  if (newEnergy>0)
  {
    fParticleChangeForGamma->ProposeMomentumDirection(aDynamicParticle->GetMomentumDirection());
    fParticleChangeForGamma->SetProposedKineticEnergy(newEnergy);
    fParticleChangeForGamma->ProposeLocalEnergyDeposit(excitationEnergy);

    const G4Track * theIncomingTrack = fParticleChangeForGamma->GetCurrentTrack();
    G4DNAChemistryManager::Instance()->CreateWaterMolecule(eExcitedMolecule,
     level, theIncomingTrack);

  }

}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

G4double TsDNAKRIonMillerGreenExcitationModel::GetPartialCrossSection(const G4Material*,
                                          G4int level,
                                          const G4ParticleDefinition* particleDefinition,
                                          G4double kineticEnergy)
{
  return PartialCrossSection(kineticEnergy, level, particleDefinition);
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

G4double TsDNAKRIonMillerGreenExcitationModel::PartialCrossSection(G4double k, G4int excitationLevel, 
                                                              const G4ParticleDefinition* particleDefinition)
{
  //                               ( ( z * aj ) ^ omegaj ) * ( t - ej ) ^ nu
  // sigma(t) = zEff^2 * sigma0 * --------------------------------------------
  //                               jj ^ ( omegaj + nu ) + t ^ ( omegaj + nu )
  //
  // where t is the kinetic energy corrected by Helium mass over proton mass for Helium ions
  //
  // zEff is:
  //  1 for protons
  //  2 for alpha++
  //  and  2 - c1 S_1s - c2 S_2s - c3 S_2p for alpha+ and He
  //
  // Dingfelder et al., RPC 59, 255-275, 2000 from Miller and Green (1973)
  // Formula (34) and Table 2

  const G4double sigma0(1.E+8 * barn);
  const G4double nu(1.);
  const G4double aj[]={876.*eV, 2084.* eV, 1373.*eV, 692.*eV, 900.*eV};
  const G4double jj[]={19820.*eV, 23490.*eV, 27770.*eV, 30830.*eV, 33080.*eV};
  const G4double omegaj[]={0.85, 0.88, 0.88, 0.78, 0.78};

  // Dingfelder's excitation levels
  const G4double Eliq[5]={ 8.17*eV, 10.13*eV, 11.31*eV, 12.91*eV, 14.50*eV};

  G4int particleTypeIndex = 0;
  G4DNAGenericIonsManager* instance;
  instance = G4DNAGenericIonsManager::Instance();

  particleTypeIndex = 0;

  G4double tCorrected;
  tCorrected = k;

  // SI - added protection
  if (tCorrected < Eliq[excitationLevel]) return 0;
  //

  G4int z = 10;

  G4double numerator = std::pow(z * aj[excitationLevel], omegaj[excitationLevel]) *
                       std::pow(tCorrected - Eliq[excitationLevel], nu);

  G4double power = omegaj[excitationLevel] + nu;

  G4double denominator = std::pow(jj[excitationLevel], power) + 
                         std::pow(tCorrected, power);

  G4double zEff = IonEffectiveCharge(particleDefinition,k,Eliq[excitationLevel]);

  G4double cross = sigma0 * zEff * zEff * numerator / denominator;


  return cross;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

G4int TsDNAKRIonMillerGreenExcitationModel::RandomSelect(G4double k,const G4ParticleDefinition* particle)
{
  G4int i = nLevels;
  G4double value = 0.;
  std::deque<G4double> values;

  G4DNAGenericIonsManager *instance;
  instance = G4DNAGenericIonsManager::Instance();

  if ( particle == instance->GetIon("lithium")   ||
       particle == instance->GetIon("lithium+")  ||
       particle == instance->GetIon("lithium++") ||
       particle == instance->GetIon("lithium++")) {
    while (i > 0) {
      i--;
      G4double partial = PartialCrossSection(k,i,particle);
      values.push_front(partial);
      value += partial;
    }

    value *= G4UniformRand();

    i = nLevels;

    while (i > 0) {
      i--;
      if (values[i] > value) return i;
      value -= values[i];
    }
  }
  return 0;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

G4double TsDNAKRIonMillerGreenExcitationModel::Sum(G4double k, const G4ParticleDefinition* particle)
{
  G4double totalCrossSection = 0.;

  for (G4int i=0; i<nLevels; i++)
  {
    totalCrossSection += PartialCrossSection(k,i,particle);
  }
  return totalCrossSection;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

G4double TsDNAKRIonMillerGreenExcitationModel::IonEffectiveCharge(const G4ParticleDefinition* particle, G4double k, G4double energyTransferred) {
    G4int Z = fAtomicNumber;
    G4int Q = particle->GetPDGCharge() / eplus;
    G4int N = Z-Q;
    G4double Mass = particle->GetPDGMass();

    if (N == 0) return Z;

    //G4double Qeff = IonSlaterEffectiveCharge(Z,N);
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

bool MillerGreenCompareAtomicShell(std::vector<G4int> orb1,std::vector<G4int> orb2) {
    if (orb1[0] < orb2[0])
        return true;
    else if ((orb1[0] == orb2[0]) and orb1[1] < orb2[1])
        return true;

    return false;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

G4double TsDNAKRIonMillerGreenExcitationModel::IonSlaterEffectiveCharge(G4int AtomicNumber, G4int NumberOfElectrons) {
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
    std::sort(OrbitalLb.begin(),OrbitalLb.end(),MillerGreenCompareAtomicShell);

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