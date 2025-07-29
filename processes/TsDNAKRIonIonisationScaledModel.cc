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
// $Id: TsDNAKRIonIonisationScaledModel.cc 96060 2016-03-11 12:58:04Z gcosmo $
// GEANT4 tag $Name:  $
//
// Modified by Z. Francis, S. Incerti to handle HZE 
// && inverse rudd function sampling 26-10-2010

#include "TsDNAKRIonIonisationScaledModel.hh"
#include "G4PhysicalConstants.hh"
#include "G4SystemOfUnits.hh"
#include "G4UAtomicDeexcitation.hh"
#include "G4LossTableManager.hh"
#include "G4DNAChemistryManager.hh"
#include "G4DNAMolecularMaterial.hh"

//SEB
#include "G4IonTable.hh"
#include "G4DNARuddAngle.hh"
#include "G4DeltaAngle.hh"
#include "G4Exp.hh"
//

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

using namespace std;

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

TsDNAKRIonIonisationScaledModel::TsDNAKRIonIonisationScaledModel(const G4ParticleDefinition*,const G4String& nam)
:G4VEmModel(nam),isInitialised(false) {
    //  nistwater = G4NistManager::Instance()->FindOrBuildMaterial("G4_WATER");
    fpWaterDensity = 0;

    lowEnergyLimit  = 0 * eV;
    highEnergyLimit = 0 * eV;

    verboseLevel= 1;
    // Verbosity scale:
    // 0 = nothing
    // 1 = warning for energy non-conservation
    // 2 = details of energy budget
    // 3 = calculation of cross sections, file openings, sampling of atoms
    // 4 = entering in methods

    if( verboseLevel>0 ) {
        G4cout << "Ts Rudd-KR ionisation model is constructed " << G4endl;
    }

    // Define default angular generator
    SetAngularDistribution(new G4DNARuddAngle());

    // Mark this model as "applicable" for atomic deexcitation
    SetDeexcitationFlag(true);
    fAtomDeexcitation = 0;
    fParticleChangeForGamma = 0;

    // Selection of stationary mode

    statCode = false;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

TsDNAKRIonIonisationScaledModel::~TsDNAKRIonIonisationScaledModel() {  
    // Cross section
    delete tableData;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

void TsDNAKRIonIonisationScaledModel::Initialise(const G4ParticleDefinition* particle, const G4DataVector&) {
    if (verboseLevel > 3)
        G4cout << "Calling TsDNAKRIonIonisationScaledModel::Initialise()" << G4endl;

    // Only Proton CS needed
    G4String fileProton  = "dna/sigma_ionisation_p_rudd";
    atomicMassNumber     = particle->GetAtomicMass();
    atomicNumber         = particle->GetLeptonNumber();

    G4int Z = atomicNumber;
    G4int Q = particle->GetPDGCharge() / eplus;
    G4int N = Z-Q;

    qEff = IonSlaterEffectiveCharge(Z,N);


    // LIMITS AND DATA

    // **********************************************************************************************

    lowEnergyLimit  = atomicNumber * 100. * eV;
    highEnergyLimit = atomicNumber * 1000 * MeV;

    // Cross section

    tableData = new G4DNACrossSectionDataSet(new G4LogLogInterpolation,atomicMassNumber*eV,1*m*m );
    tableData->LoadData(fileProton);
    
    // **********************************************************************************************

    SetLowEnergyLimit(lowEnergyLimit);
    SetHighEnergyLimit(highEnergyLimit);

    //----------------------------------------------------------------------

    if( verboseLevel>0 ) {
        G4cout << "Rudd-KR ionisation model is initialized " << G4endl
               << "Energy range: "
               << LowEnergyLimit() / eV << " eV - "
               << HighEnergyLimit() / keV << " keV for "
               << particle->GetParticleName()
               << G4endl;
    }

    // Initialize water density pointer
    fpWaterDensity    = G4DNAMolecularMaterial::Instance()->GetNumMolPerVolTableFor(G4Material::GetMaterial("G4_WATER"));
    fAtomDeexcitation = G4LossTableManager::Instance()->AtomDeexcitation();

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

G4double TsDNAKRIonIonisationScaledModel::CrossSectionPerVolume(const G4Material* material,
                                                                const G4ParticleDefinition* particleDefinition,
                                                                G4double k,G4double,G4double) {
    //SEB: particleDefinition->GetParticleName() is for eg. Fe56
    //     particleDefinition->GetPDGMass() is correct
    //     particleDefinition->GetAtomicNumber() is correct

    if (verboseLevel > 3)
        G4cout << "Calling CrossSectionPerVolume() of TsDNAKRIonIonisationScaledModel" << G4endl;

    // Calculate total cross section for model

    G4DNAGenericIonsManager *instance;
    instance = G4DNAGenericIonsManager::Instance();

    if (particleDefinition != instance->GetIon("lithium+++") &&
        particleDefinition != instance->GetIon("lithium++")  &&
        particleDefinition != instance->GetIon("lithium+")   &&
        particleDefinition != instance->GetIon("lithium")    &&
        //particleDefinition != instance->GetIon("beryllium")  &&
        particleDefinition != instance->GetIon("boron")      &&
        particleDefinition != instance->GetIon("carbon")     &&
        particleDefinition != instance->GetIon("nitrogen")   &&
        particleDefinition != instance->GetIon("oxygen")     &&
        particleDefinition != instance->GetIon("iron")       &&
        particleDefinition != instance->GetIon("neon")       &&
        particleDefinition != instance->GetIon("argon"))
        return 0;

    G4double highLim = 0;
    G4double sigma=0;

    G4double waterDensity = (*fpWaterDensity)[material->GetIndex()];

    if(waterDensity!= 0.0) {
        const G4String& particleName = particleDefinition->GetParticleName();

        if (k <= highEnergyLimit) {
            //SI : XS must not be zero otherwise sampling of secondaries method ignored

            if (k < lowEnergyLimit) k = lowEnergyLimit;

            sigma = tableData->FindValue(k);

            // Effective Charge Scalling for Ions
            G4ParticleDefinition* part = instance->GetIon(particleDefinition->GetParticleName());
            G4int ionizationShell      = RandomSelect(k,particleName);
            G4double energyTransferred = ProposedSampledEnergy(part,k,ionizationShell);
            G4double Zeff = IonEffectiveCharge(particleDefinition,k,energyTransferred);
            sigma *= Zeff*Zeff;

        } // if (k >= lowLim && k < highLim)

        if (verboseLevel > 3) {
            G4cout << "__________________________________" << G4endl;
            G4cout << "TsDNAKRIonIonisationScaledModel - XS INFO START" << G4endl;
            G4cout << "Kinetic energy(eV)=" << k/eV << " particle : " << particleDefinition->GetParticleName() << G4endl;
            G4cout << "Cross section per water molecule (cm^2)=" << sigma/cm/cm << G4endl;
            G4cout << "Cross section per water molecule (cm^-1)=" << sigma*waterDensity/(1./cm) << G4endl;
            //G4cout << " - Cross section per water molecule (cm^-1)=" << sigma*material->GetAtomicNumDensityVector()[1]/(1./cm) << G4endl;
            G4cout << "TsDNAKRIonIonisationScaledModel - XS INFO END" << G4endl;

        }

    }
    return sigma*waterDensity;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

void TsDNAKRIonIonisationScaledModel::SampleSecondaries(std::vector<G4DynamicParticle*>* fvect,
                                                         const G4MaterialCutsCouple* couple,
                                                         const G4DynamicParticle* particle,
                                                         G4double, G4double) {

    //SEB: particle->GetDefinition()->GetParticleName() is for eg. Fe56
    //     particle->GetDefinition()->GetPDGMass() is correct
    //     particle->GetDefinition()->GetAtomicNumber() is correct
    //     particle->GetDefinition()->GetAtomicMass() is correct

    if (verboseLevel > 3)
        G4cout << "Calling SampleSecondaries() of TsDNAKRIonIonisationScaledModel" << G4endl;

    G4double lowLim  = lowEnergyLimit;
    G4double highLim = highEnergyLimit;
    G4double k = particle->GetKineticEnergy();

    const G4String& particleName = particle->GetDefinition()->GetParticleName();

    if (k >= lowLim && k < highLim) {
        G4ParticleDefinition* definition = particle->GetDefinition();
        G4ParticleMomentum primaryDirection = particle->GetMomentumDirection();
        G4int ionizationShell = RandomSelect(k,particleName);

        // sample deexcitation
        // here we assume that H_{2}O electronic levels are the same as Oxygen.
        // this can be considered true with a rough 10% error in energy on K-shell,

        G4int secNumberInit    = 0;  // need to know at a certain point the energy of secondaries
        G4int secNumberFinal   = 0; // So I'll make the diference and then sum the energies
        G4double bindingEnergy = 0;
        bindingEnergy = waterStructure.IonisationEnergy(ionizationShell);

        //SI: additional protection if tcs interpolation method is modified
        if (k<bindingEnergy) { return; }

        G4int Z = 8;

        if(fAtomDeexcitation) {
            G4AtomicShellEnumerator as = fKShell;

            if (ionizationShell <5 && ionizationShell >1) {
                as = G4AtomicShellEnumerator(4-ionizationShell);
            }
            else if (ionizationShell <2) {
                as = G4AtomicShellEnumerator(3);
            }

            const G4AtomicShell* shell = fAtomDeexcitation->GetAtomicShell(Z, as);
            secNumberInit = fvect->size();
            fAtomDeexcitation->GenerateParticles(fvect, shell, Z, 0, 0);
            secNumberFinal = fvect->size();
        }

        G4double secondaryKinetic = RandomizeEjectedElectronEnergy(definition,k,ionizationShell);

	    G4ThreeVector deltaDirection = GetAngularDistribution()->SampleDirectionForShell(particle, secondaryKinetic, 
							                                                             Z, ionizationShell,
							                                                             couple->GetMaterial());
        fParticleChangeForGamma->ProposeMomentumDirection(primaryDirection);
        G4double scatteredEnergy = k-bindingEnergy-secondaryKinetic;
        G4double deexSecEnergy = 0;
        for (G4int j=secNumberInit; j < secNumberFinal; j++) {
            deexSecEnergy = deexSecEnergy + (*fvect)[j]->GetKineticEnergy();
        }

        if (!statCode) {
          fParticleChangeForGamma->SetProposedKineticEnergy(scatteredEnergy);
          fParticleChangeForGamma->ProposeLocalEnergyDeposit(k-scatteredEnergy-secondaryKinetic-deexSecEnergy);
        }
        else {
          fParticleChangeForGamma->SetProposedKineticEnergy(k);
          fParticleChangeForGamma->ProposeLocalEnergyDeposit(k-scatteredEnergy);
        }

        G4DynamicParticle* dp = new G4DynamicParticle (G4Electron::Electron(),deltaDirection,secondaryKinetic) ;
        fvect->push_back(dp);

        const G4Track * theIncomingTrack = fParticleChangeForGamma->GetCurrentTrack();
        G4DNAChemistryManager::Instance()->CreateWaterMolecule(eIonizedMolecule,
                                                               ionizationShell,
                                                               theIncomingTrack);
    }

    if (k < lowLim)
    {
        fParticleChangeForGamma->SetProposedKineticEnergy(0.);
        fParticleChangeForGamma->ProposeTrackStatus(fStopAndKill);
        fParticleChangeForGamma->ProposeLocalEnergyDeposit(k);
    }
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

G4double TsDNAKRIonIonisationScaledModel::RandomizeEjectedElectronEnergy(G4ParticleDefinition* particleDefinition, 
                                                                          G4double k,
                                                                          G4int shell)
{
    //-- Fast sampling method -----
    G4double proposed_energy;
    G4double random1;
    G4double value_sampling;
    G4double max1;

    do {
        proposed_energy = ProposedSampledEnergy(particleDefinition, k, shell); // Proposed energy by inverse function sampling

        max1=0.;
        G4double rej_val = 0;

        for(G4double en=0.; en<20.; en+=1.) {
            rej_val = RejectionFunction(particleDefinition, k, en, shell);
            if(rej_val > max1) {
                max1=rej_val;
            }
        }

        random1 = G4UniformRand()*max1;

        value_sampling = RejectionFunction(particleDefinition, k, proposed_energy, shell);

    } while(random1 > value_sampling);

    return(proposed_energy);
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

G4double TsDNAKRIonIonisationScaledModel::RejectionFunction(G4ParticleDefinition* particleDefinition, 
                                                             G4double k,
                                                             G4double proposed_ws,
                                                             G4int ionizationLevelIndex)
{
    const G4int j=ionizationLevelIndex;
    G4double Bj_energy, alphaConst;
    G4double Ry = 13.6*eV;
    const G4double Gj[5] = {0.99, 1.11, 1.11, 0.52, 1.};

    // Following values provided by M. Dingfelder (priv. comm)
    const G4double Bj[5] = {12.60*eV, 14.70*eV, 18.40*eV, 32.20*eV, 540*eV};

    if (j == 4) {
        alphaConst = 0.66;
        //---Note that the following (j==4) cases are provided by   M. Dingfelder (priv. comm)
        Bj_energy = waterStructure.IonisationEnergy(ionizationLevelIndex);
        //---
    }
    else {
        alphaConst = 0.64;
        Bj_energy = Bj[ionizationLevelIndex];
    }

    G4double energyTransfer = proposed_ws + Bj_energy;
    proposed_ws /= Bj_energy;

    G4double tau   = (electron_mass_c2 / particleDefinition->GetPDGMass()) * k;
    G4double A_ion = particleDefinition->GetAtomicMass();

    G4double v2;
    G4double beta2;

    // Clasical Energy Range
    if((tau/MeV)<5.447761194e-2) {
        v2    = tau / Bj_energy;
        beta2 = 2.*tau / electron_mass_c2;
    }

    // Relativistic Energy Range
    else {
        v2    = (electron_mass_c2 / 2. / Bj_energy) * (1. - (1./ pow( (1.+ (tau/electron_mass_c2)),2) ));
        beta2 = 1. - 1./(1.+ (tau/electron_mass_c2/A_ion))/(1.+ (tau/electron_mass_c2/A_ion));
    }

    G4double v  = std::sqrt(v2);
    G4double wc = 4.*v2 - 2.*v - (Ry/(4.*Bj_energy));
    G4double rejection_term = 1.+G4Exp(alphaConst*(proposed_ws - wc) / v);
    rejection_term = (1./rejection_term)*Gj[j];

    G4double Zeffion = IonEffectiveCharge(particleDefinition,k,energyTransfer);
    rejection_term*=Zeffion*Zeffion;

    return (rejection_term);
}


//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......


G4double TsDNAKRIonIonisationScaledModel::ProposedSampledEnergy(G4ParticleDefinition* particle, 
                                                                 G4double k,
                                                                 G4int ionizationLevelIndex)
{
    const G4int j=ionizationLevelIndex;

    G4double A1, B1, C1, D1, E1, A2, B2, C2, D2;
    G4double Bj_energy;

    // Following values provided by M. Dingfelder (priv. comm)
    const G4double Bj[5] = {12.60*eV, 14.70*eV, 18.40*eV, 32.20*eV, 540*eV};

    if (j == 4) {
        //Data For Liquid Water K SHELL from Dingfelder (Protons in Water)
        A1 = 1.25;
        B1 = 0.5;
        C1 = 1.00;
        D1 = 1.00;
        E1 = 3.00;
        A2 = 1.10;
        B2 = 1.30;
        C2 = 1.00;
        D2 = 0.00;
        //alphaConst = 0.66;
        //---Note that the following (j==4) cases are provided by   M. Dingfelder (priv. comm)
        Bj_energy = waterStructure.IonisationEnergy(ionizationLevelIndex);
        //---
    }
    else {
        //Data For Liquid Water from Dingfelder (Protons in Water)
        A1 = 1.02;
        B1 = 82.0;
        C1 = 0.45;
        D1 = -0.80;
        E1 = 0.38;
        A2 = 1.07;
        //B2 = 14.6; From Ding Paper
        // Value provided by M. Dingfelder (priv. comm)
        B2 = 11.6;
        //
        C2 = 0.60;
        D2 = 0.04;
        //alphaConst = 0.64;

        Bj_energy = Bj[ionizationLevelIndex];
    }

    G4double tau   = (electron_mass_c2 / particle->GetPDGMass()) * k;
    G4double A_ion = particle->GetAtomicMass();

    G4double v2;
    G4double beta2;

    // Classical Energy Range
    if((tau/MeV)<5.447761194e-2) {
        v2 = tau / Bj_energy;
        beta2 = 2.*tau / electron_mass_c2;
    }
    // Relativistic Energy Range
    else {
        v2 = (electron_mass_c2 / 2. / Bj_energy) * (1. - (1./ pow( (1.+ (tau/electron_mass_c2)),2) ));
        beta2 =1. - 1./(1.+ (tau/electron_mass_c2/A_ion))/(1.+ (tau/electron_mass_c2/A_ion));
    }

    G4double v = std::sqrt(v2);
    G4double L1 = (C1* std::pow(v,(D1))) / (1.+ E1*std::pow(v, (D1+4.)));
    G4double L2 = C2*std::pow(v,(D2));
    G4double H1 = (A1*std::log(1.+v2)) / (v2+(B1/v2));
    G4double H2 = (A2/v2) + (B2/(v2*v2));
    G4double F1 = L1+H1;
    G4double F2 = (L2*H2)/(L2+H2);

    // ZF. generalized & relativistic version
    G4double maximumEnergy;

    //---- maximum kinetic energy , non relativistic ------
    if( (k/MeV)/(particle->GetPDGMass()/MeV)  <= 0.1 ) {
        maximumEnergy = 4.* (electron_mass_c2 / particle->GetPDGMass()) * k;
    }
    //---- relativistic -----------------------------------
    else  {
        G4double gamma = 1./sqrt(1.-beta2);
        maximumEnergy = 2.*electron_mass_c2*(gamma*gamma-1.)/
                (1.+2.*gamma*(electron_mass_c2/particle->GetPDGMass())+pow(electron_mass_c2/particle->GetPDGMass(), 2.) );
    }

    //either it is transfered energy or secondary electron energy ...
    //maximumEnergy-=Bj_energy;

    //-----------------------------------------------------
    G4double wmax = maximumEnergy/Bj_energy;
    G4double c = wmax*(F2*wmax+F1*(2.+wmax))/(2.*(1.+wmax)*(1.+wmax));
    c=1./c; //!!!!!!!!!!! manual calculus leads to  c=1/c
    G4double randVal = G4UniformRand();
    G4double proposed_ws = F1*F1*c*c + 2.*F2*c*randVal - 2.*F1*c*randVal;
    proposed_ws  = -F1*c+2.*randVal+std::sqrt(proposed_ws);
    proposed_ws /= ( F1*c + F2*c - 2.*randVal );
    proposed_ws *= Bj_energy;

    return(proposed_ws);

}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

G4int TsDNAKRIonIonisationScaledModel::RandomSelect(G4double k, const G4String& )
{  
    G4double* valuesBuffer = new G4double[tableData->NumberOfComponents()];

    const size_t n(tableData->NumberOfComponents());
    size_t i(n);
    G4double value = 0.;

    while (i>0)
    {
        i--;
        valuesBuffer[i] = tableData->GetComponent(i)->FindValue(k);

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

    if (valuesBuffer) delete[] valuesBuffer;

    return 0;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

G4double TsDNAKRIonIonisationScaledModel::PartialCrossSection(const G4Track& track )
{
    G4double sigma = 0.;

    const G4DynamicParticle* particle = track.GetDynamicParticle();
    G4double k = particle->GetKineticEnergy();

    if (k >= lowEnergyLimit && k <= highEnergyLimit) {
        tableData->FindValue(k);
    }

    return sigma;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

G4double TsDNAKRIonIonisationScaledModel::Sum(G4double /* energy */, const G4String& /* particle */)
{
    return 0;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

G4double TsDNAKRIonIonisationScaledModel::IonEffectiveCharge(const G4ParticleDefinition* particle, G4double k, G4double energyTransferred) {
    G4int Z = atomicNumber;
    G4int Q = particle->GetPDGCharge() / eplus;
    G4int N = Z-Q;
    G4double Mass = particle->GetPDGMass();

    if (N == 0) return Z;

    G4double tElectron = 0.511*k / Mass;
    G4double H = 2.*13.60569172 * eV;
    G4double R  = (std::sqrt( 2. * tElectron / H ) / ( energyTransferred / H )) * qEff;

    // Lithium value 0.05, this value might be different for different Ions
    R *= 0.05;

    G4double ZE = EFFECTIVECHART[N][0]-(EFFECTIVECHART[N][1]*(Z-Q-1));
    G4double NU = EFFECTIVECHART[N][2]-(EFFECTIVECHART[N][3]*(Z-Q-1));
    G4double Screening = 0;
    if (R < 100)
        Screening = 1 / ((NU / ZE)*(exp(R * ZE) - 1) + 1);

    //G4cout << Z << "  " << N << "  " << Screening << "  " << R << "  " << tElectron/eV << "  " << energyTransferred/eV << "  " << ZE << "  " << NU << "  " << Z - (N * (1 - Screening)) << "  " << qEff <<G4endl;
    //G4cout << Z - (N * (1 - Screening)) << "  " << qEff << G4endl;
    return Z - (N * (1 - Screening));
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

bool CompareAtomicShell(std::vector<G4int> orb1,std::vector<G4int> orb2) {
    if (orb1[0] < orb2[0])
        return true;
    else if ((orb1[0] == orb2[0]) and orb1[1] < orb2[1])
        return true;

    return false;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

G4double TsDNAKRIonIonisationScaledModel::IonSlaterEffectiveCharge(G4int AtomicNumber, G4int NumberOfElectrons) {
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
    std::sort(OrbitalLb.begin(),OrbitalLb.end(),CompareAtomicShell);

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