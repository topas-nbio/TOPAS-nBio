// Extra Class for TsEmDNAChemistry
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
// Author: Mathieu Karamitros
//
// WARNING : This class is released as a prototype.
// It might strongly evolve or even disapear in the next releases.
//
// History:
// -----------
// 10 Oct 2011 M.Karamitros created
//
// -------------------------------------------------------------------

#include "TsDNAWaterDissociationDisplacer.hh"
#include "G4PhysicalConstants.hh"
#include "G4SystemOfUnits.hh"
#include "G4H2O.hh"
#include "G4H2.hh"
#include "G4Hydrogen.hh"
#include "G4Oxygen.hh"
#include "G4OH.hh"
#include "G4H3O.hh"
#include "G4Electron_aq.hh"
#include "G4H2O2.hh"
#include "Randomize.hh"
#include "G4Molecule.hh"
#include "G4MolecularConfiguration.hh"
#include "G4RandomDirection.hh"
#include "G4Exp.hh"
#include "G4UnitsTable.hh"
#include "G4EmParameters.hh"
#include "G4DNAOneStepThermalizationModel.hh"

using namespace std;

//------------------------------------------------------------------------------

G4CT_COUNT_IMPL(TsDNAWaterDissociationDisplacer,
                Ionisation_DissociationDecay)

G4CT_COUNT_IMPL(TsDNAWaterDissociationDisplacer,
                A1B1_DissociationDecay)

G4CT_COUNT_IMPL(TsDNAWaterDissociationDisplacer,
                B1A1_DissociationDecay)

G4CT_COUNT_IMPL(TsDNAWaterDissociationDisplacer,
                B1A1_DissociationDecay2)

G4CT_COUNT_IMPL(TsDNAWaterDissociationDisplacer,
                AutoIonisation)

G4CT_COUNT_IMPL(TsDNAWaterDissociationDisplacer,
                DissociativeAttachment)

TsDNAWaterDissociationDisplacer::TsDNAWaterDissociationDisplacer():
G4VMolecularDissociationDisplacer(),
ke(1.7*eV)
{
	dnaSubType = G4EmParameters::Instance()->DNAeSolvationSubType();
}

//------------------------------------------------------------------------------

TsDNAWaterDissociationDisplacer::~TsDNAWaterDissociationDisplacer() {;}

//------------------------------------------------------------------------------

G4ThreeVector TsDNAWaterDissociationDisplacer::GetMotherMoleculeDisplacement(const G4MolecularDissociationChannel* theDecayChannel) const
{
    G4int decayType = theDecayChannel->GetDisplacementType();
    G4double RMSMotherMoleculeDisplacement(0.);

    switch (decayType) {
        case Ionisation_DissociationDecay:
            RMSMotherMoleculeDisplacement = 1.24 * nanometer;
            break;
        case A1B1_DissociationDecay:
            RMSMotherMoleculeDisplacement = 0. * nanometer;
            break;
        case B1A1_DissociationDecay:
            RMSMotherMoleculeDisplacement = 0. * nanometer;
            break;
        case B1A1_DissociationDecay2:
            RMSMotherMoleculeDisplacement = 0. * nanometer;
            break;
        case AutoIonisation:
            RMSMotherMoleculeDisplacement = 1.24 * nanometer;
            break;
        case DissociativeAttachment:
            RMSMotherMoleculeDisplacement = 0. * nanometer;
            break;    
    }

    if (RMSMotherMoleculeDisplacement == 0)
    {
        return G4ThreeVector(0, 0, 0);
    }
    auto RandDirection = radialDistributionOfProducts(RMSMotherMoleculeDisplacement);
    return RandDirection;
}

//------------------------------------------------------------------------------

vector<G4ThreeVector> TsDNAWaterDissociationDisplacer::GetProductsDisplacement(const G4MolecularDissociationChannel* pDecayChannel) const
{
    G4int nbProducts = pDecayChannel->GetNbProducts();
    vector<G4ThreeVector> theProductDisplacementVector(nbProducts);

    typedef map<const G4MoleculeDefinition*, G4double> RMSmap;
    RMSmap theRMSmap;

    G4int decayType = pDecayChannel->GetDisplacementType();

    G4double theOHRMSDisplacement = 0.62 * nanometer;
    G4double theH3OOHDisplacement = 0.62 * nanometer; // significant
    G4double theOH2Displacement   = 0.62 * nanometer; // not significant at lowLET
    G4double theHOHDisplacement   = 0.62 * nanometer; // significant at lowLET

    switch (decayType)
    {
        case Ionisation_DissociationDecay:
        {
            if (fVerbose) {
                G4cout << "Ionisation_DissociationDecay" << G4endl;
                G4cout << "Channel's name: " << pDecayChannel->GetName() << G4endl;
            }
            G4double RdmValue = G4UniformRand();

            if (RdmValue < 0.5) {
                // H3O
                theRMSmap[G4H3O::Definition()] = 0. * nanometer;
                // OH
                theRMSmap[G4OH::Definition()] = theH3OOHDisplacement;
            } else {
                // H3O
                theRMSmap[G4H3O::Definition()] = theH3OOHDisplacement;
                // OH
                theRMSmap[G4OH::Definition()] = 0. * nanometer;
            }

            for (int i = 0; i < nbProducts; i++) {
                auto pProduct = pDecayChannel->GetProduct(i);
                G4double theRMSDisplacement = theRMSmap[pProduct->GetDefinition()];

                if (theRMSDisplacement == 0.) {
                    theProductDisplacementVector[i] = G4ThreeVector();
                } else {
                    auto RandDirection = radialDistributionOfProducts(theRMSDisplacement);
                    theProductDisplacementVector[i] = RandDirection;
                }
            }
            break;
        }
        case A1B1_DissociationDecay:
        {
            if (fVerbose) {
                G4cout << "A1B1_DissociationDecay" << G4endl;
                G4cout << "Channel's name: " << pDecayChannel->GetName() << G4endl;
            }
            G4double theRMSDisplacement = theHOHDisplacement;
            auto RandDirection = radialDistributionOfProducts(theRMSDisplacement);

            for (G4int i = 0; i < nbProducts; i++) {
                auto pProduct = pDecayChannel->GetProduct(i);

                if (pProduct->GetDefinition() == G4OH::Definition()) {
                    theProductDisplacementVector[i] = -1. / 18. * RandDirection;
                } else if (pProduct->GetDefinition() == G4Hydrogen::Definition()) {
                    theProductDisplacementVector[i] = +17. / 18. * RandDirection;
                }
            }
            break;
        }
        case B1A1_DissociationDecay:
        {
            if (fVerbose) {
                G4cout << "B1A1_DissociationDecay" << G4endl;
                G4cout << "Channel's name: " << pDecayChannel->GetName() << G4endl;
            }

            G4double theRMSDisplacement = theOH2Displacement;
            auto RandDirection = radialDistributionOfProducts(theRMSDisplacement);

            G4int NbOfOH = 0;
            for (G4int i = 0; i < nbProducts; ++i) {
                auto pProduct = pDecayChannel->GetProduct(i);
                if (pProduct->GetDefinition() == G4H2::Definition()) {
                    // H2
                    theProductDisplacementVector[i] = -16. / 18. * RandDirection;
                } else if (pProduct->GetDefinition() == G4OH::Definition()) {
                    // OH
                    G4ThreeVector OxygenDisplacement = +2. / 18. * RandDirection;
                    G4double OHRMSDisplacement = 2.0/3.0*theOHRMSDisplacement;

                    auto OHDisplacement = radialDistributionOfProducts(OHRMSDisplacement);

                    if (NbOfOH == 0) { OHDisplacement = 0.5 * OHDisplacement; }
                    else { OHDisplacement = -0.5 * OHDisplacement; }

                    theProductDisplacementVector[i] = OHDisplacement + OxygenDisplacement;

                    ++NbOfOH;
                }
            }
            break;
        }
        case B1A1_DissociationDecay2:
        {
            if(fVerbose) {
                G4cout<<"B1A1_DissociationDecay2"<<G4endl;
                G4cout<<"Channel's name: "<<pDecayChannel->GetName()<<G4endl;
            }

            G4int NbOfH = 0;
            for(G4int i =0; i < nbProducts; ++i) {
                auto pProduct = pDecayChannel->GetProduct(i);
                if(pProduct->GetDefinition() == G4Oxygen::Definition()) {
                    // O(3p)
                    theProductDisplacementVector[i] = G4ThreeVector(0,0,0);
                } else if(pProduct->GetDefinition() == G4Hydrogen::Definition()){
                    // H
                    G4double HRMSDisplacement = 1.6 * nanometer;

                    auto HDisplacement = radialDistributionOfProducts(HRMSDisplacement);

                    if(NbOfH==0) { HDisplacement = 0.5*HDisplacement; }
                    else         { HDisplacement = -0.5*HDisplacement; }
                    theProductDisplacementVector[i] = HDisplacement;

                    ++NbOfH;
                }
            }
            break;
        }
        case AutoIonisation:
        {
            if (fVerbose) {
                G4cout << "AutoIonisation" << G4endl;
                G4cout << "Channel's name: " << pDecayChannel->GetName() << G4endl;
            }

            G4double RdmValue = G4UniformRand();

            if (RdmValue < 0.5) {
                // H3O
                theRMSmap[G4H3O::Definition()] = 0. * nanometer;
                // OH
                theRMSmap[G4OH::Definition()] = theH3OOHDisplacement;
            } else {
                // H3O
                theRMSmap[G4H3O::Definition()] = theH3OOHDisplacement;
                // OH
                theRMSmap[G4OH::Definition()] = 0. * nanometer;
            }

            for (G4int i = 0; i < nbProducts; i++) {
                auto pProduct = pDecayChannel->GetProduct(i);
                auto theRMSDisplacement = theRMSmap[pProduct->GetDefinition()];

                if (theRMSDisplacement == 0) {
                    theProductDisplacementVector[i] = G4ThreeVector();
                } else {
                    auto RandDirection = radialDistributionOfProducts(theRMSDisplacement);
                    theProductDisplacementVector[i] = RandDirection;
                }
                if (pProduct->GetDefinition() == G4Electron_aq::Definition()) {
                    theProductDisplacementVector[i] = radialDistributionOfElectron();
                }
            }
            break;
        }
        case DissociativeAttachment:
        {
            if (fVerbose) {
                G4cout << "DissociativeAttachment" << G4endl;
                G4cout << "Channel's name: " << pDecayChannel->GetName() << G4endl;
            }
            G4double theRMSDisplacement = theOH2Displacement; //0.93 * nanometer;
            auto RandDirection = radialDistributionOfProducts(theRMSDisplacement);

            G4int NbOfOH = 0;
            for (G4int i = 0; i < nbProducts; ++i) {
                auto pProduct = pDecayChannel->GetProduct(i);
                if (pProduct->GetDefinition() == G4H2::Definition()) {
                    theProductDisplacementVector[i] = -16. / 18. * RandDirection;
                }
                else if (pProduct->GetDefinition() == G4OH::Definition()) {
                    G4ThreeVector OxygenDisplacement = +2. / 18. * RandDirection;
                    G4double OHRMSDisplacement = 2./3.*theOHRMSDisplacement;

                    auto OHDisplacement = radialDistributionOfProducts(OHRMSDisplacement);

                    if (NbOfOH == 0) {
                        OHDisplacement = 0.5 * OHDisplacement;
                    }
                    else {
                        OHDisplacement = -0.5 * OHDisplacement;
                    }
                    theProductDisplacementVector[i] = OHDisplacement + OxygenDisplacement;
                    ++NbOfOH;
                }
            }
            break;
        }
    }
    return theProductDisplacementVector;
}

//------------------------------------------------------------------------------

G4ThreeVector TsDNAWaterDissociationDisplacer::radialDistributionOfProducts(G4double Rrms) const
{
    static const double inverse_sqrt_3 = 1. / sqrt(3.);
    double sigma = Rrms * inverse_sqrt_3;
    double x = G4RandGauss::shoot(0., sigma);
    double y = G4RandGauss::shoot(0., sigma);
    double z = G4RandGauss::shoot(0., sigma);
    return G4ThreeVector(x, y, z);
}

//------------------------------------------------------------------------------

G4ThreeVector TsDNAWaterDissociationDisplacer::radialDistributionOfElectron() const
{
	G4ThreeVector pdf = G4ThreeVector(0,0,0);

	if(dnaSubType == fRitchie1994eSolvation) DNA::Penetration::Ritchie1994::GetPenetration(ke,pdf);
	else if(dnaSubType == fTerrisol1990eSolvation) DNA::Penetration::Terrisol1990::GetPenetration(ke,pdf);
	else if(dnaSubType == fMeesungnoensolid2002eSolvation) DNA::Penetration::Meesungnoen2002_amorphous::GetPenetration(ke,pdf);
	else if(dnaSubType == fKreipl2009eSolvation) DNA::Penetration::Kreipl2009::GetPenetration(ke,pdf);
	else DNA::Penetration::Meesungnoen2002::GetPenetration(ke,pdf);

	return pdf;
}
