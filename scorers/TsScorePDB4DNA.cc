// Scorer for pdb4dna
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

#include "TsScorePDB4DNA.hh"
#include "TsTrackInformation.hh"

#include "G4SystemOfUnits.hh"

#include <map>

TsScorePDB4DNA::TsScorePDB4DNA(TsParameterManager* pM, TsMaterialManager* mM, TsGeometryManager* gM, TsScoringManager* scM, TsExtensionManager* eM,
                               G4String scorerName, G4String quantity, G4String outFileName, G4bool isSubScorer)
:TsVNtupleScorer(pM, mM, gM, scM, eM, scorerName, quantity, outFileName, isSubScorer)

{
    SetUnit("");
    
    fPDBFileName = fPm->GetStringParameter(GetFullParmName("PDB4DNAFileName"));
    std::fstream in;
    in.open(fPDBFileName,std::ios::in);
    if (!in.is_open() ) {
        G4cerr << "Topas is exiting due to a serious error in scoring." << G4endl;
        G4cerr << "Input file: " << fPDBFileName << " cannot be opened for Scorer name: " << GetName() << G4endl;
        fPm->AbortSession(1);
    }
    in.close();
    
    fMoleculeList = NULL;
    fBarycenterList = NULL;
    
    G4int verbosity = 0;
    unsigned short int isProtein;
    fMoleculeList = fPDBlib.Load(fPDBFileName, isProtein, verbosity);
    
    if (fMoleculeList != NULL) {
        fPDBlib.ComputeNbResiduesPerChain(fMoleculeList);
        fBarycenterList = fPDBlib.ComputeResidueBarycenters(fMoleculeList);
    }
    
    if (fMoleculeList != NULL) {
        std::cout << "Read file " << fPDBFileName << std::endl;
    }
    
    // Default parameters
    fThresDistForDSB = 10;
    fThresEdepForSSB = 8.22 * eV;
    
    fNbOfAlgo = 1;
    
    if ( fPm->ParameterExists(GetFullParmName("MinimumDistanceForDSB")) )
        fThresDistForDSB = fPm->GetIntegerParameter(GetFullParmName("MinimumDistanceForDSB"));
    if ( fPm->ParameterExists(GetFullParmName("LowerEnergyForSamplingSSB")) )
        fThresEdepForSSB = fPm->GetDoubleParameter(GetFullParmName("LowerEnergyForSamplingSSB"), "Energy");
    
    // This is for variance reduction
    if ( fPm->ParameterExists(GetFullParmName("NumberOfSplit")) )
        fNbOfAlgo = fPm->GetIntegerParameter(GetFullParmName("NumberOfSplit"));
    
    fThresEdepForSSB /= eV;
    
    fNtuple->RegisterColumnI(&fEventID, "Event number");
    fNtuple->RegisterColumnI(&fHits,    "Number of hits");
    fNtuple->RegisterColumnI(&fChainNum, "Chain number");
    fNtuple->RegisterColumnI(&fResidueNum, "Residue number");
    fNtuple->RegisterColumnS(&fAtomType, "Atom type");
    
    SuppressStandardOutputHandling();
    
}


TsScorePDB4DNA::~TsScorePDB4DNA() {
}


G4bool TsScorePDB4DNA::ProcessHits(G4Step* aStep,G4TouchableHistory*)
{
    if (!fIsActive) {
        fSkippedWhileInactive++;
        return false;
    }
    
    G4double edep = aStep->GetTotalEnergyDeposit()/eV;
    
    if ( edep > 0 ) {
        G4StepPoint* preStep = aStep->GetPreStepPoint();
        G4ThreeVector pos = preStep->GetPosition();
        G4double x = pos.x()/nanometer;
        G4double y = pos.y()/nanometer;
        G4double z = pos.z()/nanometer;
        
        int chainNum = 0;
        int residueNum = 0;
        string atomType;
        unsigned short int hit = fPDBlib.ComputeMatchEdepProtein(fBarycenterList, fMoleculeList, x*10, y*10, z*10,
                                                                 chainNum, residueNum, atomType);
        
        if ( 1 == hit ) {
            fHits++;
            fChainNum = chainNum;
            fResidueNum = residueNum;
            fAtomType = atomType;
            fNtuple->Fill();
        }
        return true;
    }
    return false;
}


void TsScorePDB4DNA::UserHookForEndOfEvent() {
    fEventID = GetEventID();
    fHits = 0;
    fChainNum = 0;
    fResidueNum = 0;
    fAtomType = "";
}


// This class was taken from Geant4/examples/extended/medical/dna/pdb4dna

void TsScorePDB4DNA::ComputeStrandBreaks(G4int* sb, G4int cluster)
{
    // sb quantities
    //
    G4int ssb1=0;
    G4int ssb2=0;
    G4int dsb=0;
    
    // nucleotide id and energy deposit for each strand
    G4int nucl1;
    G4int nucl2;
    G4double edep1;
    G4double edep2;
    
    //Read strand1
    //
    while ( !fVEdepStrand1[cluster].empty() )
    {
        nucl1 = fVEdepStrand1[cluster].begin()->first;
        edep1 = fVEdepStrand1[cluster].begin()->second;
        fVEdepStrand1[cluster].erase( fVEdepStrand1[cluster].begin() );
        
        // SSB in strand1
        //
        if ( edep1 >= fThresEdepForSSB )
        {
            ssb1++;
        }
        
        // Look at strand2
        //
        if ( !fVEdepStrand2[cluster].empty() )
        {
            do
            {
                nucl2 = fVEdepStrand2[cluster].begin()->first;
                edep2 = fVEdepStrand2[cluster].begin()->second;
                if ( edep2 >= fThresEdepForSSB )
                {
                    ssb2++;
                }
                fVEdepStrand2[cluster].erase( fVEdepStrand2[cluster].begin() );
            } while ( ((nucl1-nucl2)>fThresDistForDSB) && (!fVEdepStrand2[cluster].empty()) );
            
            // no dsb
            //
            if ( nucl2-nucl1 > fThresDistForDSB )
            {
                fVEdepStrand2[cluster][nucl2]=edep2;
                if ( edep2 >= fThresEdepForSSB )
                {
                    ssb2--;
                }
            }
            
            // one dsb
            //
            if ( std::abs(nucl2-nucl1) <= fThresDistForDSB )
            {
                if ( ( edep2 >= fThresEdepForSSB ) &&
                    ( edep1 >= fThresEdepForSSB ) )
                {
                    ssb1--;
                    ssb2--;
                    dsb++;
                }
            }
        }
    }
    
    // End with not processed data
    //
    while ( !fVEdepStrand1[cluster].empty() )
    {
        nucl1 = fVEdepStrand1[cluster].begin()->first;
        edep1 = fVEdepStrand1[cluster].begin()->second;
        if ( edep1 >= fThresEdepForSSB )
        {
            ssb1++;
        }
        fVEdepStrand1[cluster].erase( fVEdepStrand1[cluster].begin() );
    }
    
    while ( !fVEdepStrand2[cluster].empty() )
    {
        nucl2 = fVEdepStrand2[cluster].begin()->first;
        edep2 = fVEdepStrand2[cluster].begin()->second;
        if ( edep2 >= fThresEdepForSSB )
        {
            ssb2++;
        }
        fVEdepStrand2[cluster].erase( fVEdepStrand2[cluster].begin() );
    }
    sb[0]=ssb1+ssb2;
    sb[1]=dsb;
}
