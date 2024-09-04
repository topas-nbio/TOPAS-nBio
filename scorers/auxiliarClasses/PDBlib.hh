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
// This example is provided by the Geant4-DNA collaboration
// Any report or published results obtained using the Geant4-DNA software
// shall cite the following Geant4-DNA collaboration publication:
// Med. Phys. 37 (2010) 4692-4708
// Delage et al. PDB4DNA: implementation of DNA geometry from the Protein Data
//                  Bank (PDB) description for Geant4-DNA Monte-Carlo
//                  simulations (submitted to Comput. Phys. Commun.)
// The Geant4-DNA web site is available at http://geant4-dna.org
//
// --------------------------------------------------------------
// Authors: E. Delage
// november 2013
// --------------------------------------------------------------
//
// $Id$
//
/// \file PDBlib.hh
/// \brief Definition of the PDBlib class

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

#ifndef PDBlib_h
#define PDBlib_h 1

#include "PDBbarycenter.hh"
#include "PDBmolecule.hh"
#include <vector>

//! PDBlib Class
/*!
 * This Class define Molecule model ... 
 */
class PDBlib
{
public:
  //! First constructor
  PDBlib();

  //! Load PDB file into memory
  Molecule * Load(const string &filename, unsigned short int &isProtein,
                  unsigned short int verbose=0);

  //! Compute nucleotide barycenter from memory
  Barycenter * ComputeResidueBarycenters(Molecule * moleculeListTemp);

  //! Compute the corresponding bounding volume parameters
  void ComputeBoundingVolumeParams(Molecule *moleculeListTemp,
      double &dX,double &dY,double &dZ,       //Dimensions for bounding volume
      double &tX,double &tY,double &tZ);      //Translation for bounding volume

  //! Compute number of nucleotide per strand
  void ComputeNbNucleotidsPerStrand(Molecule * moleculeListTemp);

  //! Compute if energy is deposited in per atom
  unsigned short int ComputeMatchEdepDNA(Barycenter *,Molecule *,
      double x, double y,double z,
      int &numStrand, int &numNucleotid, int &codeResidue);

  //! Compute if energy is deposited in per atom
  unsigned short int ComputeMatchEdepProtein(Barycenter *BarycenterList,
      Molecule *moleculeListTemp,
      double x, double y, double z,
      int &chainNum, int &residueNum, string &atomType);

  //! return distance between two 3D points
  double DistanceTwo3Dpoints(double xA,double xB,
      double yA,double yB,
      double zA,double zB);

private:
  //! Number of nucleotid per strand
  int fNbNucleotidsPerStrand;

  //! Number of residues per chain
  int fNbResiduesPerChain;
};

#endif

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
