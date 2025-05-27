# TOPAS-CellModels
1) Description:

Cell Models for TOPAS/Geant4 and the inclusion of nano particles in particle scattering simulations.

The C++ classes in this repository extend the functionality of the TOPAS (https://github.com/OpenTOPAS/OpenTOPAS) Monte-Carlo program, which is itself a wrapper of the Geant4 MCS Toolkit (http://geant4.org).


2) Installation:

The following instructions assume you are on MacOS.
Change directories as appropriate for other operating systems:
- (see quickStart guides on https://github.com/OpenTOPAS/OpenTOPAS and the README on https://github.com/topas-nbio/TOPAS-nBio)

Navigate to the Applications directory:

  cd /Applications/

Clone or download the sourcecode into your TOPAS extension directory:
 
  git clone https://github.com/BAMresearch/TOPAS-CellModels.git
 
Change to your build directory:

  cd /Applications/TOPAS/OpenTOPAS-build/

Run cmake:

  cmake ../OpenTOPAS -DTOPAS_EXTENSIONS_DIR=/Applications/TOPAS-CellModels/

Followed by:

  make -j20 install

NOTE: In order to run other TOPAS-nBio simulations (i.e. NOT nanoparticles) you would need to redo the cmake/make described in the README (https://github.com/topas-nbio/TOPAS-nBio) 

3) Description:

A simple spherical cell with nanoparticles can be generated in a fast manner.
The user has the option to include nanoparticles and different organelles in the cell, e.g. a nucleus, mitochondria, a cell membrane.
Details can be found in https://doi.org/10.1038/s41598-021-85964-2

4) Usage:

Examples can be found in the  "examples/" directory.


 
5) Bugs:

Please report bugs to hahn@physik.fu-berlin.de or on https://github.com/BAMresearch/TOPAS-CellModels


6) Literature:

If you use this extension please cite the following literature:

Hahn, M.B., Zutta Villate, J.M. "Combined cell and nanoparticle models for TOPAS to study radiation dose enhancement in cell organelles."
Sci Rep 11, 6721 (2021). 

https://doi.org/10.1038/s41598-021-85964-2


7) Etc:

Tags: topas topasmc topasmcs topas-mc topas-mcs topas-nbio mcs 
