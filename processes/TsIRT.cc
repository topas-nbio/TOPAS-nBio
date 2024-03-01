// Extra Class for TsEmDNAChemistry
#include "TsIRT.hh"
#include "TsIRTUtils.hh"
#include "TsIRTConfiguration.hh"
#include "TsParameterManager.hh"

#include "G4Step.hh"
#include "G4Track.hh"
#include "G4Molecule.hh"
#include "G4VSolid.hh"

#include "globals.hh"
#include "G4SystemOfUnits.hh"
#include "G4ThreeVector.hh"
#include "Randomize.hh"
#include "G4Timer.hh"
#include "G4UIcommand.hh"
#include "G4UnitsTable.hh"
#include "G4LogicalVolumeStore.hh"
#include "G4GeometryTolerance.hh"

#include <stdlib.h>
#include <algorithm>
#include <string>
#include <iostream>
#include <set>

TsIRT::TsIRT(TsParameterManager* pM, G4String parmName): TsVIRTProcedure(pM,parmName){;}
TsIRT::~TsIRT() {;}