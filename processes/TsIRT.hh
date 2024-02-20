#ifndef TsIRT_hh
#define TsIRT_hh

#include "TsIRTConfiguration.hh"
#include "TsVIRTProcedure.hh"

#include "globals.hh"
#include "G4ThreeVector.hh"
#include "G4Step.hh"
#include "G4VSolid.hh"

#include <vector>
#include <map>
#include <unordered_map>

class TsParameterManager;
class TsIRTUtils;

class TsIRT : public TsVIRTProcedure  {
public:
	TsIRT(TsParameterManager* pM, G4String parmName);
	~TsIRT();

	void SampleIndependantReactionTimes();
};
#endif

