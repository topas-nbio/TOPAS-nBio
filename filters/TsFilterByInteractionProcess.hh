//
// ********************************************************************
// *                                                                  *
// * This file is part of the TOPAS-nBio extensions to the            *
// *   TOPAS Simulation Toolkit.                                      *
// * The TOPAS-nBio extensions are freely available under the license *
// *   agreement set forth at: https://topas-nbio.readthedocs.io/     *
// *                                                                  *
// ********************************************************************


#ifndef TsFilterByInteractionProcess_hh
#define TsFilterByInteractionProcess_hh

#include "TsVFilter.hh"

class G4ProcessVector;

class TsFilterByInteractionProcess : public TsVFilter
{
public:
	TsFilterByInteractionProcess(G4String name, TsParameterManager* pM, TsMaterialManager* mM, TsGeometryManager* gM,
							 TsFilterManager* fM, TsVGenerator* generator, TsVScorer* scorer, TsVFilter* parentFilter);
	virtual ~TsFilterByInteractionProcess();

	void ResolveParameters();

	virtual G4bool Accept(const G4Step*) const;
	virtual G4bool AcceptTrack(const G4Track*) const;

private:
	G4ProcessVector* fProcesses;
};
#endif
