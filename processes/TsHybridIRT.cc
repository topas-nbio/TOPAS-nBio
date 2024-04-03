// Extra Class for TsEmDNAChemistry
#include "TsHybridIRT.hh"
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

TsHybridIRT::TsHybridIRT(TsParameterManager* pM, G4String parmName)
: TsVIRTProcedure(pM,parmName), fPm(pM), fName(parmName), fHighTimeScavenger(false), fVerbosity(0), fScorersInitialized(false)
{
	fReactionConf = new TsIRTConfiguration(parmName, fPm);
	fUtils = new TsIRTUtils();
	
	G4String chemistryList = fPm->GetStringParameter("Ch/ChemistryName");
	
	if ( fPm->ParameterExists("Ts/Verbosity"))
		fVerbosity = fPm->GetIntegerParameter("Ts/Verbosity");
	
	if ( fPm->ParameterExists(GetFullParmName("TimeLower")) ) {
		G4double tmin = fPm->GetDoubleParameter(GetFullParmName("TimeLower"), "Time");
		G4double tmax = fPm->GetDoubleParameter(GetFullParmName("TimeUpper"), "Time");
		G4int    tbin = fPm->GetIntegerParameter(GetFullParmName("TimeBins"));
		fStepTimes = fUtils->CreateTimeSteps(tmin, tmax, tbin, true);
		fReactionConf->SetTimeLimits(tmin, tmax);
		if ( fPm->ParameterExists(GetFullParmName("EnableHighTimeScavengers")))
			fHighTimeScavenger = fPm->GetBooleanParameter(GetFullParmName("EnableHighTimeScavengers"));
		fTimeCut = tmax; 
	} else {
		fStepTimes = fUtils->CreateTimeSteps(1.0*ps, 1.0e6*ps, 100, true);
		fReactionConf->SetTimeLimits(1.0*ps, 1.0e6*ps);
		fTimeCut = 1.0e6*ps;
	}
	
	fUseSpinScaled = false;
	if (fPm->ParameterExists(GetFullParmName("UseScaledProbabilityForSpinBehavior")))
		fUseSpinScaled = fPm->GetBooleanParameter(GetFullParmName("UseScaledProbabilityForSpinBehavior"));
	
	G4cout << "\n -------------- TOPAS IRT Warning -------------" << G4endl;
	
	if ( fUseSpinScaled )
		G4cout << " -- SpinBehavior is ScaledProbability" << G4endl;
	else
		G4cout << " -- SpinBehavior is ExplicitAssignament" << G4endl;
	
	G4cout << " ----------------------------------------------" << G4endl;
	G4cout << G4endl;
	
	
	fMoleculesIDs = fReactionConf->GetMoleculeIDs();
	std::map<G4int, G4String> moleculesNames = fReactionConf->GetMoleculeNames();
	
	if ( !fPm->ParameterExists(GetFullParmName("ReportMoleculesNamed")) ) {
		for ( int i = 1; i <= (int)moleculesNames.size(); i++ )
			fMolecules[i] = moleculesNames[i];
		
	} else {
		G4String* reportMolecules = fPm->GetStringVector(GetFullParmName("ReportMoleculesNamed"));
		G4int nbOfReportedMolecules = fPm->GetVectorLength(GetFullParmName("ReportMoleculesNamed"));
		
		for ( int i = 0; i < nbOfReportedMolecules; i++ )
			fMolecules[fMoleculesIDs[reportMolecules[i]]] = reportMolecules[i];
	}
	
	fReportDelta = false;
	if ( fPm->ParameterExists(GetFullParmName("ReportDeltaGValues"))) {
		fReportDelta = fPm->GetBooleanParameter(GetFullParmName("ReportDeltaGValues"));
	}
	
	fRCutOff = fReactionConf->GetRCutOff(1*us);
	fBinWidth = fRCutOff;
	if ( fVerbosity > 0 )
		G4cout << " -- Cut off -- " << fRCutOff/nm << " nm " << std::endl;
	
	if ( fPm->ParameterExists("Ch/"+chemistryList+"/SpaceBinningWidth") )
		fBinWidth = fPm->GetDoubleParameter("Ch/"+chemistryList+"/SpaceBinningWidth", "Length");
	
	fTestForContactReactions = false;
	if ( fPm->ParameterExists("Ch/"+chemistryList+"/TestForContactReactions"))
		fTestForContactReactions = fPm->GetBooleanParameter("Ch/"+chemistryList+"/TestForContactReactions");
	
	if ( fVerbosity > 0 && fTestForContactReactions )
		std::cout << " Going to test tracks by reaction at contact. " << G4endl;
	
	fXMin = 1e9*nm;
	fYMin = 1e9*nm;
	fZMin = 1e9*nm;
	
	fXMax = 0e0*nm;
	fYMax = 0e0*nm;
	fZMax = 0e0*nm;

	fChunk   = 0;
	fSpeciesIndex = 0;
	fIndex   = 0;
	fSampleIRTatStart = true;
	fMinTime = 0;
	fCubicVolume = 0;

	fGillespieFinalTau = true;
	if (fPm->ParameterExists("Ch/"+chemistryList+"/EndChemistryWithTau")) 
		fGillespieFinalTau = fPm->GetBooleanParameter("Ch/"+chemistryList+"/EndChemistryWithTau");
}


TsHybridIRT::~TsHybridIRT() {
	delete fReactionConf;
	delete fUtils;
}


G4bool TsHybridIRT::Inside(G4ThreeVector p ) {
	G4double delta = 1 * nm;
	G4double dist = std::max(std::max(
									  std::abs(p.x())-fDx,
									  std::abs(p.y())-fDy),
							 std::abs(p.z())-fDz);
	if (dist > delta) return false;
	return !(dist > -delta);
}


void TsHybridIRT::AddMolecule(TsIRTConfiguration::TsMolecule aMol) {
	if ( !fScorersInitialized )
		initializeScorers();

	G4ThreeVector position = aMol.position;
	
	fSpeciesOfAKind[aMol.id][fSpeciesIndex] = true;
	fChemicalSpecies[fSpeciesIndex] = aMol;
	fSpeciesIndex++;
	
	// count
	G4int tBin = fUtils->FindBin(aMol.time, fStepTimes);
	
	if ( -1 < tBin ) {
		for ( int tbin = tBin; tbin < (int)fStepTimes.size(); tbin++ ) {
			if ( fTheGvalue.find(aMol.id) == fTheGvalue.end() )
				fTheGvalue[aMol.id][tbin] = 1;
			else
				fTheGvalue[aMol.id][tbin]++;
		}
	}
	
	if ( fXMin > position.x() ) fXMin = position.x();
	if ( fYMin > position.y() ) fYMin = position.y();
	if ( fZMin > position.z() ) fZMin = position.z();
	
	if ( fXMax < position.x() ) fXMax = position.x();
	if ( fYMax < position.y() ) fYMax = position.y();
	if ( fZMax < position.z() ) fZMax = position.z();
}


void TsHybridIRT::AddMolecule(G4Step* aStep, G4double time, G4int moleculeID, G4ThreeVector offset) {
	G4int pdg = -1;
	G4ThreeVector position = aStep->GetPreStepPoint()->GetPosition();
	const G4String& name = GetMolecule(aStep->GetTrack())->GetName();
	
	pdg = fMoleculesIDs[name];
	
	if ( pdg > 0 ) {
		TsIRTConfiguration::TsMolecule aMol;
		aMol.id = pdg;
		aMol.position = position + offset;
		aMol.time = time;
		aMol.reacted = false;
		// trackID != 0 reserved for DNA, < 0 = left strand, > 0 = right strand
		aMol.trackID = moleculeID;
		aMol.isDNA = false;
		aMol.isNew = true;
		
		if ( pdg == 1 | pdg == 5 )
			aMol.spin = G4UniformRand() > 0.5 ? 1 : 0;
		else
			aMol.spin = -1;
		
		AddMolecule(aMol);
	}
}


void TsHybridIRT::AddMolecule(G4Track* aTrack, G4double time, G4int moleculeID, G4ThreeVector offset, G4bool isDNA) {
	G4int pdg = -1;
	G4ThreeVector position = aTrack->GetPosition();
	const G4String& name = GetMolecule(aTrack)->GetName();
	
	pdg = fMoleculesIDs[name];
	      
	if ( pdg > 0 ) {
		TsIRTConfiguration::TsMolecule aMol;
		aMol.id = pdg;
		aMol.position = position + offset;
		aMol.time = time;
		aMol.reacted = false;
		aMol.trackID = moleculeID;
		aMol.isDNA   = isDNA;
		aMol.isNew   = true;
		
		if ( pdg == 1 | pdg == 5 )
			aMol.spin = G4UniformRand() > 0.5 ? 1 : 0;
		else
			aMol.spin = -1;
		
		AddMolecule(aMol);
	}	
}


void TsHybridIRT::AddMolecule(G4Track* aTrack, G4double time, G4int moleculeID, G4ThreeVector offset) {
	G4int pdg = -1;
	G4ThreeVector position = aTrack->GetPosition();
	const G4String& name = GetMolecule(aTrack)->GetName();
	
	pdg = fMoleculesIDs[name];
	
	if ( pdg > 0 ) {
		TsIRTConfiguration::TsMolecule aMol;
		aMol.id = pdg;
		aMol.position = position + offset;
		aMol.time = time;
		aMol.reacted = false;
		aMol.isNew   = true;
		// trackID != 0 reserved for DNA, < 0 = left strand, > 0 = right strand
		aMol.trackID = moleculeID;

		if ( pdg == 1 | pdg == 5 )
			aMol.spin = G4UniformRand() > 0.5 ? 1 : 0;
		else
			aMol.spin = -1;
		
		AddMolecule(aMol);
	}
}


void TsHybridIRT::AddMolecule(const G4Track* aTrack, G4double time, G4int moleculeID, G4ThreeVector offset) {
	G4int pdg = -1;
	G4ThreeVector position = aTrack->GetPosition();
	const G4String& name = GetMolecule(aTrack)->GetName();
	
	pdg = fMoleculesIDs[name];
	
	if ( pdg > 0 ) {
		TsIRTConfiguration::TsMolecule aMol;
		aMol.id = pdg;
		aMol.position = position + offset;
		aMol.time = time;
		aMol.reacted = false;
		aMol.isNew   = true;
		// trackID != 0 reserved for DNA, < 0 = left strand, > 0 = right strand
		aMol.trackID = moleculeID;
		
		if ( pdg == 1 | pdg == 5 )
			aMol.spin = G4UniformRand() > 0.5 ? 1 : 0;
		else
			aMol.spin = -1;
		
		AddMolecule(aMol);
	}
}


TsIRTConfiguration::TsMolecule TsHybridIRT::ConstructMolecule(G4Track* aTrack, G4double time, G4int moleculeID, G4ThreeVector offset) {
	G4int pdg = -1;
	G4ThreeVector position = aTrack->GetPosition();
	const G4String& name = GetMolecule(aTrack)->GetName();
	pdg = fMoleculesIDs[name];

	TsIRTConfiguration::TsMolecule aMol;

	if ( pdg > 0 ) {	
		aMol.id = pdg;
		aMol.position = position + offset;
		aMol.time = time;
		aMol.reacted = false;
		// trackID != 0 reserved for DNA, < 0 = left strand, > 0 = right strand
		aMol.trackID = moleculeID;
		aMol.isDNA = false;
		aMol.isNew = true;

		if ( pdg == 1 | pdg == 5 )
			aMol.spin = G4UniformRand() > 0.5 ? 1 : 0;
		else
			aMol.spin = -1;
	}

	return aMol;
}


void TsHybridIRT::Clean() {
	fChemicalSpecies.clear();
	fTheGvalue.clear();
	fTheGvalueInVolume.clear();
	fSpaceBinned.clear();
	fSpeciesOfAKind.clear();
	fUsed.clear();

	fHomogeneousConcentrations.clear();
	fConcentrationsAtThisTime.clear();

	fGValues.clear();
	
	if (fReportDelta)
		fDeltaGValues.clear();
	
	fScorersInitialized = false;
	fSampleIRTatStart = true;
	fMinTime = 0;
	fTimeCut = fStepTimes[fStepTimes.size()-1];

	CleanIRTVariables();
}


void TsHybridIRT::FindBinIndexes(G4ThreeVector thisPos, G4double rcutOff) {
	fxiniIndex = fUtils->FindBin(fNx, fXMin, fXMax, thisPos.x()-rcutOff);
	fxendIndex = fUtils->FindBin(fNx, fXMin, fXMax, thisPos.x()+rcutOff);
	fyiniIndex = fUtils->FindBin(fNy, fYMin, fYMax, thisPos.y()-rcutOff);
	fyendIndex = fUtils->FindBin(fNy, fYMin, fYMax, thisPos.y()+rcutOff);
	fziniIndex = fUtils->FindBin(fNz, fZMin, fZMax, thisPos.z()-rcutOff);
	fzendIndex = fUtils->FindBin(fNz, fZMin, fZMax, thisPos.z()+rcutOff);
}


void TsHybridIRT::initializeScorers() {
	for ( auto& nameAndID : fMolecules ) {
		G4int id = nameAndID.first;
		for ( int i = 0; i < (int)fStepTimes.size(); i++ ) {
			fTheGvalue[id][i] = 0.0;
		}
	}
	
	if ( fReportDelta ) {
		for  (int id = 0; id < fReactionConf->GetNumberOfReactions(); id++ ) {
			for (int i = 0; i < (int)fStepTimes.size(); i++ ) {
				G4double time = fStepTimes[i];
				fDeltaGValues[id][time] = 0.0;
			}
		}
	}
	
	fScorersInitialized = true;
}


void TsHybridIRT::contactReactions(G4int i,std::unordered_map<G4int, G4bool> used) {
	fReactedByContact = false;
	fContactProducts.clear();

	if (!MoleculeExists(i)) {return;}

	// background reaction
	G4int indexOf1stOrder = fReactionConf->ContactFirstOrderAndBackgroundReactions(fChemicalSpecies[i]);
	if ( -1 < indexOf1stOrder ) {
		fReactedByContact = fReactionConf->MakeReaction(fChemicalSpecies, fSpeciesIndex, fSpaceBinned,
														fNx, fNy, fNz, fXMin, fXMax,
														fYMin, fYMax, fZMin, fZMax,
														fTheGvalue, fStepTimes, i,
														indexOf1stOrder, fChemicalSpecies[i].time, fUsed, fContactProducts, fSpeciesOfAKind);
		if ( fReactedByContact ) {
			if ( fReportDelta ) {
				G4int tBin = fUtils->FindBin(fChemicalSpecies[i].time, fStepTimes);
				if ( -1 < tBin ) {
					for ( int ti = tBin; ti < (int)fStepTimes.size(); ti++ ) {
						fDeltaGValues[indexOf1stOrder][fStepTimes[ti]]++;
					}
				}
			}
			fChemicalSpecies[i].reacted = true;
			RemoveMolecule(i);
			return;
		}
	}
	
	G4ThreeVector thisPos = fChemicalSpecies[i].position;
	FindBinIndexes(thisPos, 0.0);
	
	for ( int ii = fxiniIndex; ii <= fxendIndex; ii++ ) {
		for ( int jj = fyiniIndex; jj <= fyendIndex; jj++ ) {
			for ( int kk = fziniIndex; kk <= fzendIndex; kk++ ) {
				for (auto& IndexAndAlive:fSpaceBinned[ii][jj][kk]) {
					G4int j = IndexAndAlive.first;

					if (!MoleculeExists(j) || fReactedByContact) {continue;}
					if (!fChemicalSpecies[j].isNew) {continue;}
					if (used.count(j) >= 1) {continue;}
					
					if ( j == i )
						continue;
					
					G4int indexOfReaction = fReactionConf->GetReactionIndex(fChemicalSpecies[j].id,
																			fChemicalSpecies[i].id);
					if ( -1 < indexOfReaction ) {
						TsIRTConfiguration::TsMolecularReaction binReaction = fReactionConf->GetReaction(indexOfReaction);
						if (binReaction.index < 0) {continue;} 
						G4double r = binReaction.effectiveReactionRadius;
						if ( binReaction.reactionType == 4 )
							r = binReaction.effectiveTildeReactionRadius;
						
						G4double p = binReaction.probabilityOfReaction;
						G4double timeJ = fChemicalSpecies[j].time;
						
						if ( fChemicalSpecies[i].time != timeJ ) continue; /// Only zero timers<---Best match Clifford et al.
						
						if ( (fChemicalSpecies[i].position - fChemicalSpecies[j].position).mag2() < r*r ) {
							fReactedByContact = fReactionConf->MakeReaction(fChemicalSpecies, fSpeciesIndex, fSpaceBinned,
																			fNx, fNy, fNz, fXMin, fXMax,
																			fYMin, fYMax, fZMin, fZMax,
																			fTheGvalue, fStepTimes, i,
																			j, indexOfReaction, fChemicalSpecies[i].time, p, fUsed, fContactProducts, fSpeciesOfAKind);
							
							if ( fReactedByContact ) {
								if ( fReportDelta ) {
									G4int tdBin = fUtils->FindBin(fChemicalSpecies[i].time, fStepTimes);
									if ( -1 < tdBin ) {
										for ( int ti = tdBin; ti < (int)fStepTimes.size(); ti++ ) {
											fDeltaGValues[indexOfReaction][fStepTimes[ti]]++;
										}
									}
								}
								fChemicalSpecies[i].reacted = true;
								fChemicalSpecies[j].reacted = true;
								RemoveMolecule(i);
								RemoveMolecule(j);
								return;
							}
						}
					}
				}
			}
		}
	}
}


void TsHybridIRT::sampleReactions(G4int i) {
	if ( !MoleculeExists(i) )
		return;
	
	if ( fChemicalSpecies[i].isDNA )
		return;
	
	if ( !fChemicalSpecies[i].isNew )
		return;

	G4ThreeVector thisPos = fChemicalSpecies[i].position;
	FindBinIndexes(thisPos, fRCutOff);

	// Type I - V sampling
	fUsed[i] = true;
	for ( int ii = fxiniIndex; ii <= fxendIndex; ii++ ) {
		for ( int jj = fyiniIndex; jj <= fyendIndex; jj++ ) {
			for ( int kk = fziniIndex; kk <= fzendIndex; kk++ ) {
				for (auto& IndexAndAlive:fSpaceBinned[ii][jj][kk]) {
					G4int j = IndexAndAlive.first;
					
					if ( !MoleculeExists(j) || (fUsed[j] && fSampleIRTatStart))
						continue;

					if ( !fChemicalSpecies[j].isNew )
						return;

					if (j == i) continue;

					G4int indexOfReaction = fReactionConf->GetReactionIndex(fChemicalSpecies[i].id,
																			fChemicalSpecies[j].id);
					if( -1 < indexOfReaction ) {
						G4double timeI = fChemicalSpecies[i].time;
						G4double timeJ = fChemicalSpecies[j].time;
						G4double dt = fChemicalSpecies[i].time - timeJ;
						
						G4ThreeVector origPositionI = fChemicalSpecies[i].position;
						G4ThreeVector origPositionJ = fChemicalSpecies[j].position;
						G4bool resampleA = false, resampleB = false;
						if ( 0 < dt ) {
							fReactionConf->Diffuse(fChemicalSpecies[j],dt);
							resampleB = true;
						}
						
						if ( 0 > dt ) {
							fReactionConf->Diffuse(fChemicalSpecies[i],-dt);
							resampleA = true;
						}
						
						G4double irt = fReactionConf->GetIndependentReactionTime(fChemicalSpecies[i],fChemicalSpecies[j],indexOfReaction);
						G4double gTime = fChemicalSpecies[i].time;

						if ( (0 < irt) && (irt + gTime <= fTimeCut) && (fMinTime <= irt + gTime) && (irt + gTime <= fTransCut)) {
							irt += gTime;
							AddToIRT(irt,indexOfReaction,i,j,fChemicalSpecies[i].time,fChemicalSpecies[i].position,fChemicalSpecies[j].position,false);
				
						}

						G4double delta_pos = (fChemicalSpecies[i].position - fChemicalSpecies[j].position).mag();

						if ( resampleB ) {
							fChemicalSpecies[j].position = origPositionJ;
							fChemicalSpecies[j].time = timeJ;
						}
						
						if ( resampleA ) {
							fChemicalSpecies[i].position = origPositionI;
							fChemicalSpecies[i].time = timeI;
						}
					}
				}
			}
		}
	}
}


void TsHybridIRT::runIRT(G4double startTime, G4double finalTime, G4double transTime, G4bool isContinuation) {
	G4String component = fPm->GetStringParameter(GetFullParmName("Component"));
	fCubicVolume = 1000 * G4LogicalVolumeStore::GetInstance()->GetVolume(component,false)->GetSolid()->GetCubicVolume() / (m3);

	fMinTime = 0;
	fTimeCut = fStepTimes[fStepTimes.size()-1];
	if (startTime != -1)
		fMinTime = startTime;
	if (finalTime != -1)
		fTimeCut = finalTime;
	if (transTime != -1)
		fTransCut = transTime;
	else
		fTransCut = DBL_MAX;

	fInitialBackgroundConcentrations = fReactionConf->GetBackgroundConcentrations();
	if (!isContinuation) {fCurrentBackgroundConcentrations=fInitialBackgroundConcentrations;}
		
	// Voxelize and Sort the Chemical Species Space
	VoxelizeAndSortSpace();

	// Test For initial Contact Reactions; Needed for QA
	TestForContactReactions();

	// Initial Sample of the IRT
	SampleIndependantReactionTimes();

	// Make the reactions according to the IRT
	ConductReactions();

	// Update the GValue containers for output
	UpdateGValues();
}


void TsHybridIRT::VoxelizeAndSortSpace() {
	fNx = G4int((fXMax-fXMin)/fBinWidth) == 0 ? 1 : G4int((fXMax-fXMin)/fBinWidth);
	fNy = G4int((fYMax-fYMin)/fBinWidth) == 0 ? 1 : G4int((fYMax-fYMin)/fBinWidth);
	fNz = G4int((fZMax-fZMin)/fBinWidth) == 0 ? 1 : G4int((fZMax-fZMin)/fBinWidth);
	
	if ( 1 < fVerbosity )
		std::cout << " Going to bin the space containing the track into "
		<< fNx << " x, " << fNy << " y and " << fNz << " z bins of width "
		<< fBinWidth/nm << " nm with limits (" << fXMin/nm << "," << fXMax/nm << "), ("
		<< fYMin/nm << ", " << fYMax/nm << "), (" << fZMin/nm << ", " << fZMax/nm << ") nm" << std::endl;

	for (auto& IndexAndMol:fChemicalSpecies) {
		G4int t   = IndexAndMol.first;
		fUsed[t] = false;
		G4int I = fUtils->FindBin(fNx, fXMin, fXMax, fChemicalSpecies[t].position.x());
		G4int J = fUtils->FindBin(fNy, fYMin, fYMax, fChemicalSpecies[t].position.y());
		G4int K = fUtils->FindBin(fNz, fZMin, fZMax, fChemicalSpecies[t].position.z());
		fSpaceBinned[I][J][K][t] = true;
	}

	if ( 1 < fVerbosity )
		std::cout << " Processing " << fChemicalSpecies.size() << " species " << std::endl;
}

void TsHybridIRT::TestForContactReactions() {
	if (!fTestForContactReactions) { return; }

	std::unordered_map<G4int, TsIRTConfiguration::TsMolecule> MockMap = fChemicalSpecies;
	std::unordered_map<G4int,G4bool> AlreadyContactSampled;

	for (auto& IndexAndMol:MockMap) {
		G4int i = IndexAndMol.first;
		if ( fChemicalSpecies[i].isDNA ) { continue; }
		if ( !fChemicalSpecies[i].isNew) { continue; }
		contactReactions(i,AlreadyContactSampled);
		AlreadyContactSampled[i] = true;
	}
	
	if ( 1 < fVerbosity )
		std::cout << "     -- After processing remained "
		<< fChemicalSpecies.size() << " species " << std::endl;
}

void TsHybridIRT::SampleIndependantReactionTimes() {
	fIRTIndex.clear();
	fIRTValues.clear();
	fIRTMolA.clear();
	fIRTMolB.clear();
	fIRTIsBackground.clear();
	fIRTOrigTime.clear();
	fIRTPositions.clear();

	if ( 2 < fVerbosity ) {
		fChunk = (G4int)fChemicalSpecies.size()/10;
		G4cout << "     -- Begin sampling of IRT for " << fChemicalSpecies.size() << " species " << G4endl;
	}
	
	for (auto& IndexAndMol:fChemicalSpecies) {
		G4int i = IndexAndMol.first;

		sampleReactions(i);
		
		if ( 2 < fVerbosity && ( i >= fChunk && i % fChunk == 0 )) {
			G4cout << "     ---- current sampled times " << i << " out " << (G4int)fChemicalSpecies.size() <<  G4endl;
		}
	}

	if ( 2 < fVerbosity ) {
		fChunk = (int)fIRTValues.size()/10;
		G4cout << "     -- Begin performing of reactions " << G4endl;
	}
}

void TsHybridIRT::SortIndependantReactionTimes() {
	if (!fSortByTime)
		std::sort(fIRTValues.begin(), fIRTValues.end());
}


void TsHybridIRT::ConductReactions() {
	fSampleIRTatStart = false;
	G4double currentTime = 1*ps;
	if (fMinTime != 0) 
		currentTime = fMinTime;

	G4int Count = 0;
	while (fIRTValues.size() > 0) {
		G4int currentIRT = 0;
		Count++;
		if (!IsNextIRTStillPossible()) {
			RemoveFirstIRTElement();
			continue;
		}

		if (currentTime >= fTransCut) { break; }

		G4int dnaVolumeID = -1;
		G4int dnaBaseID   = -1;
		G4int dnaStrandID = -1;

		G4double irt = fIRTValues[currentIRT].first;

		UpdateGillespieDirect(currentTime);
		if ((fHomogeneousTimeStep + currentTime < irt)      &&
			(fHomogeneousTimeStep + currentTime < fTimeCut) &&
			(fHomogeneousReactIndex >= 0)) {
			G4bool Success = DoGillespieDirect(fHomogeneousTimeStep+currentTime);
			if (Success) {
				currentTime  += fHomogeneousTimeStep;
				continue;
			}
		}

		if ( irt < fTimeCut ) {
			G4int idx = fIRTValues[currentIRT].second;
			G4int iM = fIRTMolA[idx];
			G4int jM = fIRTMolB[idx];

			currentTime = irt;
			G4int indexOfReaction = fIRTIndex[idx];
			
			// 5.1 set the position of the reactants to those used to sample the irt
			TsIRTConfiguration::TsMolecularReaction binReaction = fReactionConf->GetReaction(indexOfReaction);
			if (binReaction.index < 0) {continue;}
			if ( binReaction.reactionType == 5 ) {
				if (fUseSpinScaled) {
					if (fChemicalSpecies[iM].spin == fChemicalSpecies[jM].spin ) {
						RemoveFirstIRTElement();
						continue;
					}
					else if ( 0.5 > G4UniformRand() ) {
						RemoveFirstIRTElement();
						continue;
					}
				}
				else {
					if ( !(0 == fChemicalSpecies[iM].spin  && fChemicalSpecies[iM].spin == fChemicalSpecies[jM].spin  )) {
						RemoveFirstIRTElement();
						continue;
					}
				}
			}
			
			G4int tBin = fUtils->FindBin(irt, fStepTimes);
			
			fChemicalSpecies[iM].position = fIRTPositions[idx].first;
			fChemicalSpecies[iM].time = fIRTOrigTime[idx];
			
			std::vector<G4ThreeVector> positions;
			std::vector<G4int> products;
			G4bool reactionIsInside = false;
			
			// 5.2 sample the position of the reactants using position approach Clifford et al 1986
			if ( !fIRTIsBackground[idx] ) {
				fChemicalSpecies[jM].position = fIRTPositions[idx].second;
				fChemicalSpecies[jM].time = fIRTOrigTime[idx];
				
				fReactionConf->ResampleReactantsPosition(fChemicalSpecies[iM], fChemicalSpecies[jM],indexOfReaction, irt);
				
				positions = fReactionConf->GetPositionOfProducts(fChemicalSpecies[iM],fChemicalSpecies[jM], indexOfReaction);
				
				TsIRTConfiguration::TsMolecularReaction binReaction = fReactionConf->GetReaction(indexOfReaction);
				if (binReaction.index < 0) {continue;}
				products = binReaction.products;
				
				if ( 0 <= tBin ) {
					for ( int ti = tBin; ti < (int)fStepTimes.size(); ti++ ) {
						fTheGvalue[fChemicalSpecies[iM].id][ti]--;
						fTheGvalue[fChemicalSpecies[jM].id][ti]--;
						if ( fReportDelta )
							fDeltaGValues[indexOfReaction][fStepTimes[ti]]++;
						
					}
				}
				
				if (fChemicalSpecies[iM].volumeID >= 0) {
					dnaVolumeID = fChemicalSpecies[iM].volumeID;
					dnaBaseID   = fChemicalSpecies[iM].baseID;
					dnaStrandID = fChemicalSpecies[iM].strandID;
				}
				else if (fChemicalSpecies[jM].volumeID >= 0) {
					dnaVolumeID = fChemicalSpecies[jM].volumeID;
					dnaBaseID   = fChemicalSpecies[jM].baseID;
					dnaStrandID = fChemicalSpecies[jM].strandID;
				}

				fChemicalSpecies[iM].reacted = true;
				fChemicalSpecies[jM].reacted = true;
				RemoveMolecule(iM);
				RemoveMolecule(jM);
				
			} else {
				fReactionConf->Diffuse(fChemicalSpecies[iM],irt-fChemicalSpecies[iM].time);
				positions = fReactionConf->GetBackgroundPositionOfProducts(fChemicalSpecies[iM],indexOfReaction);
				
				TsIRTConfiguration::TsMolecularReaction binReaction = fReactionConf->GetReaction(indexOfReaction);
				if (binReaction.index < 0) {continue;}
				products = binReaction.products;
				if ( 0 <= tBin ) {
					for ( int ti = tBin; ti < (int)fStepTimes.size(); ti++ ) {
						fTheGvalue[fChemicalSpecies[iM].id][ti]--;
						if ( fReportDelta )
							fDeltaGValues[indexOfReaction][fStepTimes[ti]]++;
					}
				}

				if (fChemicalSpecies[iM].volumeID >= 0) {
					dnaVolumeID = fChemicalSpecies[iM].volumeID;
					dnaBaseID   = fChemicalSpecies[iM].baseID;
					dnaStrandID = fChemicalSpecies[iM].strandID;
				}
				
				fChemicalSpecies[iM].reacted = true;
				RemoveMolecule(iM);
			}
			
			// 5.3 create new products at positions resampled in previous step
			for ( size_t u = 0; u < products.size(); u++ ) {
				TsIRTConfiguration::TsMolecule aProd;
				aProd.id = products[u];
				aProd.position = positions[u];
				aProd.time = irt;
				aProd.reacted = false;
				aProd.trackID = 0;
				aProd.isDNA = false;
				aProd.isNew = true;

				if (fMolecules[aProd.id] == "") {continue;}

				if (dnaVolumeID >= 0) {
					aProd.volumeID = dnaVolumeID;
					aProd.baseID   = dnaBaseID;
					aProd.strandID = dnaStrandID;
				}
				
				if ( 1 == products[u] || 5 == products[u] )
					aProd.spin = G4UniformRand() > 0.5 ? 1 : 0;
				else
					aProd.spin = -1;
				
				G4int newID = fSpeciesIndex;
				
				G4int I = fUtils->FindBin(fNx, fXMin, fXMax, positions[u].x());
				G4int J = fUtils->FindBin(fNy, fYMin, fYMax, positions[u].y());
				G4int K = fUtils->FindBin(fNz, fZMin, fZMax, positions[u].z());
				
				// Merge this specie with the track
				fSpaceBinned[I][J][K][fSpeciesIndex] = true;
				fChemicalSpecies[fSpeciesIndex] = aProd;
				fUsed[fSpeciesIndex] = false;
				fSpeciesOfAKind[products[u]][fSpeciesIndex] = true;
				fSpeciesIndex++;
				
				if ( 0 <= tBin ) {
					for ( int ti = tBin; ti < (int)fStepTimes.size(); ti++ )
						fTheGvalue[aProd.id][ti]++;
				}
				
				std::unordered_map<G4int, G4bool> emptyUsed;
				contactReactions(newID,emptyUsed);
				if (!fReactedByContact) {
					// Sample the New Molecule if it didn't reacted by contact
					sampleReactions(newID);
				}
				else {
					// Sample the products of the Molecule
					for (size_t jProds = 0; jProds < fContactProducts.size(); jProds++) {
						sampleReactions(fContactProducts[jProds]);
					}
				}
			}
		} else {
			break;
		}
		
		if (fVerbosity > 2 && ( currentIRT >= fChunk && currentIRT % fChunk == 0 )) {
			G4cout << "     ---- current reactions performed " << currentIRT << " out " << fChunk*10 << G4endl;
		}
		RemoveFirstIRTElement();
	}

	if (fGillespieFinalTau) {
		// Convert Chemical Species to Homogeneous
		for (auto& IndexAndConcentrations:fSpeciesOfAKind) {
			G4int Index = IndexAndConcentrations.first;
			fHomogeneousConcentrations[Index] += fSpeciesOfAKind[Index].size() + fConcentrationsAtThisTime[Index].size();
		}
		fSpeciesOfAKind.clear();
		fConcentrationsAtThisTime.clear();

		G4int Tries = 0;
		while (currentTime < fTimeCut) {
			UpdateGillespieTauLeaping(currentTime);
			G4bool Success = DoGillespieTauLeaping(currentTime + fHomogeneousTimeStep);
			if (Success) {
				currentTime  += fHomogeneousTimeStep; 
				Tries = 0; 
			}
			else {Tries++;}
			if (Tries == 100) {break;}
		}
	}

	else {
		while (currentTime < fTimeCut) {
			UpdateGillespieDirect(currentTime);
			if ((fHomogeneousTimeStep + currentTime < fTimeCut) && (fHomogeneousReactIndex >= 0)) {
				G4bool Success = DoGillespieDirect(fHomogeneousTimeStep+currentTime);
				if (Success) {
					currentTime  += fHomogeneousTimeStep;
					continue;
				}
			}
		}
	}
}

void TsHybridIRT::UpdateGValues() {
	for ( auto& nameAndID : fMolecules ) {
		G4String name = nameAndID.second;
		G4int id = nameAndID.first;
		for ( auto& timeAndGvalue : fTheGvalue[id] ) {
			G4int tBin = timeAndGvalue.first;
			G4double time = fStepTimes[tBin];
			fGValues[name][time] = timeAndGvalue.second;
		}
	}
}

void TsHybridIRT::CleanIRTVariables() {
	fIRTIndex.clear();
	fIRTValues.clear();
	fIRTMolA.clear();
	fIRTMolB.clear();
	fIRTIsBackground.clear();
	fIRTOrigTime.clear();
	fIRTPositions.clear();
	fTheGvalue.clear();
	fTheGvalueInVolume.clear();
	fSpaceBinned.clear();
	fUsed.clear();
	
	fSampleIRTatStart = true;
	fXMin = 1e9*nm;
	fYMin = 1e9*nm;
	fZMin = 1e9*nm;
	
	fXMax = 0e0*nm;
	fYMax = 0e0*nm;
	fZMax = 0e0*nm;

	fTimeCut = 1E12*s;
	fSpeciesIndex = 0;
	fIndex = 0;
}


void TsHybridIRT::AddToIRT(G4double Time, G4int Reaction, G4int molA, G4int molB, G4double OrigTime, G4ThreeVector aPos, G4ThreeVector bPos, G4bool isBack) {
	AddIRTinAscendantOrder(Time);
	fIRTIndex[fIndex] = Reaction;
	fIRTMolA[fIndex]  = molA;
	fIRTMolB[fIndex]  = molB;
	fIRTOrigTime[fIndex]     = OrigTime;
	fIRTPositions[fIndex]    = std::make_pair(aPos, bPos);
	fIRTIsBackground[fIndex] = isBack;
	fIndex++;	
}


void TsHybridIRT::AddIRTinAscendantOrder(G4double Time) {
	size_t N = fIRTValues.size();
	std::pair<G4double,G4int> newPair = std::make_pair(Time,fIndex);
	if ( N == 0) {
		fIRTValues.push_back(newPair);
	}
	else if (Time < fIRTValues[0].first) {
		fIRTValues.insert(fIRTValues.begin(),newPair);
	}
	else if (Time >= fIRTValues[N-1].first) {
		fIRTValues.push_back(newPair);
	}
	else {
		std::vector<std::pair<G4double,G4int>>::iterator it = std::lower_bound(fIRTValues.begin(),fIRTValues.end(),newPair);
		fIRTValues.insert(it,newPair);
	}
}


void TsHybridIRT::RemoveFirstIRTElement() {
	G4int index = fIRTValues[0].second;
	fIRTValues.erase(fIRTValues.begin());
	fIRTIndex.erase(index);
	fIRTMolA.erase(index);
	fIRTMolB.erase(index);
	fIRTOrigTime.erase(index);
	fIRTPositions.erase(index);
	fIRTIsBackground.erase(index);
}


void TsHybridIRT::RemoveMolecule(G4int Index) {
	G4int id = fChemicalSpecies[Index].id;
	G4int I = fUtils->FindBin(fNx, fXMin, fXMax, fChemicalSpecies[Index].position.x());
	G4int J = fUtils->FindBin(fNy, fYMin, fYMax, fChemicalSpecies[Index].position.y());
	G4int K = fUtils->FindBin(fNz, fZMin, fZMax, fChemicalSpecies[Index].position.z());
	fSpaceBinned[I][J][K].erase(Index);
	fChemicalSpecies.erase(Index);
	fSpeciesOfAKind[id].erase(Index);
	fUsed.erase(Index);

	for (size_t i = 0; i < fConcentrationsAtThisTime[id].size(); i++) {
		if (fConcentrationsAtThisTime[id][i] == Index) {
			fConcentrationsAtThisTime[id].erase(fConcentrationsAtThisTime[id].begin()+i);
			break;
		}
	}
}


G4String TsHybridIRT::GetFullParmName(G4String suffix ) {
	return "Sc/" + fName + "/" + suffix;
}

void TsHybridIRT::SetGValues(std::map<G4String, std::map<G4double, G4int>> gval) {
	fGValues   = gval;
	for ( auto& nameAndID : fMolecules ) {
		G4String name = nameAndID.second;
		G4int id = nameAndID.first;
		G4int bin = 0;
		for ( auto& timeAndGvalue : gval[name] ) {
			fTheGvalue[id][bin] = timeAndGvalue.second;
			bin++;
		}
	}
}

std::vector<G4int> TsHybridIRT::GetSurvivingMolecules(G4int parentID) {
	std::vector<G4int> AliveMolecules;
	for (auto& IndexAndMol: fChemicalSpecies) {
		G4int i = IndexAndMol.first;
		if (!fChemicalSpecies[i].reacted && fChemicalSpecies[i].trackID != 0) {
			if (parentID == -1)
				AliveMolecules.push_back(fChemicalSpecies[i].trackID);
			else if (fChemicalSpecies[i].parentID == i) {
				AliveMolecules.push_back(fChemicalSpecies[i].trackID);
			}
		}
	}
	return AliveMolecules;
}


std::vector<TsIRTConfiguration::TsMolecule> TsHybridIRT::GetSurvivingMoleculesWithMolID(G4int molID) {
	std::vector<TsIRTConfiguration::TsMolecule> AliveMolecules;
	for (auto& IndexAndMol: fChemicalSpecies) {
		G4int i = IndexAndMol.first;
		if (!fChemicalSpecies[i].reacted) {
			if (fChemicalSpecies[i].id == molID) {
				AliveMolecules.push_back(fChemicalSpecies[i]);
			}
		}
	}
	return AliveMolecules;
}


G4bool TsHybridIRT::MoleculeExists(G4int Index) {
	std::unordered_map<G4int,TsIRTConfiguration::TsMolecule>::iterator IT = fChemicalSpecies.find(Index);
	if (IT == fChemicalSpecies.end()) {
		return false;
	}

	return !fChemicalSpecies[Index].reacted;
}


G4int TsHybridIRT::CountSurvivingMolecules() {
	G4int Result = 0;
	for (auto& IndexAndMolecule:fChemicalSpecies) {
		if (IndexAndMolecule.second.reacted){
			Result++;
		}
	}
	return Result;
}


std::map<G4int,G4int>  TsHybridIRT::GetEscapeYields() {
	std::map<G4int,G4int> Result;
	for (auto& MoleculeTimeAndCount:fGValues) {
		G4String Name = MoleculeTimeAndCount.first;
		if (!(MoleculeTimeAndCount.second).empty()){
			G4double Time = std::prev(MoleculeTimeAndCount.second.end())->first;
			for (auto& IndexAndName:fMolecules) {
				if (IndexAndName.second == Name){
					Result[IndexAndName.first] = MoleculeTimeAndCount.second[Time];
				}
			}
		}
	}
	return Result;
}


std::map<G4int,G4int> TsHybridIRT::GetDeltaYields() {
	std::map<G4int, G4int> Result;
	for (auto& IndexTimeAndDelta:fDeltaGValues) {
		G4int Index = IndexTimeAndDelta.first;
		std::map<G4double,G4int> TimeAndDelta = IndexTimeAndDelta.second;
		if (!TimeAndDelta.empty()) {
			G4double Time = std::prev(TimeAndDelta.end())->first;
			Result[Index] = TimeAndDelta[Time];
		}
	}
	return Result;
}


void TsHybridIRT::UpdateGillespieDirect(G4double time) {
	// Update Propensities
	fHomogeneousReactIndex     = -1;
	fHomogeneousTimeStep       = -1;
	fTotalPropensityAtThisTime = 0;

	G4double M = 1.0e3*mole/(m*m*m);
	//G4double fCubicVolume = 27E-15;

	std::map<G4int, TsIRTConfiguration::TsMolecularReaction> reactions = fReactionConf->GetReactions();
	std::map<G4int,G4String> MolNames = fReactionConf->GetMoleculeNames();
	G4int tBin = fUtils->FindBin(time, fStepTimes);

	// Update the concentration Container
	for (auto& IndexAndName:MolNames) {
		G4int Index = IndexAndName.first;
		std::vector<G4int> NewMoleculesOfThisKind = GetAllSpeciesOfAKindAtATime(time,Index);
		
		for (size_t i = 0; i < NewMoleculesOfThisKind.size() ; i++) {
			fConcentrationsAtThisTime[Index].push_back(NewMoleculesOfThisKind[i]);
		}
	}

	// Re-calculate Propensities
	fPropensityAtThisTime.clear();
	for (auto& IndexAndReaction:reactions) {
		G4int Index = IndexAndReaction.first;
		TsIRTConfiguration::TsMolecularReaction reaction = IndexAndReaction.second;
		G4int type    = reaction.reactionType;
		G4int reactA  = reaction.reactorA;
		G4int reactB  = reaction.reactorB;
		G4double kobs = reaction.kobs * M / (6.022e23 * fCubicVolume);
		G4double concentration = reaction.concentration * 6.022e23 * fCubicVolume / M;

		G4int TotalReactA = fConcentrationsAtThisTime[reactA].size() + fHomogeneousConcentrations[reactA];
		G4int TotalReactB = fConcentrationsAtThisTime[reactB].size() + fHomogeneousConcentrations[reactB];

		if (type != 6) {
			if (reactA != reactB) {
				fPropensityAtThisTime[Index] = kobs * TotalReactA * TotalReactB;
			}
			else {
				fPropensityAtThisTime[Index] = kobs * TotalReactA * (TotalReactB-1);
			}
		}
		else {
			if (concentration > 0) {
				fPropensityAtThisTime[Index] = kobs * TotalReactA * concentration;
			}
			else {
				fPropensityAtThisTime[Index] = kobs * TotalReactA * (6.022e23 * fCubicVolume);
			}
		}

		if (fPropensityAtThisTime[Index] >= 0)
			fTotalPropensityAtThisTime  += fPropensityAtThisTime[Index];
	}

	if (fTotalPropensityAtThisTime <= 0) {
		return;
	}

	// Sample Reaction and Time
	G4double rand = G4UniformRand();
	G4double PartialProp = 0.0;

	for(auto& IndexAndReactions:reactions) {
		G4int Index = IndexAndReactions.first;
		if (fPropensityAtThisTime[Index] >= 0)
			PartialProp += fPropensityAtThisTime[Index];

		if (PartialProp > rand * fTotalPropensityAtThisTime) {
			fHomogeneousReactIndex = Index;
			fHomogeneousTimeStep   = (1.0/fTotalPropensityAtThisTime) * std::log(1.0 / G4UniformRand());
			return;
		}
	}
}


G4bool TsHybridIRT::DoGillespieDirect(G4double currentTime) {
	TsIRTConfiguration::TsMolecularReaction reaction = fReactionConf->GetReaction(fHomogeneousReactIndex);
	G4int ReactA  = reaction.reactorA;
	G4int ReactB  = reaction.reactorB;
	G4int Type    = reaction.reactionType;
	G4bool Random = !reaction.positionSensitive;
	std::vector<G4int> Products = reaction.products;
	G4int tBin = fUtils->FindBin(currentTime, fStepTimes);

	// Find the reactors
	G4int IndexA = -1;
	G4int IndexB = -1;

	// Random: Faster but less detailed, perfect for concentration aimed studies
	// Distance based: Slower but offers more detail, perfect for position sensitive studies like DNA damage
	if (Random) {
		IndexA = ChooseReactantRandomly(fConcentrationsAtThisTime[ReactA],fHomogeneousConcentrations[ReactA],currentTime);
		if (Type != 6) {
			IndexB = ChooseReactantRandomly(fConcentrationsAtThisTime[ReactB],fHomogeneousConcentrations[ReactB],currentTime);
		}
	}
	else {
		std::pair<G4int, G4int> Reactants = ChooseReactantsBasedOnDistance(fConcentrationsAtThisTime[ReactA],fConcentrationsAtThisTime[ReactB],fHomogeneousConcentrations[ReactA],fHomogeneousConcentrations[ReactB],Type,currentTime);
		IndexA = Reactants.first;
		IndexB = Reactants.second;
	}

	// Check if the reaction is Possible
	if ((IndexA == -1 || IndexB == -1) && Type != 6) {return false;}
	else if (IndexA == -1 && Type == 6) {return false;}

	// Decrease number of molecules from containers
	if (IndexA == -2) { fHomogeneousConcentrations[ReactA]--;}
	if (IndexB == -2) { fHomogeneousConcentrations[ReactB]--;}
	if ( 0 <= tBin ) {
		for ( size_t ti = tBin; ti < fStepTimes.size(); ti++ ) {
			fTheGvalue[ReactA][ti]--;
			if (Type != 6) {
				fTheGvalue[ReactB][ti]--;
			}

			if (fReportDelta) {
				fDeltaGValues[fHomogeneousReactIndex][fStepTimes[ti]]++;
			}
		}
	}

	G4bool CreateMolecule = false;
	G4int Scenario = 0;
	std::vector<G4ThreeVector> Positions;
	G4int dnaVolumeID = -1;
	G4int dnaBaseID   = -1;
	G4int dnaStrandID = -1;
	if ((IndexA >= 0 && IndexB >= 0) && !Random) {
		CreateMolecule = true;
		Positions = fReactionConf->GetPositionOfProducts(fChemicalSpecies[IndexA],fChemicalSpecies[IndexB], fHomogeneousReactIndex);
		Scenario = 1;

		if (fChemicalSpecies[IndexA].volumeID >= 0) {
			dnaVolumeID = fChemicalSpecies[IndexA].volumeID;
			dnaBaseID   = fChemicalSpecies[IndexA].baseID;
			dnaStrandID = fChemicalSpecies[IndexA].strandID;
		}

		else if (fChemicalSpecies[IndexB].volumeID >= 0) {
			dnaVolumeID = fChemicalSpecies[IndexB].volumeID;
			dnaBaseID   = fChemicalSpecies[IndexB].baseID;
			dnaStrandID = fChemicalSpecies[IndexB].strandID;
		}
	}
	else if (((IndexA >= 0 && IndexB == -1) || (IndexA >= 0 && IndexB == -2)) && !Random) {
		CreateMolecule = true;
		Positions = fReactionConf->GetBackgroundPositionOfProducts(fChemicalSpecies[IndexA],fHomogeneousReactIndex);
		Scenario = 2;
		if (fChemicalSpecies[IndexA].volumeID >= 0) {
			dnaVolumeID = fChemicalSpecies[IndexA].volumeID;
			dnaBaseID   = fChemicalSpecies[IndexA].baseID;
			dnaStrandID = fChemicalSpecies[IndexA].strandID;
		}		
	}
	else if ((IndexA == -2 && IndexB >= 0) && !Random) {
		CreateMolecule = true;
		Positions = fReactionConf->GetBackgroundPositionOfProducts(fChemicalSpecies[IndexB],fHomogeneousReactIndex);
		Scenario = 3;
		if (fChemicalSpecies[IndexB].volumeID >= 0) {
			dnaVolumeID = fChemicalSpecies[IndexA].volumeID;
			dnaBaseID   = fChemicalSpecies[IndexA].baseID;
			dnaStrandID = fChemicalSpecies[IndexA].strandID;
		}	
	}

	// Conduct the Reaction
	for ( size_t u = 0; u < Products.size(); u++ ) {
		if (fMolecules[Products[u]] == "") {continue;}

		// Increase number of molecules from containers
		if (CreateMolecule) {
			TsIRTConfiguration::TsMolecule aProd;
			aProd.id       = Products[u];
			aProd.position = Positions[u];
			aProd.time     = currentTime+fHomogeneousTimeStep;
			aProd.reacted  = false;
			aProd.trackID  = 0;
			aProd.isDNA    = false;
			aProd.isNew    = true;
			aProd.volumeID = dnaVolumeID;
			aProd.baseID   = dnaBaseID;
			aProd.strandID = dnaStrandID;
			aProd.spin     = -1;
			aProd.chemAlgo = 0; // Assigning 0 for instance of direct Gillespie
			G4int I = fUtils->FindBin(fNx, fXMin, fXMax, Positions[u].x());
			G4int J = fUtils->FindBin(fNy, fYMin, fYMax, Positions[u].y());
			G4int K = fUtils->FindBin(fNz, fZMin, fZMax, Positions[u].z());
			
			// Merge this specie with the track
			fUsed[fSpeciesIndex]                 = false;
			fSpaceBinned[I][J][K][fSpeciesIndex] = true;
			fChemicalSpecies[fSpeciesIndex]      = aProd;
			fSpeciesOfAKind[Products[u]][fSpeciesIndex] = true;
			fSpeciesIndex++;
		}
		else {
			fHomogeneousConcentrations[Products[u]]++;
		}
		if ( 0 <= tBin ) {
			for ( size_t ti = tBin; ti < fStepTimes.size(); ti++ )
				fTheGvalue[Products[u]][ti]++;
		}
	}

	if (IndexA >= 0) { RemoveMolecule(IndexA); }
	if (IndexB >= 0) { RemoveMolecule(IndexB); }

	return true;
}


void TsHybridIRT::UpdateGillespieTauLeaping(G4double time) {
	G4double eps = 0.03;

	// Retrieve Reaction Information
	G4double M = 1.0e3*mole/(m*m*m);
	//G4int tBin = fUtils->FindBin(time, fStepTimes);

	std::map<G4int, TsIRTConfiguration::TsMolecularReaction> reactions = fReactionConf->GetReactions();
	std::map<G4int,G4String> MolNames = fReactionConf->GetMoleculeNames();

	// Clean step containers
	fPropensityAtThisTime.clear();
	fMuForMolecule.clear();
	fSigmaForMolecule.clear();
	fHORForMolecule.clear();
	fGFactorForMolecule.clear();
	fLimitsForMolecule.clear();

	// Calculate Auxiliary Quantities
	for (auto& IndexAndReaction:reactions) {
		G4int Index   = IndexAndReaction.first;
		G4int type    = IndexAndReaction.second.reactionType;
		G4int reactA  = IndexAndReaction.second.reactorA;
		G4int reactB  = IndexAndReaction.second.reactorB;
		G4int concentrationA = fHomogeneousConcentrations[reactA]; //fTheGvalue[reactA][tBin];
		G4int concentrationB = fHomogeneousConcentrations[reactB]; //fTheGvalue[reactB][tBin];

		G4double OrderOfReaction = 0;
		std::vector<G4int> prods = IndexAndReaction.second.products;
		G4double kobs = IndexAndReaction.second.kobs * (M*s) / (6.022e23 * fCubicVolume);
		G4double concentration = IndexAndReaction.second.concentration * 6.022e23 * fCubicVolume / M;

		G4double partialProp = 0;
		G4int vij = 1;

		if (type != 6) {
			if (reactA != reactB) {
				partialProp = kobs * concentrationA * concentrationB;
				OrderOfReaction = 2;
			}
			else {
				partialProp = kobs * concentrationA * (concentrationB-1);
				OrderOfReaction = 2 + (1 / (concentrationA-1));
				vij = 2;
			}
		}
		else {
			if (concentration > 0) {
				partialProp = kobs * concentrationA * concentration;
				OrderOfReaction = 2;
			}
			else {
				partialProp = kobs * concentrationA * (6.022e23 * fCubicVolume);
				OrderOfReaction = 1;
			}
		}

		if (partialProp < 0) {
			continue;
		}

		fPropensityAtThisTime[Index] = partialProp;
		if (vij == 1) {
			fMuForMolecule[reactA]    += fPropensityAtThisTime[Index];
			fMuForMolecule[reactB]    += fPropensityAtThisTime[Index];
			fSigmaForMolecule[reactA] += fPropensityAtThisTime[Index];
			fSigmaForMolecule[reactB] += fPropensityAtThisTime[Index];
		}

		else if (vij == 2) {
			fMuForMolecule[reactA]    += 2 * fPropensityAtThisTime[Index];
			fSigmaForMolecule[reactA] += 4 * fPropensityAtThisTime[Index];	
		}

		for (size_t i = 0; i < prods.size(); i++){
			fMuForMolecule[prods[i]]    -= fPropensityAtThisTime[Index];
			fSigmaForMolecule[prods[i]] += fPropensityAtThisTime[Index];
		}

		if (OrderOfReaction >= fHORForMolecule[reactA]) {
			fHORForMolecule[reactA]     = OrderOfReaction;
			fGFactorForMolecule[reactA] = OrderOfReaction;
		}
		if (OrderOfReaction >= fHORForMolecule[reactB]) {
			fHORForMolecule[reactB]     = OrderOfReaction;
			fGFactorForMolecule[reactB] = OrderOfReaction;
		}
	}

	// Get The Tau Leap
	fTotalChanges = 0;
	fMaxForMolecule.clear();
	fMinForMolecule.clear();
	fChangesForMolecule.clear();
	G4double Tau = DBL_MAX;
	for (auto& IndexAndName:MolNames) {
		G4int Index = IndexAndName.first;
		G4double temp = 0;
		fMaxForMolecule[Index] = 1;
		fMinForMolecule[Index] = 1;
		
		if (fGFactorForMolecule[Index] > 0) {
			G4double temp = eps * fHomogeneousConcentrations[Index] / fGFactorForMolecule[Index];
			fMaxForMolecule[Index] = temp > 1 ? temp : 1;
		}

		if (fMuForMolecule[Index] == 0 || fSigmaForMolecule[Index] == 0) { continue; }

		G4double temp1 = fMaxForMolecule[Index] / fabs(fMuForMolecule[Index]);
		G4double temp2 = fMaxForMolecule[Index]*fMaxForMolecule[Index] / fSigmaForMolecule[Index];

		fMinForMolecule[Index] = temp1 < temp2 ? temp1 : temp2;

		if (Tau > fMinForMolecule[Index]) { Tau = fMinForMolecule[Index]; }
	}

	Tau *= 3.5;

	// Calculate Changes By Reactions
	for (auto& IndexAndReaction:reactions) {
		G4int Index  = IndexAndReaction.first;
		G4int ReactA = IndexAndReaction.second.reactorA;
		G4int ReactB = IndexAndReaction.second.reactorB;
		G4int Type   = IndexAndReaction.second.reactionType;
		G4int concentrationA = fHomogeneousConcentrations[ReactA];
		G4int concentrationB = fHomogeneousConcentrations[ReactB];

		fChangesForMolecule[Index] = CLHEP::RandPoisson::shoot(fPropensityAtThisTime[Index]*Tau);

		if (Type != 6) {
			fLimitsForMolecule[Index]  = concentrationA < concentrationB ? concentrationA : concentrationB;
		}
		else
			fLimitsForMolecule[Index]  = concentrationA;

		fChangesForMolecule[Index] = fChangesForMolecule[Index] < fLimitsForMolecule[Index] ? fChangesForMolecule[Index] : fLimitsForMolecule[Index];
		fTotalChanges += fChangesForMolecule[Index];
	}

	fHomogeneousTimeStep = Tau*s;
}


G4bool TsHybridIRT::DoGillespieTauLeaping(G4double time) {
	// Retrieve Reaction Information
	G4int tBin = fUtils->FindBin(time, fStepTimes);
	std::map<G4int, TsIRTConfiguration::TsMolecularReaction> reactions = fReactionConf->GetReactions();

	if (fTotalChanges <= 0) { return false; }

	for(auto& IndexAndReaction:reactions) {
		G4int Index   = IndexAndReaction.first;
		fHomogeneousReactIndex = Index;
		G4int ReactA  = IndexAndReaction.second.reactorA;
		G4int ReactB  = IndexAndReaction.second.reactorB;
		G4int Type    = IndexAndReaction.second.reactionType;
		std::vector<G4int> Prods = IndexAndReaction.second.products;
		G4int Changes = fChangesForMolecule[Index];
		G4bool PositionSensitive = IndexAndReaction.second.positionSensitive;
		std::map<G4int, G4String> moleculesNames = fReactionConf->GetMoleculeNames();

		if (Changes <= 0) { continue; }
		else if (fHomogeneousConcentrations[ReactA] <= 0) { continue; }
		else if (Type != 6 && fHomogeneousConcentrations[ReactB] <= 0) { continue;}

		if (ReactA == ReactB && 2*Changes > fHomogeneousConcentrations[ReactA]) {Changes = fHomogeneousConcentrations[ReactA]/2;}

		if (fHomogeneousConcentrations[ReactA] < Changes) {Changes = fHomogeneousConcentrations[ReactA];}
		else if (Type != 6 && fHomogeneousConcentrations[ReactB] < Changes) {Changes = fHomogeneousConcentrations[ReactB];}

		if (fHomogeneousConcentrations[ReactA] - Changes >= 0) {fHomogeneousConcentrations[ReactA] -= Changes;}
		else {fHomogeneousConcentrations[ReactA] = 0;}

		if (Type != 6) {
			if (fHomogeneousConcentrations[ReactB] - Changes >= 0) {fHomogeneousConcentrations[ReactB] -= Changes;}
			else {fHomogeneousConcentrations[ReactB] = 0;}
		}

		for (size_t i = 0; i < Prods.size(); i++) {
			fHomogeneousConcentrations[Prods[i]] += Changes;		
		}

		for ( size_t ti = tBin; ti < fStepTimes.size(); ti++ ) {
			if (Type != 6) {
				fTheGvalue[ReactA][ti] = fHomogeneousConcentrations[ReactA];
				fTheGvalue[ReactB][ti] = fHomogeneousConcentrations[ReactB];
			}
			else {
				fTheGvalue[ReactA][ti] = fHomogeneousConcentrations[ReactA];
			}

			for (size_t i = 0; i < Prods.size(); i++) {
				fTheGvalue[Prods[i]][ti] = fHomogeneousConcentrations[Prods[i]];
			}
			if ( fReportDelta ) {
				fDeltaGValues[Index][fStepTimes[ti]] += Changes;
			}
		}
	}

	return true;
}


G4bool TsHybridIRT::IsNextIRTStillPossible() {
	G4int idx = fIRTValues[0].second;
	G4int iM = fIRTMolA[idx];
	G4int jM = fIRTMolB[idx];

	if (!MoleculeExists(iM) || !MoleculeExists(jM)) {
		return false;
	}
	return true;
}

std::vector<G4int> TsHybridIRT::GetAllSpeciesOfAKindAtATime(G4double time, G4int react) {
	std::vector<G4int> count;
	for(auto& IndexAandTrue:fSpeciesOfAKind[react]) {
		G4int Index = IndexAandTrue.first;
		if (MoleculeExists(Index)) {
			if (time >= fChemicalSpecies[Index].time) { count.push_back(Index); }
		}
	}

	for (size_t i = 0; i < count.size(); i++) {
		fSpeciesOfAKind[react].erase(count[i]);
	}

	return count;
}

std::pair<G4int,G4int> TsHybridIRT::ChooseReactantsBasedOnDistance(std::vector<G4int> ReactA, std::vector<G4int> ReactB, G4int HomogeneousA, G4int HomogeneousB,G4int Type, G4double currentTime) {
	G4int IndexA = -1;
	G4int IndexB = -1;

	// Check if the reaction happens with the homogeneous background
	if (HomogeneousA > 0) {
		G4bool isHomogeneousA = G4UniformRand() < (HomogeneousA/(ReactA.size() + HomogeneousA));
		if (isHomogeneousA) IndexA = -2;
	}

	if (HomogeneousB > 0 || Type != 6) {
		G4bool isHomogeneousB = G4UniformRand() < (HomogeneousB/(ReactB.size() + HomogeneousB));
		if (isHomogeneousB) IndexB = -2;
	}

	if (IndexA == -2 || IndexB == -2) {
		if (IndexA != -2) {
			while (true) {
				IndexA = ReactA[ReactA.size()*G4UniformRand()];
				if (MoleculeExists(IndexA)) break;
			}
		}

		if (IndexB != -2 && Type != 6) {
			while (true) {
				IndexB = ReactB[ReactB.size()*G4UniformRand()];
				if (MoleculeExists(IndexB)) break;
			}
		}

		return std::make_pair(IndexA, IndexB);
	}

	// If both reactants are heterogeneous
	G4double DistanceOrTime = DBL_MAX;
	std::pair<G4int,G4int> Reactants;

	if (Type != 6) {
		G4double GoodTimeA = 0;
		G4double GoodTimeB = 0;
		G4ThreeVector GoodPosA = G4ThreeVector();
		G4ThreeVector GoodPosB = G4ThreeVector();
		for(size_t i = 0; i < ReactA.size(); i++) {
			G4int  indexA      = ReactA[i];
			if (!MoleculeExists(indexA)) {continue;}
			for (size_t j = 0; j < ReactB.size(); j++) {
				G4int  indexB      = ReactB[j];
				if (!MoleculeExists(indexB)) {continue;}
				else if (indexA == indexB) {continue;}
				G4double timeA     = fChemicalSpecies[indexA].time;
				G4double timeB     = fChemicalSpecies[indexB].time;
				G4ThreeVector posA = fChemicalSpecies[indexA].position;
				G4ThreeVector posB = fChemicalSpecies[indexB].position;

				G4double dtA = (currentTime+fHomogeneousTimeStep)-timeA;
				G4double dtB = (currentTime+fHomogeneousTimeStep)-timeB;

				if (dtA > 0)
					fReactionConf->Diffuse(fChemicalSpecies[IndexA],dtA);
				if (dtB > 0)
					fReactionConf->Diffuse(fChemicalSpecies[IndexB],dtB);

				G4double distance = (fChemicalSpecies[indexA].position-fChemicalSpecies[indexB].position).mag();
				if (distance < DistanceOrTime) {
					IndexA = indexA;
					IndexB = indexB;
					DistanceOrTime = distance;
					GoodTimeA = currentTime+fHomogeneousTimeStep;
					GoodTimeB = currentTime+fHomogeneousTimeStep;
					GoodPosA  = fChemicalSpecies[indexA].position;
					GoodPosB  = fChemicalSpecies[indexB].position;
				}
				fChemicalSpecies[indexA].time = timeA;
				fChemicalSpecies[indexB].time = timeB;
				fChemicalSpecies[indexA].position = posA;
				fChemicalSpecies[indexB].position = posB;
			}
		}
		fChemicalSpecies[IndexA].time = GoodTimeA;
		fChemicalSpecies[IndexB].time = GoodTimeB;
		fChemicalSpecies[IndexA].position = GoodPosA;
		fChemicalSpecies[IndexB].position = GoodPosB;
	}
	else {
		for(size_t i = 0; i < ReactA.size(); i++) {
			G4int  indexA = ReactA[i];
			if (!MoleculeExists(indexA)) {continue;}
			G4double timeA = fChemicalSpecies[indexA].time;
			if (timeA < DistanceOrTime) {
				IndexA = indexA;
				DistanceOrTime = timeA;

			}
		}
		G4double dt = (currentTime+fHomogeneousTimeStep)-DistanceOrTime;
		if (dt > 0)
			fReactionConf->Diffuse(fChemicalSpecies[IndexA],(currentTime+fHomogeneousTimeStep)-DistanceOrTime);
	}

	return std::make_pair(IndexA,IndexB);
}


G4int TsHybridIRT::ChooseReactantRandomly(std::vector<G4int> React, G4int Homogeneous,G4double Time) {
	G4int Heterogeneous = React.size();

	if ((Heterogeneous + Homogeneous) <= 0) {
		return -1;
	}

	G4double HomogeneousContribution = Homogeneous/(Heterogeneous + Homogeneous);
	if (G4UniformRand() < HomogeneousContribution) {
		return -2;
	}

	while (true) {
		G4int Index = React[React.size()*G4UniformRand()];
		if (MoleculeExists(Index)) {
			fReactionConf->Diffuse(fChemicalSpecies[Index],(Time+fHomogeneousTimeStep)-fChemicalSpecies[Index].time);
			return Index;
		}
	}
}


void TsHybridIRT::SetContainersForNextPulse() {
	fSpeciesIndex = 0;
	fChemicalSpecies.clear();
	fSpaceBinned.clear();
	fSpeciesOfAKind.clear();
	fUsed.clear();
	fSampleIRTatStart = true;

	fIRTIndex.clear();
	fIRTValues.clear();
	fIRTMolA.clear();
	fIRTMolB.clear();
	fIRTIsBackground.clear();
	fIRTOrigTime.clear();
	fIRTPositions.clear();
	fIndex = 0;
}

