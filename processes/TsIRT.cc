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

TsIRT::TsIRT(TsParameterManager* pM, G4String parmName)
: fPm(pM), fName(parmName), fHighTimeScavenger(false), fVerbosity(0), fScorersInitialized(false)
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
	
	G4cout << "\n -------------- OpenTOPAS-nBio IRT Warning -------------" << G4endl;
	
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
	
	fCurrentTimeScale = 1*us;
	fRCutOff = fReactionConf->GetRCutOff(fCurrentTimeScale);
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
	
	fTestIsInside = false;
	if (fPm->ParameterExists(GetFullParmName("OnlyIncludeChemicalSpeciesFromVirtualRegion/HLX"))) {
		fDx = fPm->GetDoubleParameter(GetFullParmName("OnlyIncludeChemicalSpeciesFromVirtualRegion/HLX"),"Length");
		fDy = fPm->GetDoubleParameter(GetFullParmName("OnlyIncludeChemicalSpeciesFromVirtualRegion/HLY"),"Length");
		fDz = fPm->GetDoubleParameter(GetFullParmName("OnlyIncludeChemicalSpeciesFromVirtualRegion/HLZ"),"Length");
		fTestIsInside = true;
	}

	
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
}


TsIRT::~TsIRT() {
	delete fReactionConf;
	delete fUtils;
}


G4bool TsIRT::Inside(G4ThreeVector p ) {
	G4double delta = 1 * nm;
	G4double dist = std::max(std::max(
									  std::abs(p.x())-fDx,
									  std::abs(p.y())-fDy),
							 std::abs(p.z())-fDz);
	if (dist > delta) return false;
	return !(dist > -delta);
}


void TsIRT::AddMolecule(TsIRTConfiguration::TsMolecule aMol) {
	if ( !fScorersInitialized )
		initializeScorers();

	G4ThreeVector position = aMol.position;
	
	fConcentrations[aMol.id]++;
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
		
		if ( fTestIsInside ) {
			if ( Inside(aMol.position)) {
				for ( int tbin = tBin; tbin < (int)fStepTimes.size(); tbin++ ) {
					if ( fTheGvalueInVolume.find(aMol.id) == fTheGvalueInVolume.end() )
						fTheGvalueInVolume[aMol.id][tbin] = 1;
					else
						fTheGvalueInVolume[aMol.id][tbin]++;
				}
			}
		}
	}
	
	if ( fXMin > position.x() ) fXMin = position.x();
	if ( fYMin > position.y() ) fYMin = position.y();
	if ( fZMin > position.z() ) fZMin = position.z();
	
	if ( fXMax < position.x() ) fXMax = position.x();
	if ( fYMax < position.y() ) fYMax = position.y();
	if ( fZMax < position.z() ) fZMax = position.z();
}


void TsIRT::AddMolecule(G4Step* aStep, G4double time, G4int moleculeID, G4ThreeVector offset) {
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


void TsIRT::AddMolecule(G4Track* aTrack, G4double time, G4int moleculeID, G4ThreeVector offset, G4bool isDNA) {
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


void TsIRT::AddMolecule(G4Track* aTrack, G4double time, G4int moleculeID, G4ThreeVector offset) {
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


void TsIRT::AddMolecule(const G4Track* aTrack, G4double time, G4int moleculeID, G4ThreeVector offset) {
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


TsIRTConfiguration::TsMolecule TsIRT::ConstructMolecule(G4Track* aTrack, G4double time, G4int moleculeID, G4ThreeVector offset) {
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


void TsIRT::Clean() {
	fChemicalSpecies.clear();
	fConcentrations.clear();
	fTheGvalue.clear();
	fTheGvalueInVolume.clear();
	fSpaceBinned.clear();
	fUsed.clear();
	
	fGValues.clear();
	fReactedDNA.clear();
	
	if (fReportDelta)
		fDeltaGValues.clear();
	
	if (fTestIsInside)
		fGValuesInVolume.clear();
	
	fScorersInitialized = false;
	fSampleIRTatStart = true;
	fCurrentTimeScale = 1*us;
	fMinTime = 0;
	fTimeCut = fStepTimes[fStepTimes.size()-1];
}


void TsIRT::FindBinIndexes(G4ThreeVector thisPos, G4double rcutOff) {
	fxiniIndex = fUtils->FindBin(fNx, fXMin, fXMax, thisPos.x()-rcutOff);
	fxendIndex = fUtils->FindBin(fNx, fXMin, fXMax, thisPos.x()+rcutOff);
	fyiniIndex = fUtils->FindBin(fNy, fYMin, fYMax, thisPos.y()-rcutOff);
	fyendIndex = fUtils->FindBin(fNy, fYMin, fYMax, thisPos.y()+rcutOff);
	fziniIndex = fUtils->FindBin(fNz, fZMin, fZMax, thisPos.z()-rcutOff);
	fzendIndex = fUtils->FindBin(fNz, fZMin, fZMax, thisPos.z()+rcutOff);
}


void TsIRT::initializeScorers() {
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
	
	if ( fTestIsInside ) {
		for ( auto& nameAndID : fMolecules ) {
			G4int id = nameAndID.first;
			for ( int i = 0; i < (int)fStepTimes.size(); i++ ) {
				fTheGvalueInVolume[id][i] = 0.0;
			}
		}
	}
	
	fScorersInitialized = true;
}


void TsIRT::contactReactions(G4int i,std::unordered_map<G4int, G4bool> used) {
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
														indexOf1stOrder, fChemicalSpecies[i].time, fUsed, fContactProducts);
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
			fConcentrations[fChemicalSpecies[i].id]--;
			RemoveMolecule(i);
			for (size_t prods = 0; prods < fContactProducts.size(); prods++) {
				if (fMolecules[fContactProducts[prods]] != "")
					fConcentrations[fContactProducts[prods]]++;
			}
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
						G4double r = binReaction.reactionRadius;
						
						G4double p = binReaction.probabilityOfReaction;
						G4double timeJ = fChemicalSpecies[j].time;
						
						if ( fChemicalSpecies[i].time != timeJ ) continue; /// Only zero timers<---Best match Clifford et al.
						
						if ( (fChemicalSpecies[i].position - fChemicalSpecies[j].position).mag2() < r*r ) {
							if ( !fTestIsInside )
								fReactedByContact = fReactionConf->MakeReaction(fChemicalSpecies, fSpeciesIndex, fSpaceBinned,
																				fNx, fNy, fNz, fXMin, fXMax,
																				fYMin, fYMax, fZMin, fZMax,
																				fTheGvalue, fStepTimes, i,
																				j, indexOfReaction, fChemicalSpecies[i].time, p, fUsed, fContactProducts);
							else
								fReactedByContact = fReactionConf->MakeReaction(fChemicalSpecies, fSpeciesIndex, fSpaceBinned,
																				fNx, fNy, fNz, fXMin, fXMax,
																				fYMin, fYMax, fZMin, fZMax,
																				fTheGvalue, fTheGvalueInVolume, fStepTimes, i,
																				j, indexOfReaction, fChemicalSpecies[i].time, p, fUsed, fContactProducts);
							
							if ( fReactedByContact ) {
								if ( fChemicalSpecies[j].isDNA )
									fReactedDNA[fChemicalSpecies[j].trackID] = std::make_pair(fChemicalSpecies[j].id,fChemicalSpecies[i].id);
								
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
								fConcentrations[fChemicalSpecies[i].id]--;
								fConcentrations[fChemicalSpecies[j].id]--;
								RemoveMolecule(i);
								RemoveMolecule(j);
								for (size_t prods = 0; prods < fContactProducts.size(); prods++) {
									if (fMolecules[fContactProducts[prods]] != "")
										fConcentrations[fContactProducts[prods]]++;
								}
								return;
							}
						}
					}
				}
			}
		}
	}
}


void TsIRT::sampleReactions(G4int i) {
	if ( !MoleculeExists(i) )
		return;
	
	if ( fChemicalSpecies[i].isDNA )
		return;
	
	if ( !fChemicalSpecies[i].isNew )
		return;

	G4ThreeVector thisPos = fChemicalSpecies[i].position;

	FindBinIndexes(thisPos, fRCutOff);

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
						// if ( dt < 0 ) continue; // No efect if all initial times are the same, but affects for flash
						
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
						
						G4double irt = fReactionConf->GetIndependentReactionTime(fChemicalSpecies[i],
																				 fChemicalSpecies[j],
																				 indexOfReaction);
						G4double gTime = fChemicalSpecies[i].time;

						if ( (0 < irt) && (irt + gTime <= fTimeCut) && (fMinTime <= irt + gTime)) {
							irt += gTime;
							AddToIRT(irt,indexOfReaction,i,j,fChemicalSpecies[i].time,fChemicalSpecies[i].position,fChemicalSpecies[j].position,false);
							
						}
						
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
	
	std::pair<G4int, G4double> tnIscav = fReactionConf->SampleIRTFirstOrderAndBackgroundReactions(fChemicalSpecies[i]);
	G4double tscav = tnIscav.second;
	G4double gTime = fChemicalSpecies[i].time;

	if (gTime + tscav < fMinTime) {
		gTime = fMinTime;
	}

	if ((fHighTimeScavenger && 0 < tscav) && (tscav+gTime <= fTimeCut)) {
		tscav += gTime;
		AddToIRT(tscav,tnIscav.first,i,i,fChemicalSpecies[i].time,fChemicalSpecies[i].position,fChemicalSpecies[i].position,true);
	}
}


void TsIRT::runIRT(G4double startTime, G4double finalTime) {
	fMinTime = 0;
	fTimeCut = fStepTimes[fStepTimes.size()-1];
	if (startTime != -1)
		fMinTime = startTime;
	if (finalTime != -1)
		fTimeCut = finalTime;
		
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


void TsIRT::VoxelizeAndSortSpace() {
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

void TsIRT::TestForContactReactions() {
	if (!fTestForContactReactions ) { return ; }

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

void TsIRT::SampleIndependantReactionTimes() {
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

void TsIRT::SortIndependantReactionTimes() {
	if (!fSortByTime)
		std::sort(fIRTValues.begin(), fIRTValues.end());
}


void TsIRT::ConductReactions() {
	fSampleIRTatStart = false;

	while (fIRTValues.size() > 0) {
		G4int currentIRT = 0;
		G4double irt = fIRTValues[currentIRT].first;

		if ( irt < fTimeCut ) {
			G4int idx = fIRTValues[currentIRT].second;
			G4int iM = fIRTMolA[idx];
			G4int jM = fIRTMolB[idx];

			if (!MoleculeExists(iM) || !MoleculeExists(jM)) {
				RemoveFirstIRTElement();
				continue;
			}
			
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
				if ( !fChemicalSpecies[iM].isDNA && !fChemicalSpecies[jM].isDNA ) { // None is DNA
					fChemicalSpecies[jM].position = fIRTPositions[idx].second;
					fChemicalSpecies[jM].time = fIRTOrigTime[idx];
					
					fReactionConf->ResampleReactantsPosition(fChemicalSpecies[iM], fChemicalSpecies[jM],
															 indexOfReaction, irt);
					
					positions = fReactionConf->GetPositionOfProducts(fChemicalSpecies[iM],
																	 fChemicalSpecies[jM], indexOfReaction);
				} else { // at least one is DNA. Set product positions as the DNA molecule position.
					// Score the base pair ID, and molecule that caused the damage.
					if ( fChemicalSpecies[iM].isDNA ) {
						for ( int ip = 0; ip < 3; ip++ )
							positions.push_back(fChemicalSpecies[iM].position);
						
						fReactedDNA[fChemicalSpecies[iM].trackID] = std::make_pair(fChemicalSpecies[iM].id,fChemicalSpecies[jM].id);
					} else {
						for ( int ip = 0; ip < 3; ip++ )
							positions.push_back(fChemicalSpecies[jM].position);
						
						fReactedDNA[fChemicalSpecies[jM].trackID] = std::make_pair(fChemicalSpecies[jM].id,fChemicalSpecies[iM].id);
					}
				}
				
				binReaction = fReactionConf->GetReaction(indexOfReaction);
				if (binReaction.index < 0) {continue;}
				products = binReaction.products;
				
				if ( 0 <= tBin ) {
					for ( int ti = tBin; ti < (int)fStepTimes.size(); ti++ ) {
						fTheGvalue[fChemicalSpecies[iM].id][ti]--;
						fTheGvalue[fChemicalSpecies[jM].id][ti]--;
						if ( fReportDelta )
							fDeltaGValues[indexOfReaction][fStepTimes[ti]]++;
						
					}
					
					if ( fTestIsInside ) {
						if ( Inside(fChemicalSpecies[iM].position)) {
							for ( int ti = tBin; ti < (int)fStepTimes.size(); ti++ ) {
								fTheGvalueInVolume[fChemicalSpecies[iM].id][ti]--;
								fTheGvalueInVolume[fChemicalSpecies[jM].id][ti]--;
							}
							reactionIsInside = true;
						}
					}
				}
				
				fChemicalSpecies[iM].reacted = true;
				fChemicalSpecies[jM].reacted = true;
				fConcentrations[fChemicalSpecies[iM].id]--;
				fConcentrations[fChemicalSpecies[jM].id]--;
				RemoveMolecule(iM);
				RemoveMolecule(jM);
				
			} else {
				fReactionConf->Diffuse(fChemicalSpecies[iM],irt-fChemicalSpecies[iM].time);
				for ( int ip = 0; ip < 3; ip++ )
					positions.push_back(fChemicalSpecies[iM].position);
				
				binReaction = fReactionConf->GetReaction(indexOfReaction);
				if (binReaction.index < 0) {continue;}
				products = binReaction.products;
				if ( 0 <= tBin ) {
					for ( int ti = tBin; ti < (int)fStepTimes.size(); ti++ ) {
						fTheGvalue[fChemicalSpecies[iM].id][ti]--;
						if ( fReportDelta )
							fDeltaGValues[indexOfReaction][fStepTimes[ti]]++;
					}
					if ( fTestIsInside ) {
						if ( Inside(fChemicalSpecies[iM].position)) {
							for ( int ti = tBin; ti < (int)fStepTimes.size(); ti++ ) {
								fTheGvalueInVolume[fChemicalSpecies[iM].id][ti]--;
							}
							reactionIsInside = true;
						}
					}
				}
				
				fChemicalSpecies[iM].reacted = true;
				fConcentrations[fChemicalSpecies[iM].id]--;
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
				fConcentrations[aProd.id]++;

				if (fMolecules[aProd.id] == "") {continue;}
				
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
				fSpeciesIndex++;
				
				if ( 0 <= tBin ) {
					for ( int ti = tBin; ti < (int)fStepTimes.size(); ti++ )
						fTheGvalue[aProd.id][ti]++;
					
					if ( fTestIsInside ) {
						if ( reactionIsInside) {
							for ( int ti = tBin; ti < (int)fStepTimes.size(); ti++ )
								fTheGvalueInVolume[aProd.id][ti]++;
						}
					}
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
}

void TsIRT::UpdateGValues() {
	for ( auto& nameAndID : fMolecules ) {
		G4String name = nameAndID.second;
		G4int id = nameAndID.first;
		for ( auto& timeAndGvalue : fTheGvalue[id] ) {
			G4int tBin = timeAndGvalue.first;
			G4double time = fStepTimes[tBin];
			fGValues[name][time] = timeAndGvalue.second;
		}
		if ( fTestIsInside ) {
			for ( auto& timeAndGvalue : fTheGvalueInVolume[id] ) {
				G4int tBin = timeAndGvalue.first;
				G4double time = fStepTimes[tBin];
				fGValuesInVolume[name][time] = timeAndGvalue.second;
			}
		}
	}
}

void TsIRT::CleanIRTVariables() {
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
	
	fCurrentTimeScale = 1*us;
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


void TsIRT::AddToIRT(G4double Time, G4int Reaction, G4int molA, G4int molB, G4double OrigTime, G4ThreeVector aPos, G4ThreeVector bPos, G4bool isBack) {
	AddIRTinAscendantOrder(Time);
	fIRTIndex[fIndex] = Reaction;
	fIRTMolA[fIndex]  = molA;
	fIRTMolB[fIndex]  = molB;
	fIRTOrigTime[fIndex]     = OrigTime;
	fIRTPositions[fIndex]    = std::make_pair(aPos, bPos);
	fIRTIsBackground[fIndex] = isBack;
	fIndex++;	
}


void TsIRT::AddIRTinAscendantOrder(G4double Time) {
	G4int N = fIRTValues.size();
	std::pair<G4double,G4int> newPair = std::make_pair(Time,fIndex);
	if (N == 0) {
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


void TsIRT::RemoveFirstIRTElement() {
	G4int index = fIRTValues[0].second;
	fIRTValues.erase(fIRTValues.begin());
	fIRTIndex.erase(index);
	fIRTMolA.erase(index);
	fIRTMolB.erase(index);
	fIRTOrigTime.erase(index);
	fIRTPositions.erase(index);
	fIRTIsBackground.erase(index);
}


void TsIRT::RemoveMolecule(G4int Index) {
	G4int I = fUtils->FindBin(fNx, fXMin, fXMax, fChemicalSpecies[Index].position.x());
	G4int J = fUtils->FindBin(fNy, fYMin, fYMax, fChemicalSpecies[Index].position.y());
	G4int K = fUtils->FindBin(fNz, fZMin, fZMax, fChemicalSpecies[Index].position.z());
	fSpaceBinned[I][J][K].erase(Index);
	fChemicalSpecies.erase(Index);
	fUsed.erase(Index);
}


G4String TsIRT::GetFullParmName(G4String suffix ) {
	return "Sc/" + fName + "/" + suffix;
}

void TsIRT::SetGValues(std::map<G4String, std::map<G4double, G4int>> gval) {
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


void TsIRT::SetDeltaGValues(std::map<G4int, std::map<G4double, G4int>> dgval) {
	fDeltaGValues = dgval;
}


std::map<G4String, std::map<G4double, G4int>> TsIRT::GetGValues() {
	if (!fTestIsInside)
		return fGValues;
	else 
		return fGValuesInVolume;
}


std::map<G4String, std::map<G4double, G4int>> TsIRT::GetGValuesInVolume() {
	return fGValuesInVolume;
}


std::map<G4int, std::map<G4double, G4int>> TsIRT::GetDeltaGValues() {
	return fDeltaGValues;
}


std::map<G4int, std::pair<G4int,G4int>> TsIRT::GetReactedDNA() {
	return fReactedDNA;
}


std::pair<G4String, G4String> TsIRT::GetReactants(G4int ReactionIndex) {
	return fReactionConf->GetReactants(ReactionIndex);
}


std::vector<G4String> TsIRT::GetProducts(G4int ReactionIndex) {
	return fReactionConf->GetProducts(ReactionIndex);
}

std::vector<G4int> TsIRT::GetSurvivingMolecules() {
	std::vector<G4int> AliveMolecules;
	for (auto& IndexAndMol: fChemicalSpecies) {
		G4int i = IndexAndMol.first;
		if (!fChemicalSpecies[i].reacted && fChemicalSpecies[i].trackID != 0) {
			AliveMolecules.push_back(fChemicalSpecies[i].trackID);
		}
	}
	return AliveMolecules;
}

G4bool TsIRT::MoleculeExists(G4int Index) {
	std::unordered_map<G4int,TsIRTConfiguration::TsMolecule>::iterator IT = fChemicalSpecies.find(Index);
	if (IT == fChemicalSpecies.end()) {
		return false;
	}

	return !fChemicalSpecies[Index].reacted;
}

G4int TsIRT::CountSurvivingMolecules() {
	G4int Result = 0;
	for (auto& IndexAndMolecule:fChemicalSpecies) {
		if (IndexAndMolecule.second.reacted){
			Result++;
		}
	}
	return Result;
}

void TsIRT::PrintConcentrations() {
	G4cout << "Concentrations: " << fConcentrations.size() << G4endl;
	for (auto& MoleculeTimeAndCount:fGValues) {
		G4String Name = MoleculeTimeAndCount.first;
		if (!(MoleculeTimeAndCount.second).empty()){
			G4double Time = std::prev(MoleculeTimeAndCount.second.end())->first;
			G4cout << Name << ":" << MoleculeTimeAndCount.second[Time] << G4endl;
		}
	}
}


std::map<G4int,G4int>  TsIRT::GetEscapeYields() {
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


std::map<G4int,G4int> TsIRT::GetDeltaYields() {
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
