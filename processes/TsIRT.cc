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

void TsIRT::runIRT(G4double, G4double, G4double, G4bool) {	
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

						G4double r = binReaction.reactionRadius;
						
						G4double p = binReaction.probabilityOfReaction;
						G4double timeJ = fChemicalSpecies[j].time;
						
						if ( fChemicalSpecies[i].time != timeJ ) continue; /// Only zero timers<---Best match Clifford et al.
						
						if ( (fChemicalSpecies[i].position - fChemicalSpecies[j].position).mag2() < r*r ) {
							fReactedByContact = fReactionConf->MakeReaction(fChemicalSpecies, fSpeciesIndex, fSpaceBinned,
																			fNx, fNy, fNz, fXMin, fXMax,
																			fYMin, fYMax, fZMin, fZMax,
																			fTheGvalue, fStepTimes, i,
																			j, indexOfReaction, fChemicalSpecies[i].time, p, fUsed, fContactProducts);
							
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

			G4int dnaVolumeID = -1;
			G4int dnaBaseID   = -1;
			G4int dnaStrandID = -1;
			
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
				products = fReactionConf->GetReactionProducts(indexOfReaction);//binReaction.products;
				
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
				products = fReactionConf->GetReactionProducts(indexOfReaction);//binReaction.products;

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

				if (fChemicalSpecies[iM].volumeID >= 0) {
					dnaVolumeID = fChemicalSpecies[iM].volumeID;
					dnaBaseID   = fChemicalSpecies[iM].baseID;
					dnaStrandID = fChemicalSpecies[iM].strandID;
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

void TsIRT::Clean() {
	fChemicalSpecies.clear();
	fConcentrations.clear();
	fTheGvalue.clear();
	fTheGvalueInVolume.clear();
	fSpaceBinned.clear();
	fUsed.clear();
	fSpeciesOfAKind.clear();
	
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

void TsIRT::RemoveMolecule(G4int Index) {
	G4int I = fUtils->FindBin(fNx, fXMin, fXMax, fChemicalSpecies[Index].position.x());
	G4int J = fUtils->FindBin(fNy, fYMin, fYMax, fChemicalSpecies[Index].position.y());
	G4int K = fUtils->FindBin(fNz, fZMin, fZMax, fChemicalSpecies[Index].position.z());
	fSpaceBinned[I][J][K].erase(Index);
	fChemicalSpecies.erase(Index);
	fUsed.erase(Index);
}

void TsIRT::SetContainersForNextPulse() {;}
