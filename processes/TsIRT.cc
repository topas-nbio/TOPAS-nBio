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

TsIRT::TsIRT(TsParameterManager* pM, G4String parmName): TsVIRTProcedure(pM,parmName){
	fReactionConf = new TsIRTConfiguration(parmName, pM);

	fSpeciesIndex = 0;
	fMoleculesName = fReactionConf->GetMoleculeNames();

	fChemVerbosity = 0;
	if ( pM->ParameterExists("Ts/ChemistryVerbosity"))
		fChemVerbosity =pM->GetIntegerParameter("Ts/ChemistryVerbosity");

	fHighTimeScavenger = true;
}


TsIRT::~TsIRT() {;}


void TsIRT::runIRT(G4double, G4double, G4double, G4bool) {
	Sampling();
	ConductReactions();
	UpdateGValues();
}


void TsIRT::Sampling(){
	if ( !fScorersInitialized )
		initializeScorers();

	fNx = G4int((fXMax-fXMin)/fBinWidth) == 0 ? 1 : G4int((fXMax-fXMin)/fBinWidth);
	fNy = G4int((fYMax-fYMin)/fBinWidth) == 0 ? 1 : G4int((fYMax-fYMin)/fBinWidth);
	fNz = G4int((fZMax-fZMin)/fBinWidth) == 0 ? 1 : G4int((fZMax-fZMin)/fBinWidth);

	for(auto it = fChemicalSpecies.begin(); it != fChemicalSpecies.end(); ++it){
		TsIRTConfiguration::TsMolecule aMol = (*it).second;
		G4ThreeVector position = aMol.position;

		G4int t = (*it).first; 

		G4int tBin = fUtils->FindBin(aMol.time, fStepTimes);
		if ( -1 < tBin ) {
			for ( int tbin = tBin; tbin < (int)fStepTimes.size(); tbin++ ) {
				if ( fTheGvalue.find(aMol.id) == fTheGvalue.end() ) {
					fTheGvalue[aMol.id][tbin] = 1;
				} else {
					fTheGvalue[aMol.id][tbin]++; 
				}
			}
		}
		G4int I = fUtils->FindBin(fNx, fXMin, fXMax, aMol.position.x());
		G4int J = fUtils->FindBin(fNy, fYMin, fYMax, aMol.position.y());
		G4int K = fUtils->FindBin(fNz, fZMin, fZMax, aMol.position.z());

		fSpaceBinned[I][J][K][t] = true;

		sampleReactions(t);
	}
}


void TsIRT::AddMolecule(TsIRTConfiguration::TsMolecule aMol) {
    G4ThreeVector position = aMol.position;
    
    if ( fXMin > position.x() ) fXMin = position.x();
    if ( fYMin > position.y() ) fYMin = position.y();
    if ( fZMin > position.z() ) fZMin = position.z();
    
    if ( fXMax < position.x() ) fXMax = position.x();
    if ( fYMax < position.y() ) fYMax = position.y();
    if ( fZMax < position.z() ) fZMax = position.z();

    fChemicalSpecies[fSpeciesIndex] = aMol;
    fSpeciesIndex++;
}


void TsIRT::AddMolecule(G4Track* aTrack, G4double time, G4int moleculeID, G4ThreeVector offset){
	TsIRTConfiguration::TsMolecule aMol = ConstructMolecule(aTrack, time, 0, G4ThreeVector());

	G4ThreeVector position = aMol.position;

	if ( fXMin > position.x() ) fXMin = position.x();
	if ( fYMin > position.y() ) fYMin = position.y();
	if ( fZMin > position.z() ) fZMin = position.z();

	if ( fXMax < position.x() ) fXMax = position.x();
	if ( fYMax < position.y() ) fYMax = position.y();
	if ( fZMax < position.z() ) fZMax = position.z();

	fChemicalSpecies[fSpeciesIndex] = aMol;
	fSpeciesIndex++;
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
		// trackID != 0 reserved for DNA, < 0 = left strand, > 0 = right strand
		aMol.trackID = moleculeID;
		aMol.isDNA = false;

		if ( pdg == 1 | pdg == 5 ) aMol.spin = G4UniformRand() > 0.5 ? 1 : 0;
		else aMol.spin = -1;
	}

	return aMol;
}


void TsIRT::sampleReactions(G4int i) {

    if ( fChemicalSpecies[i].isDNA )
        return;

	if ( !MoleculeExists(i) )
		return;

	// background contact reaction
	G4int indexOf1stOrder = fReactionConf->ContactFirstOrderAndBackgroundReactions(fChemicalSpecies[i]);
	if ( -1 < indexOf1stOrder ) {
		AddToIRT(fChemicalSpecies[i].time,indexOf1stOrder,i,i,fChemicalSpecies[i].time,fChemicalSpecies[i].position,fChemicalSpecies[i].position,true);
		return;
	}

	G4ThreeVector thisPos = fChemicalSpecies[i].position;

	FindBinIndexes(thisPos, fRCutOff);

	for ( int ii = fxiniIndex; ii <= fxendIndex; ii++ ) {
		for ( int jj = fyiniIndex; jj <= fyendIndex; jj++ ) {
			for ( int kk = fziniIndex; kk <= fzendIndex; kk++ ) {
				for (auto& IndexAndAlive:fSpaceBinned[ii][jj][kk]) {
					G4int j = IndexAndAlive.first;

					if (j == i) continue;

					if ( !MoleculeExists(j) )
						continue;

					G4int indexOfReaction = fReactionConf->GetReactionIndex(fChemicalSpecies[i].id,
																			fChemicalSpecies[j].id);
					if( -1 < indexOfReaction ) {
						G4double timeI = fChemicalSpecies[i].time;
						G4double timeJ = fChemicalSpecies[j].time;
						G4double dt = fChemicalSpecies[i].time - timeJ;
						// if ( dt < 0 ) continue; // No efect if all initial times are the same, but affects for flash

						G4ThreeVector origPositionI = fChemicalSpecies[i].position;
						G4ThreeVector origPositionJ = fChemicalSpecies[j].position;

						// sampling contact reaction
						G4double r = fReactionConf->GetReaction(indexOfReaction).reactionRadius;
						if ( (timeI == timeJ) && (origPositionI - origPositionJ).mag2() < r*r ) {
							G4double p = fReactionConf->GetReaction(indexOfReaction).probabilityOfReaction;
							if(G4UniformRand() < p)	AddToIRT(timeI,indexOfReaction,i,j,fChemicalSpecies[i].time,fChemicalSpecies[i].position,fChemicalSpecies[j].position,false);
						}


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
					
					positions = fReactionConf->ResampleReactantsPosition(fChemicalSpecies[iM], fChemicalSpecies[jM],
															 indexOfReaction, irt);
					
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
				if (binReaction.index < 0) {continue;} // <-- why this would be < 0 if indexOfReactione exists?
				products = fReactionConf->GetReactionProducts(indexOfReaction);
				
				if(fChemVerbosity > 0){
				if ( fMoleculesName[fChemicalSpecies[iM].id]  == "OH^0" || fMoleculesName[fChemicalSpecies[jM].id]  == "OH^0" )
					G4cout<<"At time: "<<irt-fMinTime<<" ns"<<'\t'<<fMoleculesName[fChemicalSpecies[iM].id] << "(" << iM << ")" 
								<<" + "<<fMoleculesName[fChemicalSpecies[jM].id]<< "(" << jM << ")" << " => " << G4endl;
				}

				if ( 0 <= tBin ) {
					for ( int ti = tBin; ti < (int)fStepTimes.size(); ti++ ) {
						fTheGvalue[fChemicalSpecies[iM].id][ti]--;
						fTheGvalue[fChemicalSpecies[jM].id][ti]--;
						if ( fReportDelta )
							fDeltaGValues[indexOfReaction][fStepTimes[ti]]++; // <-- there should be a DeltaG for GvalueInVolume too
						
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
				
				fConcentrations[fChemicalSpecies[iM].id]--;
				fConcentrations[fChemicalSpecies[jM].id]--;
				RemoveMolecule(iM);
				RemoveMolecule(jM);
				
			} else {
				fReactionConf->Diffuse(fChemicalSpecies[iM],irt-fChemicalSpecies[iM].time);
				for ( int ip = 0; ip < 3; ip++ )
					positions.push_back(fChemicalSpecies[iM].position);
				
				binReaction = fReactionConf->GetReaction(indexOfReaction);
				if (binReaction.index < 0) {continue;} // why this would be < 0 if indexOfReaction exists?
				products = fReactionConf->GetReactionProducts(indexOfReaction);

				if ( 0 <= tBin ) {
					for ( int ti = tBin; ti < (int)fStepTimes.size(); ti++ ) {
						fTheGvalue[fChemicalSpecies[iM].id][ti]--;
						if ( fReportDelta )
							fDeltaGValues[indexOfReaction][fStepTimes[ti]]++; // there should be a DeltaG for GValueinVolume too
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

				if(fChemVerbosity > 0){
					G4int molB = binReaction.reactorB;
                                        if ( fMoleculesName[fChemicalSpecies[iM].id]  == "OH^0" )
					G4cout<<"  VI: At rest: "<<irt-fMinTime<<" ns"<<'\t'<<fMoleculesName[fChemicalSpecies[iM].id]<<" ("<< iM << ") (" << fChemicalSpecies[iM].id<<") "
								<<" + "<<fMoleculesName[molB]<<" ("<<molB<<") => " << G4endl;
				}

				if (fChemicalSpecies[iM].volumeID >= 0) {
					dnaVolumeID = fChemicalSpecies[iM].volumeID;
					dnaBaseID   = fChemicalSpecies[iM].baseID;
					dnaStrandID = fChemicalSpecies[iM].strandID;
				}
				
				fConcentrations[fChemicalSpecies[iM].id]--;
				RemoveMolecule(iM);
			}
			
			// 5.3 create new products at positions resampled in previous step
			for ( size_t u = 0; u < products.size(); u++ ) {
				TsIRTConfiguration::TsMolecule aProd;
				aProd.id = products[u];
				aProd.position = positions[u];
				aProd.time = irt;
				aProd.trackID = 0;
				aProd.isDNA = false;
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

                if (fChemVerbosity && u != 0)
                {
                    G4cout << " + ";
                }
                if (fChemVerbosity)
                {
                    G4cout << fMoleculesName[aProd.id];
                }

				sampleReactions(newID);

			}
            
            if (fChemVerbosity) {
                G4cout << G4endl;
            }

			if(fChemVerbosity >= 2){
				if (products.size() == 0)
				{
					G4cout << "No product"<<G4endl;
				}else{
	                G4cout << G4endl;
				}
			}

		} else {
			break;
		}
		
		if (fChemVerbosity > 2 && ( currentIRT >= fChunk && currentIRT % fChunk == 0 )) {
			G4cout << "     ---- current reactions performed " << currentIRT << " out " << fChunk*10 << G4endl;
		}
		RemoveFirstIRTElement();
	}
}


std::vector<TsIRTConfiguration::TsMolecule> TsIRT::GetSurvivingMoleculesWithMolID(G4int molID) {
    std::vector<TsIRTConfiguration::TsMolecule> aliveMolecules;
    for (auto& indexAndMol: fChemicalSpecies) {
        G4int i = indexAndMol.first;
        if (!fChemicalSpecies[i].reacted) {
            if (fChemicalSpecies[i].id == molID) {
                aliveMolecules.push_back(fChemicalSpecies[i]);
            }
        }
    }
    return aliveMolecules;
}


void TsIRT::Clean() {
	fChemicalSpecies.clear();
	fConcentrations.clear();
	fTheGvalue.clear();
	fTheGvalueInVolume.clear();
	fSpaceBinned.clear();
	
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

G4bool TsIRT::MoleculeExists(G4int Index) {
	std::unordered_map<G4int,TsIRTConfiguration::TsMolecule>::iterator IT = fChemicalSpecies.find(Index);
	if (IT == fChemicalSpecies.end()) {
		return false;
	}

	return true;
}


void TsIRT::RemoveMolecule(G4int Index) {
	G4int I = fUtils->FindBin(fNx, fXMin, fXMax, fChemicalSpecies[Index].position.x());
	G4int J = fUtils->FindBin(fNy, fYMin, fYMax, fChemicalSpecies[Index].position.y());
	G4int K = fUtils->FindBin(fNz, fZMin, fZMax, fChemicalSpecies[Index].position.z());
	fSpaceBinned[I][J][K].erase(Index);
	fChemicalSpecies.erase(Index);
}


void TsIRT::SetContainersForNextPulse() {;}
