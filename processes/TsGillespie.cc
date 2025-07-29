// Extra Class for TsEmDNAChemistry

#include "G4SystemOfUnits.hh"
#include "G4UIcommand.hh"

#include "TsGillespie.hh"
#include "TsParameterManager.hh"
#include "TsIRTUtils.hh"

TsGillespie::TsGillespie(TsParameterManager* pM, G4String parmName, TsVIRTProcedure* irt, G4double initTime, G4double endTime, G4int bins, G4double vol, G4double DR, G4double t)
: fPm(pM), fName(parmName)
{
	fIRT = irt;
	fMoleculeNames   = fIRT->GetIRTConfiguration()->GetMoleculeNames();

	fInitialTime = initTime;
	fFinalTime   = endTime;

	fUtils = fIRT->GetUtils();
	fTimeSteps = fIRT->GetStepTimes();

	fTotalNumberOfMolecules = 0;
	fpVolume = vol;

	fNoneIndex = -1000;
	for (auto& IndexAndName:fMoleculeNames) {
		G4int Index   = IndexAndName.first;
		G4String Name = IndexAndName.second;

		if (Name == "None") {
			fNoneIndex = Index;
			break;
		}
	}

	fNumberOfSavedStates  = 0;
	fPercentageDifference = 0.01;
	fNumberOfStatedToEquilibrium = 5000;
	fContainersInitialized = false;
	fTotalNumberOfMoleculesSum1 = 0;
	fTotalNumberOfMoleculesSum2 = 0;
	fPrintConcentrations = false;
	prop.clear();



	G4int Nreaction = fIRT->GetIRTConfiguration()->GetNumberOfReactions();
	fLastReactionIndex    = DBL_MAX;
	for(G4int i= 0; i<Nreaction;i++){
		prop.push_back(0);
	}

	G4int Nmolecules = fIRT->GetIRTConfiguration()->GetMoleculeNames().size();
	for(G4int i=0; i<Nmolecules; i++){
		fReverse.push_back(0);
	}

	fCurrentTime = 0;
	fPropensity = 0;
	fSpeciesIndex = 0;

	fCubicVolume = vol;
	fDoserate = DR;
	fIrradiationTime = t;
	fInsertPropensity = 0;
}

TsGillespie::~TsGillespie() {;}

void TsGillespie::Initialize() {

	fSpeciesIndex     = fIRT->GetSpeciesIndex();
	fMoleculesAtTime  = fIRT->GetGValues();
	fDeltaGPerReactionPerTime = fIRT->GetDeltaGValues();
//	fConcentrations   =  fIRT->GetEscapeYields();
/*
	for (auto& IndexAndConcentrations:fConcentrations) {
		//G4int Index = IndexAndConcentrations.first;
		G4int Conc  = IndexAndConcentrations.second;
		fTotalNumberOfMolecules += Conc;
	}*/

	// 3. Atomic Reaction Rate and Concentration Conversion
	RecoverReactionData();
}


void TsGillespie::RecoverReactionData() {
//	fCubicVolume = fpVolume;
//	fCubicVolume = (fCubicVolume / m3) * 1000;

	// unit remove
	fReactions     = fIRT->GetIRTConfiguration()->GetReactions();
	fMoleculeNames = fIRT->GetIRTConfiguration()->GetMoleculeNames();
	std::map<G4int, std::vector<std::pair<G4int,G4int>>> reactability = fIRT->GetIRTConfiguration()->GetReactability();

	for(auto it_reaction = fReactions.begin(); it_reaction != fReactions.end();){
		G4int Index  = it_reaction->first;
		if (fReactions[Index].reactionType != 6) {
			fReactions[Index].kobs /= fCubicVolume * CLHEP::Avogadro;
		}
		else {
			fReactions[Index].kobs = fReactions[Index].scavengingCapacity;
//			fReactions[Index].concentration = 0;

			/*
			if (fReactions[Index].concentration != 0 && fReactions[Index].kobs != 0
//					&& // first order reactions
//					fReactions[Index].reactorB != 6 && fReactions[Index].reactorB != 7 // acid based reactions
					) {

				fReactions[Index].concentration *= fCubicVolume;
				fReactions[Index].kobs /= fCubicVolume * CLHEP::Avogadro;

				if(fMoleculeNames[fReactions[Index].reactorB] != "None" && fMoleculeNames[fReactions[Index].reactorB] != "None2"){
					fConcentrations[fReactions[Index].reactorB] = fReactions[Index].concentration;
				}
				fReactions[Index].concentration = 0;
				fReactions[Index].kobs = 0;

				it_reaction = fReactions.erase(it_reaction);
				continue;

			}
//			else if(fReactions[Index].reactorB == 6 || fReactions[Index].reactorB == 7){
//				fReactions[Index].kobs *= fReactions[Index].concentration;
//			}
			else {
				fReactions[Index].kobs = fReactions[Index].scavengingCapacity;
				fReactions[Index].concentration = 0;

				if(fReactions[Index].kobs == 0) {
					it_reaction = fReactions.erase(it_reaction);
					continue;
				}
			}*/
		}

		fMoleculeCanReactWith[-fReactions[Index].reactorA].push_back(Index);
		if(fReactions[Index].reactorA != fReactions[Index].reactorB) fMoleculeCanReactWith[-fReactions[Index].reactorB].push_back(Index);

		std::vector<G4int> lastReact = fReactions[Index].products;
		lastReact.push_back(fReactions[Index].reactorA);

		if(fReactions[Index].reactionType != 6) lastReact.push_back(fReactions[Index].reactorB);

		for(auto it = lastReact.begin(); it != lastReact.end(); ++it){
			auto it2 = reactability[(*it)].begin();
			for(; it2 != reactability[(*it)].end(); ++it2){
				G4int id = it2->second;
				fMoleculeCanReactWith[Index].push_back(id);

			}
		}
		std::sort(fMoleculeCanReactWith[Index].begin(), fMoleculeCanReactWith[Index].end());
		fMoleculeCanReactWith[Index].erase(std::unique(fMoleculeCanReactWith[Index].begin(), fMoleculeCanReactWith[Index].end()), fMoleculeCanReactWith[Index].end());
		++it_reaction;

	}

/*
	for (auto& IndexAndReaction:fReactions) {
		G4int Index  = IndexAndReaction.first;
		G4cout <<Index<<" | Reaction |" << fMoleculeNames[fReactions[Index].reactorA] << " + " << fMoleculeNames[fReactions[Index].reactorB]<<" => ";
		std::vector<G4int> products = fReactions[Index].products;
		for(auto pro = products.begin(); pro != products.end(); ++pro){
			G4cout<<fMoleculeNames[(*pro)]<<'\t';
		}
		G4cout<<products.size()<<" | Type " << fReactions[Index].reactionType << " |" << G4endl;
		if (fReactions[Index].reactionType != 6) {
//			G4cout << "---- Molar Kobs=      " << fReactions[Index].kobs * (6.022e23 * fCubicVolume) * s << " /M/s " << G4endl;
			G4cout << "---- Molecular Kobs = " << fReactions[Index].kobs * s << " /s" << G4endl;
			G4cout << "---- Molar Kobs=      " << fReactions[Index].kobs * fCubicVolume * CLHEP::Avogadro * (1e-6 * mm3) * s << " /M/s " << G4endl;
		}
		else {
			G4cout << "---- | Dissociation Rate = " << fReactions[Index].kobs * s<< " /s" << G4endl;
		}


	}

	G4cout<<"Reactability: "<<G4endl;

	auto test = fMoleculeCanReactWith.begin();
	for(; test != fMoleculeCanReactWith.end(); ++test){
		G4int Index1 = test->first;

		if(Index1 < 0) G4cout<<"Last input: "<<Index1<<'\t'<<fMoleculeNames[-Index1]<<G4endl;
		else G4cout<<"Last reaction: "<<Index1<<'\t'<<fMoleculeNames[fReactions[Index1].reactorA]<<" + "<<fMoleculeNames[fReactions[Index1].reactorB]<<" => ";
		std::vector<G4int> products = fReactions[Index1].products;
		for(auto pro = products.begin(); pro != products.end(); ++pro){
			G4cout<<fMoleculeNames[(*pro)]<<'\t';
		}
		G4cout<<G4endl;

		auto test2 = test->second.begin();
		for(; test2 != test->second.end(); ++test2){
			G4int Index2 = (*test2);
			G4cout<<"Affected reactions: "<<Index2<<'\t'<<fMoleculeNames[fReactions[Index2].reactorA]<<'\t'<<fMoleculeNames[fReactions[Index2].reactorB]<<G4endl;
		}
	}*/
}


void TsGillespie::InitializeContainers() {
	if (fContainersInitialized) {return;}
	// 4.1 Number of molecules (GValues)
	for (auto& IndexAndName:fMoleculeNames) {
		G4String Name = IndexAndName.second;
		for (size_t i = 0; i < fTimeSteps.size(); i++) {
			fMoleculesAtTime[Name][fTimeSteps[i]] = 0;
		}
	}

	// 4.2. Delta G
	std::map<G4int, std::map<G4double, G4int>> DeltaG = fIRT->GetDeltaGValues();
	for (auto& IndexAndReaction:fReactions) {
		G4int Index = IndexAndReaction.first;
		for (size_t i = 0; i < fTimeSteps.size(); i++) {
			G4double Time = fTimeSteps[i];
			fDeltaGPerReactionPerTime[Index][Time] = 0;
		}
	}
	fContainersInitialized = true;
}


void TsGillespie::Clean() {
	fCurrentTime = fInitialTime;
	fTotalNumberOfMolecules = 0;
	fConcentrations.clear();
	fMoleculesAtTime.clear();
	fMoleculesFromIRT.clear();
	fDeltaGPerReactionPerTime.clear();
	fReactions.clear();
	fSpeciesIndex = 0;
	fPulseTimes.clear();
	fEscapeYieldsAtTime.clear();
	fDeltaYieldsAtTime.clear();
	fConcentrationsSum1.clear();
	fConcentrationsSum2.clear();
	fTotalNumberOfMoleculesSum1 = 0;
	fTotalNumberOfMoleculesSum2 = 0;
	fContainersInitialized = false;
}

std::pair<G4int, G4double> TsGillespie::Propensity() {
	if(fCurrentTime < fIrradiationTime)	fPropensity = fInsertPropensity;
	else fPropensity = 0;

	if(fLastReactionIndex > 10000){ // full sampling
		for (auto& IndexAndReaction:fReactions) {
			G4int i = IndexAndReaction.first;
			G4int Type    = fReactions[i].reactionType;
			G4int ReactA  = fReactions[i].reactorA;
			G4int ReactB  = fReactions[i].reactorB;
			G4double ConA = fConcentrations[ReactA];
			G4double ConB = fConcentrations[ReactB];
			G4double Kobs = fReactions[i].kobs;

			if (Type != 6) {
				if (ReactA == ReactB) {
					if(ConA < 2) prop[i] = 0;
					else prop[i] = Kobs * ConA * (ConA-1) * 0.5; // * (0.5);
				}
				else {
					if(ConA < 1 || ConB < 1) prop[i] = 0;
					else prop[i] =  Kobs * ConA * ConB;
				}
			}
			else {
//				prop[i] = 0;
				if(ConA < 1) prop[i] = 0;
				else prop[i] = Kobs * ConA; // Kobs * ConA;
			}

			fPropensity += prop[i];

//			G4cout<<"Full sampling: "<< i<<'\t'<<
//					fMoleculeNames[ReactA]<<" ("<<ConA<<") "<<
//					fMoleculeNames[ReactB]<<" ("<<ConB<<") "<<
//					Kobs<<'\t'<<
//					prop[i]<<'\t'<<fPropensity<<'\t'<<fInsertPropensity<<G4endl;
		}

	}else{ // partial sampling
		auto it = fMoleculeCanReactWith[fLastReactionIndex].begin();
		for(; it != fMoleculeCanReactWith[fLastReactionIndex].end(); ++it){
			G4int i = (*it);
			G4int Type    = fReactions[i].reactionType;
			G4int ReactA  = fReactions[i].reactorA;
			G4int ReactB  = fReactions[i].reactorB;
			G4double ConA = fConcentrations[ReactA];
			G4double ConB = fConcentrations[ReactB];
			G4double Kobs = fReactions[i].kobs;

			if (Type != 6) {
				if (ReactA == ReactB) {
					if(ConA < 2) prop[i] = 0;
					else prop[i] = Kobs * ConA * (ConA-1) * 0.5; // * (0.5);
				}
				else {
					if(ConA < 1 || ConB < 1) prop[i] = 0;
					else prop[i] =  Kobs * ConA * ConB;
				}
			}
			else {
//				prop[i] = 0;
				if(ConA < 1) prop[i] = 0;
				else prop[i] = Kobs * ConA; // Kobs * ConA;
			}

//			G4cout<<"Partial sampling: "<< i<<'\t'<<
//					fMoleculeNames[ReactA]<<" ("<<ConA<<") "<<
//					fMoleculeNames[ReactB]<<" ("<<ConB<<") "<<
//					Kobs<<'\t'<<
//					prop[i]<<'\t'<<fInsertPropensity<<G4endl;

//		    G4cout<<"Sampling: "<<i<<'\t'<<fReactions[i].reactionType<<'\t'<<fMoleculeNames[ReactA]<<" ("<<fConcentrations[ReactA]<<") + "<<fMoleculeNames[ReactB]<<" ("<<fConcentrations[ReactB]<<") "<<prop[i]<<'\t'<<fPropensity<<G4endl;
		}

		for(auto& IndexAndReaction: fReactions){
			G4int i = IndexAndReaction.first;
			fPropensity += prop[i];
		}
	}

	G4double rand = G4UniformRand() * fPropensity;
	G4double PartialProp = 0.0;
	G4int i= 0;

	for (auto& IndexAndReaction:fReactions) {
		i = IndexAndReaction.first;
		PartialProp += prop[i];
		if( rand < PartialProp) return std::make_pair(i, fPropensity);
	}

	for (auto& IndexAndInput:fMolPerTime){
		PartialProp += IndexAndInput.second;
		if( rand < PartialProp) return std::make_pair(-IndexAndInput.first, fPropensity);
	}

	return std::make_pair(-1, -1);
}

G4double TsGillespie::TimeIncrement(G4double Prop) {
	return (1.0/Prop) * std::log(1.0 / G4UniformRand());
}


G4double TsGillespie::GetNextEscapeYieldTime() {
	return fPulseTimes[0]/s;
}


G4int TsGillespie::GetMoleculeIndex(G4int Molecule) {
	for (auto& IndexAndSpecies:fMoleculesFromIRT) {
		G4int Index = IndexAndSpecies.first;
		if (fMoleculesFromIRT[Index].id == Molecule && fMoleculesFromIRT[Index].reacted == false) {
			return Index;
		}
	}
	return -1;
}


void TsGillespie::AddMolecule(G4int Molecule, G4int MolA, G4int MolB) {
	G4double Time = fMoleculesFromIRT[MolA].time;
	G4ThreeVector Position1 = fMoleculesFromIRT[MolA].position;
	G4ThreeVector Position2 = fMoleculesFromIRT[MolB].position;
	G4ThreeVector Position = (Position1 + Position2)/2;
	TsIRTConfiguration::TsMolecule aMol;
	aMol.id       = Molecule;
	aMol.position = Position;
	aMol.time     = fCurrentTime*s;
	aMol.reacted  = false;
	aMol.trackID  = 0;
	aMol.isDNA    = false;
	if ( Molecule == 1 | Molecule == 5 )
			aMol.spin = G4UniformRand() > 0.5 ? 1 : 0;
		else
			aMol.spin = -1;

	G4double dT = std::abs(fCurrentTime*s - Time);
	fIRT->GetIRTConfiguration()->Diffuse(aMol,dT);

	fMoleculesFromIRT[fSpeciesIndex] = aMol;
	fSpeciesIndex++;
}


void TsGillespie::AddMolecule(G4int Molecule) {
	TsIRTConfiguration::TsMolecule aMol;
	aMol.id = Molecule;
	aMol.position = G4ThreeVector();
	aMol.time     = fInitialTime;
	aMol.reacted  = false;
	aMol.trackID  = 0;
	aMol.isDNA    = false;
	aMol.isNew    = false;
	aMol.spin     = -1;

	fMoleculesFromIRT[fSpeciesIndex] = aMol;
	fSpeciesIndex++;
}


void TsGillespie::SampleNewMolecules() {
	fMoleculesFromIRT.clear();
	fSpeciesIndex = 0;
	for (auto& IndexAndCount:fConcentrations){
		G4int Index = IndexAndCount.first;
		G4int Count = IndexAndCount.second;
		for(G4int i = 0; i < Count;i++){
			AddMolecule(Index);
		}
	}
}

void TsGillespie::DoReaction(G4double Prop, G4int Index) {
	if (Prop <= 0) {return;}
	G4int index = Index;
    fLastReactionIndex = index;
	G4int tBin = fUtils->FindBin(fCurrentTime, fTimeSteps);

    if (index < 0){
    	G4int delta = fReverse[-index];
    	fConcentrations[-index] += delta;
    	fTotalNumberOfMolecules += delta;

    	fDeltaGPerReactionPerTime[index][fTimeSteps[tBin]]++;
//        G4cout<<"Input: "<<index<<'\t'<<fCurrentTime<<'\t'<<fMoleculeNames[-index]<<" ("<<fConcentrations[-index]<<")"<<G4endl;
        return;
    }

    G4int ReactA = fReactions[index].reactorA;
    G4int ReactB = fReactions[index].reactorB;
    std::vector<G4int> Prod = fReactions[index].products;

//    G4cout<<"At rest: "<<index<<'\t'<<fCurrentTime<<'\t'<<fReactions[index].reactionType<<'\t'<<fMoleculeNames[ReactA]<<" ("<<fConcentrations[ReactA]<<") + "<<fMoleculeNames[ReactB]<<" ("<<fConcentrations[ReactB]<<") => ";

    for (size_t i = 0; i < Prod.size(); i++) {
    	if (Prod[i] == fNoneIndex) {continue;}
        fConcentrations[Prod[i]] ++;
    	fTotalNumberOfMolecules ++;

//    	G4cout<<fMoleculeNames[Prod[i]]<<" ("<<fConcentrations[Prod[i]]<<") "<<'\t';
    }
//    G4cout<<G4endl;
    fConcentrations[ReactA] --;
    fTotalNumberOfMolecules --;

    if (fReactions[index].reactionType != 6){
    	fConcentrations[ReactB]--;
    	fTotalNumberOfMolecules--;
    }

	fDeltaGPerReactionPerTime[index][fTimeSteps[tBin]]++;

}


void TsGillespie::UpdateConcentrations() {
	G4int TotalNumberOfMoleculesInPulse = 0;
	G4double Time = fPulseTimes[0];
	std::map<G4int,G4int> MoleculesAndYields = fEscapeYieldsAtTime[0];
	if (Time/s <= fCurrentTime) {
		for (auto& MolindexAndYields:MoleculesAndYields) {
			G4int MolIndex = MolindexAndYields.first;
			G4int Yields   = MolindexAndYields.second;
			if (Yields == 0 && MolIndex > 100) {continue;}
			fConcentrations[MolIndex] += Yields;
			fTotalNumberOfMolecules   += Yields;
			TotalNumberOfMoleculesInPulse += Yields;
		}
		fEscapeYieldsAtTime.erase(fEscapeYieldsAtTime.begin());
		fPulseTimes.erase(fPulseTimes.begin());
		return;
	}
}


void TsGillespie::RunAlt(G4double iniT, G4double finT) {

	fLastReactionIndex    = 1000000;

	// unit remove
	if (iniT > 0) {
		fCurrentTime = iniT;
		fInitialTime = iniT;
	}
	else {
		fCurrentTime = fTimeSteps[0];
	}

	if (finT > 0) {
		fFinalTime   = finT;
	}
	else {
		fFinalTime = fTimeSteps[fTimeSteps.size() - 1];
	}

/*	if (iniT > 0) {
		fCurrentTime = iniT/s;
		fInitialTime = iniT;
	}
	else {
		fCurrentTime = fTimeSteps[0]/s;
	}

	if (finT > 0) {
		fFinalTime   = finT;
		fFinalTimeInSeconds = finT/s;
	}
	else {
		fFinalTime = fTimeSteps[fTimeSteps.size() - 1];
		fFinalTimeInSeconds = fFinalTime/s;
	}*/
	G4int CurrentPulse = 1;

	fReactions = fIRT->GetIRTConfiguration()->GetReactions();

	RecoverReactionData();
	InitializeContainers();

	G4int tBin = 0;
	G4double timeStep = 0;

	G4int max_time = fTimeSteps.size() - 1;

	while(fCurrentTime <= fFinalTime) {
		auto step = Propensity();
		G4double StepPropensity = step.second;
		G4int StepIndex = step.first;

		if(StepPropensity > 0) timeStep = TimeIncrement(StepPropensity);
		else break;

		fCurrentTime += timeStep;
		DoReaction(StepPropensity, StepIndex);

		while(fCurrentTime > fTimeSteps[tBin] && tBin <= max_time){
			for(auto it = fMoleculeNames.begin(); it != fMoleculeNames.end(); ++it){
				G4String Mol = it->second;
				G4int index = it->first;
				fMoleculesAtTime[Mol][fTimeSteps[tBin]] = fConcentrations[index];
			}

			G4cout<<fTimeSteps[tBin] / s<<" s proceed with "<<fTotalNumberOfMolecules<<" remaining molecules"<<G4endl;
			tBin++;
		}
	}

	while(tBin < fTimeSteps.size()){
		for(auto it = fMoleculeNames.begin(); it != fMoleculeNames.end(); ++it){
			G4String Mol = it->second;
			G4int index = it->first;
			fMoleculesAtTime[Mol][fTimeSteps[tBin]] = fConcentrations[index];
		}

		G4cout<<fTimeSteps[tBin] / s<<" s proceed with "<<fTotalNumberOfMolecules<<" remaining molecules"<<G4endl;
		tBin++;
	}
	fConcentrationsSum1.clear();
	fConcentrationsSum2.clear();
}

void TsGillespie::DiffuseAllMolecules() {
	for(auto& IndexAndMolecule:fMoleculesFromIRT) {
		G4int Index = IndexAndMolecule.first;
		G4double dT = std::abs(fCurrentTime - fMoleculesFromIRT[Index].time);
		fIRT->GetIRTConfiguration()->Diffuse(fMoleculesFromIRT[Index],dT);
	}
}


void TsGillespie::AddEscapeYields(std::map<G4int,G4int> Yields,G4double Time) {

	if (fPulseTimes.size() > 0) {
		G4double LastTime = fPulseTimes[fPulseTimes.size()-1];
		if (Time < LastTime) {
			Time = LastTime + (LastTime/fPulseTimes.size());
		}
	}
	fReactions = fIRT->GetIRTConfiguration()->GetReactions();
	InitializeContainers();
	G4int tBin = fUtils->FindBin(Time,fTimeSteps);
	for (size_t i = tBin; i < fTimeSteps.size(); i++) {
		G4double time = fTimeSteps[i];
		for (auto& IndexAndYields:Yields) {
			G4int Index = IndexAndYields.first;
			G4int Yield = IndexAndYields.second;
			G4String Name = fMoleculeNames[Index];
			if (Yield == 0 && Index > 100) {continue;}
			fMoleculesAtTime[Name][time] += Yield;
		}
	}

	fPulseTimes.push_back(Time);
	fEscapeYieldsAtTime.push_back(Yields);
}

void TsGillespie::SetMolPerTime(G4int index, G4double Gvalue){
	G4double DR = fDoserate;

	if(fMoleculeNames[index] == "None" || fMoleculeNames[index] == "None2") return;

	if(Gvalue < 0) return;
	else if(Gvalue == 0) return;
	else fReverse[index] = 1;

/*
	G4cout<<index<<'\t'<<fMoleculeNames[index]<<'\t'<<"G-value: "<<Gvalue<<" #/100 eV"<<'\t';
	fMolPerTime[index] = abs(Gvalue) / CLHEP::Avogadro / 100 / 1.602e-19 * DR;
	G4cout<<" Molar per time: "<<fReverse[index] * fMolPerTime[index]<<" M/s in "<<DR/gray*s<<" Gy/s"<<'\t';
	fMolPerTime[index] *= fCubicVolume * CLHEP::Avogadro;
	G4cout<<"converted to "<<fReverse[index] * fMolPerTime[index]<<" molecules/s in "<<fCubicVolume/m3<<" m3"<<G4endl;
*/
	G4double Gvalue_number = abs(Gvalue) / (100* eV); // # / MeV
	G4cout<<index<<'\t'<<fMoleculeNames[index]<<'\t'<<"G-value: "<<Gvalue_number * (100 * eV)<<" #/100 eV"<<'\t';

	G4double mass = fCubicVolume * g/cm3; // mm3 / 1e6 dm3/mm3 * kg = kg = J*s*s/m/m
	G4double molPerDose = Gvalue_number / CLHEP::Avogadro / fCubicVolume * mass; // mol/mm3 / (MeV/kg)
	G4cout<<" Molar per gray: "<<molPerDose /(mole/1e6/mm3)*gray<<" M/Gy"<<'\t';

	fMolPerTime[index] = molPerDose * DR;
	G4cout<<" Molar per time: "<<fReverse[index] * fMolPerTime[index] /(mole/1e6/mm3)*s<<" M/s in "<<DR/gray*s<<" Gy/s"<<'\t';

	fMolPerTime[index] *= fCubicVolume * CLHEP::Avogadro;
	G4cout<<"converted to "<<fReverse[index] * fMolPerTime[index] * s<<" molecules/s in "<<fCubicVolume/mm3<<" mm3 "<<mass/kg<<" kg"<<'\t';

	fInsertPropensity += fMolPerTime[index];
	G4cout<<"Total propensity: "<<fInsertPropensity<<G4endl;
}


void TsGillespie::AddDeltaYields(std::map<G4int,G4int> Yields, G4double Time) {
	fReactions = fIRT->GetIRTConfiguration()->GetReactions();
	InitializeContainers();
	//fDeltaYieldsAtTime[Time] = Yields;
	G4int tBin = fUtils->FindBin(Time, fTimeSteps);
	for(size_t i = tBin; i < fTimeSteps.size(); i++) {
		G4double time = fTimeSteps[i];
		for(auto& IndexAndDelta:Yields) {
			G4int Index = IndexAndDelta.first;
			G4int Delta = IndexAndDelta.second;
			//fDeltaYieldsAtTime[time][Index] += Delta;
			fDeltaGPerReactionPerTime[Index][time] += Delta;
		}
	}
}


void TsGillespie::PrintConcentrations() {
	G4cout << "Concentrations: " << fConcentrations.size() << G4endl;
	for (auto& MoleculesAndConcentrations:fConcentrations) {
		G4int ID = MoleculesAndConcentrations.first;
		G4String Name = fMoleculeNames[ID];
		G4double Conc = MoleculesAndConcentrations.second;

		G4cout << Name << ": " << Conc << G4endl;
	}
}

void TsGillespie::PrintMoleculesAtTime() {
	for (auto& NameAndTimesValue:fMoleculesAtTime) {
		G4String Name = NameAndTimesValue.first;
		for (auto& TimesAndValue:NameAndTimesValue.second) {
			G4double Time  = TimesAndValue.first;
			G4double Value = TimesAndValue.second;
			G4cout << Name << "  " << Time/ps << "  " << Value << G4endl;
		}
	}
}


void TsGillespie::PrintReactionInfo(G4int Index) {
	if (Index >= 0) {
		G4String MolA = fMoleculeNames[fReactions[Index].reactorA];
		G4String MolB = fMoleculeNames[fReactions[Index].reactorB];
		G4double Kobs = fReactions[Index].kobs;
		G4cout << Index << " | " << MolA << " + " << MolB << " | " << "kObs = " << Kobs << G4endl;
	}
	else {
		for (auto& IndexAndReaction:fReactions) {
			G4int index   = IndexAndReaction.first;
			G4String MolA = fMoleculeNames[fReactions[index].reactorA];
			G4String MolB = fMoleculeNames[fReactions[index].reactorB];
			G4double Kobs = fReactions[index].kobs;
			G4cout << index << " | " << MolA << " + " << MolB << " | " << "kObs = " << Kobs << G4endl;
		}
	}
}


void TsGillespie::PrintIndividualPropensities() {
	auto Step = Propensity();
	G4double StepPropensity = Step.second;
	G4int StepIndex = Step.first;
//	G4double StepPropensity = Propensity();
	for (auto& IndexAndReaction:fReactions) {
		G4int Index = IndexAndReaction.first;
		G4double ConA = fConcentrations[fReactions[Index].reactorA];
		G4double ConB = fConcentrations[fReactions[Index].reactorB];
		G4int Type    = fReactions[Index].reactionType;
		G4String MolA = fMoleculeNames[fReactions[Index].reactorA];
		G4String MolB = fMoleculeNames[fReactions[Index].reactorB];
		G4double Kobs = fReactions[Index].kobs;
		G4double Prop = Kobs * ConA;
		if (Type!=6) {
			if (ConA == ConB)
				Prop *= ConB-1;
			else
				Prop *= ConB;
		}
		else {
			ConB = fReactions[Index].concentration;
			if (ConB != 0)
				Prop*=ConB;
		} 
		G4cout << Index << " | " << Type << " | " << MolA << " + " << MolB << " | ConA =" <<  ConA << " | ConB =" << ConB << " | Kobs = " << Kobs << " | Prop = " << Prop << "  " << Prop/StepPropensity <<G4endl;
	}
}


void TsGillespie::SaveState() {
	fTotalNumberOfMoleculesSum1 += fTotalNumberOfMolecules;
	fTotalNumberOfMoleculesSum2 += fTotalNumberOfMolecules*fTotalNumberOfMolecules;
	for(auto& IndexAndConcentrations:fConcentrations) {
		G4int Index = IndexAndConcentrations.first;
		G4int Conct = IndexAndConcentrations.second;
		fConcentrationsSum1[Index] += Conct;
		fConcentrationsSum2[Index] += Conct*Conct;
	}
	fNumberOfSavedStates++;
}


