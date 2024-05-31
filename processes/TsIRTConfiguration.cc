// Extra Class for TsIRT
#include "TsIRTConfiguration.hh"
#include "TsIRTUtils.hh"
#include "TsParameterManager.hh"

#include "G4SystemOfUnits.hh"
#include "G4RandomDirection.hh"
#include "G4PhysicalConstants.hh"
#include "Randomize.hh"
#include "G4UIcommand.hh"
#include "G4IosFlagsSaver.hh"

#include <vector>
#include <map>
#include <fstream>
#include <cmath>
#include <algorithm>
#include <functional>
#include <sstream>
#include <random>

TsIRTConfiguration::TsIRTConfiguration(G4String name, TsParameterManager* pM)
: fPm(pM), fName(name), fReactionID(0), fTotalBinaryReaction(0), fLastMoleculeID(0),
fKick(false), fAllTotallyDiffusionControlled(false)
{
	G4String chemistryList = fPm->GetStringParameter("Ch/ChemistryName");
	
	fUtils = new TsIRTUtils();
	fLowerTime = 1.0e-13*s;
	fUpperTime = 1.0e-6*s;
	
	fExistingMolecules.clear();
	
	fExistingMolecules["hydrogen"]         = "H^0";
	fExistingMolecules["hydroxyl"]         = "OH^0";
	fExistingMolecules["hydrogenperoxide"] = "H2O2^0";
	fExistingMolecules["dyhydrogen"]       = "H_2^0";
	fExistingMolecules["solvatedelectron"] = "e_aq^-1";
	fExistingMolecules["hydronium"]        = "H3O^1";
	fExistingMolecules["hydroxide"]        = "OH^-1";
	fExistingMolecules["oxygen"]           = "O2^0";
	fExistingMolecules["superoxideanion"]  = "O2^-1";
	fExistingMolecules["hydroperoxy"]      = "HO2^0";
	fExistingMolecules["dioxidanide"]      = "HO2^-1";
	fExistingMolecules["atomicoxygen"]     = "O3P^0";
	fExistingMolecules["oxyde"]            = "O^-1";
	fExistingMolecules["trioxide"]         = "O3^-1";
	fExistingMolecules["ozone"]            = "O3^0";
	fExistingMolecules["none"]             = "None";
	
	AddMolecule("H^0",     7.0e9*nm*nm/s,   0, 0.19*nm);
	AddMolecule("OH^0",    2.2e9*nm*nm/s,   0, 0.22*nm);
	AddMolecule("H2O2^0",  2.3e9*nm*nm/s,   0, 0.21*nm);
	AddMolecule("H_2^0",   4.8e9*nm*nm/s,   0, 0.14*nm);
	AddMolecule("e_aq^-1", 4.9e9*nm*nm/s,  -1, 0.50*nm);
	AddMolecule("H3O^1",   9.46e9*nm*nm/s,  1, 0.25*nm);
	AddMolecule("OH^-1",   5.3e9*nm*nm/s,  -1, 0.33*nm);
	AddMolecule("O2^0",    2.4e9*nm*nm/s,   0, 0.17*nm);
	AddMolecule("O2^-1",   1.75e9*nm*nm/s, -1, 0.22*nm);
	AddMolecule("HO2^0",   2.3e9*nm*nm/s,  0, 0.21*nm);
	AddMolecule("HO2^-1",  1.4e9*nm*nm/s, -1, 0.25*nm);
	AddMolecule("O3P^0",   2.0e9*nm*nm/s,  0, 0.20*nm);
	AddMolecule("O^-1",    2.0e9*nm*nm/s, -1, 0.25*nm);
	AddMolecule("O3^-1",   2.0e9*nm*nm/s, -1, 0.20*nm);
	AddMolecule("O3^0",    2.0e9*nm*nm/s,  0, 0.20*nm);
	AddMolecule("None",    0.0e9*nm*nm/s,  0, 0.00*nm);
	
	// Re-set diffusion coefficient, radius and/or charge.
	std::vector<G4String>* moleculeNames = new std::vector<G4String>;
	fPm->GetParameterNamesStartingWith("Mo/", moleculeNames);
	G4int numberOfMolecules = moleculeNames->size();
	std::vector<G4String> moleculesDontExist;
	
	for ( int i = 0; i < numberOfMolecules; i++ ) {
		G4String fullName = (*moleculeNames)[i];
		G4StrUtil::to_lower(fullName);
		G4bool moleculeExists = false;
		
		if ( G4StrUtil::contains(fullName, "diffusioncoefficient") ) {
			G4String molName = fullName.substr(3, fullName.find("diffusioncoefficient")-4);
			if ( fExistingMolecules.find(molName) != fExistingMolecules.end() ) {
				G4int molID =  fMoleculesID[fExistingMolecules[molName]];
				fMoleculesDefinition[molID].diffusionCoefficient = fPm->GetDoubleParameter(fullName,"surface perTime");
				G4cout << molName << ": Re-set diffusion coefficient to "
				<< fMoleculesDefinition[molID].diffusionCoefficient/(nm*nm/s) << " nm2/s" << std::endl;
				moleculeExists = true;
			}
		}
		
		if (G4StrUtil::contains(fullName, "radius") ) {
			G4String molName = fullName.substr(3, fullName.find("radius")-4);
			if ( fExistingMolecules.find(molName) != fExistingMolecules.end() ) {
				G4int molID =  fMoleculesID[fExistingMolecules[molName]];
				fMoleculesDefinition[molID].radius = fPm->GetDoubleParameter(fullName, "Length");
				G4cout << molName << ": Re-set radius to "
				<< fMoleculesDefinition[molID].radius/nm << " nm" << std::endl;
				moleculeExists = true;
			}
		}
		
		if (G4StrUtil::contains(fullName, "charge")) {
			G4String molName = fullName.substr(3, fullName.find("charge")-4);
			if ( fExistingMolecules.find(molName) != fExistingMolecules.end() ) {
				G4int molID =  fMoleculesID[fExistingMolecules[molName]];
				fMoleculesDefinition[molID].charge = fPm->GetUnitlessParameter(fullName);
				G4cout << molName << ": Re-set charge to "
				<< fMoleculesDefinition[molID].radius/nm << " nm" << std::endl;
				moleculeExists = true;
			}
		}
		
		if (G4StrUtil::contains(fullName, "symbol") && !moleculeExists) {
			moleculesDontExist.push_back(fullName.substr(3, fullName.find("symbol")-4));
		}
	}
	
	// Creates user-defined molecules.
	for ( size_t u = 0; u < moleculesDontExist.size(); u++ ) {
		G4StrUtil::to_lower(moleculesDontExist[u]);
		G4double charge = fPm->GetUnitlessParameter("Mo/" + moleculesDontExist[u] + "/Charge");
		G4double radius = fPm->GetDoubleParameter("Mo/" + moleculesDontExist[u] + "/Radius", "Length");
		G4double diffusionCoefficient = fPm->GetDoubleParameter("Mo/" + moleculesDontExist[u] +
																"/DiffusionCoefficient","surface perTime");
		G4String symbol = fPm->GetStringParameter("Mo/" + moleculesDontExist[u] + "/Symbol");
		
		if ( fPm->ParameterExists("Mo/" + moleculesDontExist[u] + "/AssignMoleculeID"))
			AddMolecule(symbol, fPm->GetIntegerParameter("Mo/" + moleculesDontExist[u] + "/AssignMoleculeID"),
						diffusionCoefficient, charge, radius);
		else
			AddMolecule(symbol, diffusionCoefficient, charge, radius);
		
		fExistingMolecules[moleculesDontExist[u]] = symbol;
	}
	
	G4String parName = "Ch/" + chemistryList + "/SetAllReactionsTotallyDiffusionControlled";
	if ( fPm->ParameterExists(parName) )
		fAllTotallyDiffusionControlled = fPm->GetBooleanParameter(parName);
	
	// Declare And Insert Binary Reactions
	std::vector<G4String>* reactionNames = new std::vector<G4String>;
	fPm->GetParameterNamesBracketedBy("Ch/" , "Products", reactionNames);
	G4int numberOfReactions = reactionNames->size();
	G4int prefixLength = G4String("Ch/" + chemistryList + "/Reaction/").length();
	
	for ( int i = 0; i < numberOfReactions; i++ ) {
		G4String aparName = (*reactionNames)[i];
		if ( G4StrUtil::contains(aparName, "BackgroundReaction"))
			continue;
		
		if ( fPm->ParameterExists(aparName.substr(0,aparName.find("Products")-1) + "/Active") &&
			!fPm->GetBooleanParameter(aparName.substr(0,aparName.find("Products")-1) + "/Active") )
			continue;
		
		G4String reactions = aparName.substr(prefixLength, aparName.find("Products")-prefixLength-1);
		G4String reactorA = reactions.substr(0, reactions.find("/"));
		G4String reactorB = reactions.substr(reactions.find("/") + 1);
		G4StrUtil::to_lower(reactorA);
		G4StrUtil::to_lower(reactorB);
		G4String* product = fPm->GetStringVector(aparName);
		G4int nbOfProduct = fPm->GetVectorLength(aparName);
		std::vector<G4String> vProduct;
		
		for ( int j = 0; j < nbOfProduct; j++ ) {
			G4StrUtil::to_lower(product[j]);
			vProduct.push_back(product[j]);
		}
		
		G4StrUtil::to_lower(reactorA);
		G4StrUtil::to_lower(reactorB);
		
		G4int reactionType = fPm->GetIntegerParameter(aparName.substr(0,aparName.find("Products")-1) + "/ReactionType");

		G4double reactionRate = fPm->GetDoubleParameter(aparName.substr(0,aparName.find("Products")-1) +
														"/ReactionRate","perMolarConcentration perTime");
		
		if ( nbOfProduct == 1 )
			InsertReaction(reactorA, reactorB,
						   product[0],"None","None",reactionRate,reactionType);
		else if ( nbOfProduct == 2 )
			InsertReaction(reactorA, reactorB,product[0],product[1],"None",reactionRate,reactionType);
		else if (nbOfProduct == 3 )
			InsertReaction(reactorA, reactorB,product[0],product[1],product[2],reactionRate,reactionType);
		else 
			InsertReaction(reactorA, reactorB,vProduct,reactionRate,reactionType);

		if (fPm->ParameterExists(aparName.substr(0,aparName.find("Products")-1) + "/ActivationRate") ||
			fPm->ParameterExists(aparName.substr(0,aparName.find("Products")-1) + "/DiffusionRate")){
			if (fReactions[fReactionID-1].reactionType == 2 || fReactions[fReactionID-1].reactionType == 4) {
				G4double activationRate = fPm->GetDoubleParameter(aparName.substr(0,aparName.find("Products")-1) +
															"/ActivationRate","perMolarConcentration perTime");
				G4double diffusionRate  = fPm->GetDoubleParameter(aparName.substr(0,aparName.find("Products")-1) +
															"/DiffusionRate","perMolarConcentration perTime");
				fReactions[fReactionID-1].kdif = diffusionRate;
				fReactions[fReactionID-1].kact = activationRate;
				G4cout << "Advance Reaction Mode for reaction" << fReactionID-1 << " : " 
				       << fExistingMolecules[reactorA] << " + " << fExistingMolecules[reactorB] << G4endl;
			}
			else {
				G4cout << "WARNING: Advance Reaction Options Only available for reations type II or IV" << G4endl;
			}
		}
	}
	
	// Declare And Insert Background Reactions
	prefixLength = G4String("Ch/" + chemistryList + "/BackgroundReaction/").length();
	for ( int i = 0; i < numberOfReactions; i++ ) {
		G4String aparName = (*reactionNames)[i];
		if ( G4StrUtil::contains(aparName, "/Reaction/"))
			continue;
		
		if ( fPm->ParameterExists(aparName.substr(0,aparName.find("Products")-1) + "/Active") &&
			!fPm->GetBooleanParameter(aparName.substr(0,aparName.find("Products")-1) + "/Active") )
			continue;
		
		G4String reactions = aparName.substr(prefixLength, aparName.find("Products")-prefixLength-1);
		G4String reactorA = reactions.substr(0, reactions.find("/"));
		G4String reactorB = reactions.substr(reactions.find("/") + 1);
		G4StrUtil::to_lower(reactorA);
		G4StrUtil::to_lower(reactorB);
		G4String* product = fPm->GetStringVector(aparName);
		G4int nbOfProduct = fPm->GetVectorLength(aparName);
		std::vector<G4String> vProduct;
		
		for ( int j = 0; j < nbOfProduct; j++ ) {
			G4StrUtil::to_lower(product[j]);
			vProduct.push_back(product[j]);

		}
		
		G4StrUtil::to_lower(reactorA);
		G4StrUtil::to_lower(reactorB);
		
		if ( fPm->ParameterExists(aparName.substr(0,aparName.find("Products")-1) + "/ScavengingCapacity") ) {
			G4double scavengingCapacity = fPm->GetDoubleParameter(aparName.substr(0,aparName.find("Products")-1) +
																  "/ScavengingCapacity","perTime");
			if ( nbOfProduct == 1 )
				InsertBackgroundReaction(reactorA, reactorB,product[0],"None","None",scavengingCapacity,false);
			else if ( nbOfProduct == 2 )
				InsertBackgroundReaction(reactorA, reactorB,product[0],product[1],"None",scavengingCapacity,false);
			else if (nbOfProduct == 3)
				InsertBackgroundReaction(reactorA, reactorB,product[0],product[1],product[2],scavengingCapacity,false);
			else
				InsertBackgroundReaction(reactorA, reactorB,vProduct,scavengingCapacity,false);
		} else {
			G4double concentration = fPm->GetDoubleParameter(aparName.substr(0,aparName.find("Products")-1) +
															 "/Concentration","molar concentration");
			G4double reactionRate = fPm->GetDoubleParameter(aparName.substr(0,aparName.find("Products")-1) +
															"/ReactionRate","perMolarConcentration perTime");
			
			G4String scavengingExponentialModel = "exponentialsinglefactor";
			if ( fPm->ParameterExists(aparName.substr(0,aparName.find("Products")-1) + "/ScavengingModel"))
				scavengingExponentialModel = fPm->GetStringParameter(aparName.substr(0,aparName.find("Products")-1) + "/ScavengingModel");
			
			G4StrUtil::to_lower(scavengingExponentialModel);
			if ( scavengingExponentialModel == "exponentialsinglefactor" ) {
				if ( nbOfProduct == 1 )
					InsertBackgroundReaction(reactorA, reactorB,product[0],"None","None",reactionRate,concentration,false);
				else if ( nbOfProduct == 2 )
					InsertBackgroundReaction(reactorA, reactorB,product[0],product[1],"None",reactionRate,concentration,false);
				else if (nbOfProduct == 3)
					InsertBackgroundReaction(reactorA, reactorB,product[0],product[1],product[2],reactionRate,concentration,false);
				else
					InsertBackgroundReaction(reactorA, reactorB,vProduct,reactionRate,concentration,false);
				
			} else if (scavengingExponentialModel == "exponentialdoublefactor" ) {
				if ( nbOfProduct == 1 )
					InsertBackgroundReaction(reactorA, reactorB,product[0],"None","None",reactionRate,concentration,true);
				else if ( nbOfProduct == 2 )
					InsertBackgroundReaction(reactorA, reactorB,product[0],product[1],"None",reactionRate,concentration,true);
				else if (nbOfProduct == 3)
					InsertBackgroundReaction(reactorA, reactorB,product[0],product[1],product[2], reactionRate,concentration,true);
				else
					InsertBackgroundReaction(reactorA, reactorB,vProduct, reactionRate,concentration,true);
			} else {
				Quit(aparName.substr(0,aparName.find("Products")-1) + "/ScavengingModel",
					 "Scavenging model does not exists in TOPAS-nBio database\n Use ExponentialSingleFactor or ExponentialDoubleFactor");
			}
		}
	}
	
	// Reactions for quality check of the IRT loop.
	parName = "Ch/" + chemistryList + "/TestIRTForQualityAssurance";
	fQualityAssurance = false;
	if (fPm->ParameterExists(parName) && fPm->GetBooleanParameter(parName) ) {
		fQualityAssurance = true;
		G4double kUnit = fPm->GetUnitValue("/M/s");
		// Fake molecules for quality assurance.
		AddMolecule("A",       17, 5.0e9*nm*nm/s,  0, 1.00*nm);
		AddMolecule("B",       18, 5.0e9*nm*nm/s,  0, 1.00*nm);
		AddMolecule("AA",      19, 5.0e9*nm*nm/s,  0, 1.00*nm);
		AddMolecule("AB",      20, 5.0e9*nm*nm/s,  0, 1.00*nm);
		AddMolecule("BB",      21, 5.0e9*nm*nm/s,  0, 1.00*nm);
		AddMolecule("AAB",     22, 5.0e9*nm*nm/s,  0, 1.00*nm);
		AddMolecule("ABB",     23, 5.0e9*nm*nm/s,  0, 1.00*nm);
		// Reactions for quality assurance.
		InsertReaction("A", "A", "AA", "None", "None", 5.0e9 * kUnit, 1);
		InsertReaction("A", "B", "AB", "None", "None", 5.0e9 * kUnit, 1);
		InsertReaction("B", "B", "BB", "None", "None", 5.0e9 * kUnit, 1);
		InsertReaction("AB", "A", "AAB", "None", "None", 5.0e9 * kUnit, 1);
		InsertReaction("AB", "B", "ABB", "None", "None", 5.0e9 * kUnit, 1);
	}
	
	// Re-scale chemistry parameters based on temperature
	parName = "Ch/" + chemistryList + "/Temperature";
	fScaleForTemperature = false;
	if ( fPm->ParameterExists(parName) ) {
		fTemperature = fPm->GetUnitlessParameter(parName);
		fScaleForTemperature = true;
		
		parName = "Ch/" + chemistryList + "/ApplyCorrectionScalingForTemperature";
		fKick = false;
		if ( fPm->ParameterExists(parName) )
			fKick = fPm->GetBooleanParameter(parName);
		
		AdjustReactionAndDiffusionRateForTemperature();
	}
	
	ResolveReactionRateCoefficients();
	CalculateContactProbabilities();
	ResolveRemainerReactionParameters();
	
	if ( fPm->ParameterExists("Ch/"+chemistryList+"/ModelAcidPropertiesFromSubstance") ) {
		fpHSolventConcentration = 0.0;
		fpHValue                = 7.1;
		
		fpHSolvent = fPm->GetStringParameter("Ch/"+chemistryList+"/ModelAcidPropertiesFromSubstance");
		G4StrUtil::to_lower(fpHSolvent);
		if (fPm->ParameterExists("Ch/"+chemistryList+"/ModelAcidPropertiesWithConcentration") &&
			fPm->ParameterExists("Ch/"+chemistryList+"/ModelAcidPropertiesWithpH")) {
			G4String message = "Cannot be defined when parameter: Ch/" + chemistryList + "/ModelAcidPropertiesWithConcentration is used.";
			Quit("Ch/"+chemistryList+"/ModelAcidPropertiesWithpH",message);
		}
		
		if (fPm->ParameterExists("Ch/"+chemistryList+"/ModelAcidPropertiesWithConcentration") ) {
			fpHSolventConcentration = fPm->GetDoubleParameter("Ch/"+chemistryList+"/ModelAcidPropertiesWithConcentration","molar concentration");
			AdjustReactionRateForPH("Concentration");
		}
		
		if (fPm->ParameterExists("Ch/"+chemistryList+"/ModelAcidPropertiesWithpH") ) {
			fpHValue = fPm->GetUnitlessParameter("Ch/"+chemistryList+"/ModelAcidPropertiesWithpH");
			AdjustReactionRateForPH("PH");
		}
	}

	PrintReactionsInformation();
}


TsIRTConfiguration::~TsIRTConfiguration()
{
	delete fUtils;
}


void TsIRTConfiguration::AddMolecule(G4String name, G4int moleculeID, G4double diffusionCoefficient,
									 G4double charge, G4double radius) {
	TsMoleculeDefinition aMolecule;
	aMolecule.diffusionCoefficient = diffusionCoefficient;
	aMolecule.charge = charge;
	aMolecule.radius = radius;
	
	fMoleculesDefinition[moleculeID] = aMolecule;
	fMoleculesID[name] = moleculeID;
	fMoleculesName[moleculeID] = name;
	
	fLastMoleculeID = moleculeID;
	
	G4bool found = false;
	for (auto it = fExistingMolecules.begin(); it != fExistingMolecules.end(); ++it) {
		if (it->second == name) {
			found = true;
			break;
		}
	}
	
	if ( !found )
		fExistingMolecules[name] = name;
	
}


void TsIRTConfiguration::AddMolecule(G4String name, G4double diffusionCoefficient,
									 G4double charge, G4double radius) {
	fLastMoleculeID++;
	TsMoleculeDefinition aMolecule;
	aMolecule.diffusionCoefficient = diffusionCoefficient;
	aMolecule.charge = charge;
	aMolecule.radius = radius;
	
	fMoleculesDefinition[fLastMoleculeID] = aMolecule;
	fMoleculesID[name] = fLastMoleculeID;
	fMoleculesName[fLastMoleculeID] = name;
	
	G4bool found = false;
	for (auto it = fExistingMolecules.begin(); it != fExistingMolecules.end(); ++it) {
		if (it->second == name) {
			found = true;
			break;
		}
	}
	
	if ( !found )
		fExistingMolecules[name] = name;
	
}


void TsIRTConfiguration::AddMolecule(G4String name) {
	fLastMoleculeID++;
	AddMolecule(name, fLastMoleculeID, 0.0, 0.0, 0.0);
	return;
}


G4bool TsIRTConfiguration::MoleculeExists(G4String name) {
	if ( fMoleculesID.find(name) == fMoleculesID.end() )
		return false;
	return true;
}


void TsIRTConfiguration::QuitIfMoleculeNotFound(G4String mol) {
	if (fExistingMolecules.find(mol) == fExistingMolecules.end()) {
		G4cerr << "TOPAS is exiting due to a fatal error in IRT Chemistry setup!" << G4endl;
		G4cerr << "--- Molecule " << mol << " does not exists in database" << G4endl;
		fPm->AbortSession(1);
	}
	return;
}


G4double TsIRTConfiguration::GetMoleculeRadius(G4int moleculeID) {
	return fMoleculesDefinition[moleculeID].radius;
}


G4int TsIRTConfiguration::GetMoleculeCharge(G4int moleculeID) {
	return fMoleculesDefinition[moleculeID].charge;
}


TsIRTConfiguration::TsMolecularReaction TsIRTConfiguration::GetReaction(G4int index) {
	if (fReactions.count(index) == 1)
		return fReactions[index];
	else {
		TsMolecularReaction NoReaction;
		NoReaction.index = -1;
		return NoReaction;
	}
}


void TsIRTConfiguration::Diffuse(TsMolecule& mol, G4double dt) {
	G4double sigma, x, y, z;
	G4double diffusionCoefficient = fMoleculesDefinition[mol.id].diffusionCoefficient;
	
	sigma = std::sqrt(2.0 * diffusionCoefficient * dt);
	
	x = G4RandGauss::shoot(0., 1.0)*sigma;
	y = G4RandGauss::shoot(0., 1.0)*sigma;
	z = G4RandGauss::shoot(0., 1.0)*sigma;
	mol.position += G4ThreeVector(x, y, z);
	mol.time += dt;
}


void TsIRTConfiguration::PrintMoleculesInformation() {
	for ( auto& molecules : fMoleculesID ) {
		G4String name = molecules.first;
		G4int id = molecules.second;
		TsMoleculeDefinition aMolecule = fMoleculesDefinition[id];
		std::cout << "----" << std::endl;
		std::cout << "Molecule name           : " << name << std::endl;
		std::cout << "  Molecule ID           : " << id << std::endl;
		std::cout << "  Diffusion coefficient : " << aMolecule.diffusionCoefficient/(nm*nm/s)
		<< " nm2/s " << std::endl;
		std::cout << "  Charge                : " << aMolecule.charge << " e " << std::endl;
		std::cout << "  Radius                : " << aMolecule.radius/nm << " nm " << std::endl;
		std::cout << "----" << std::endl;
		std::cout << "" << std::endl;
	}
}


void TsIRTConfiguration::PrintReactionsInformation() {
	G4IosFlagsSaver iosfs(G4cout);
	std::map<size_t, std::vector<TsMolecularReaction> > temporal;
	for ( size_t i = 0; i < fReactions.size(); i++) {
		G4int type = fReactions[i].reactionType;
		temporal[type].push_back(fReactions[i]);
	}
	
	G4String* outputReactionParentA = new G4String[fReactions.size()];
	G4String* outputReactionParentB = new G4String[fReactions.size()];
	G4String* outputReactionProducts = new G4String[fReactions.size()];
	G4String* outputKobs = new G4String[fReactions.size()];
	G4String* outputKdif = new G4String[fReactions.size()];
	G4String* outputKact = new G4String[fReactions.size()];
	G4String* outputScav = new G4String[fReactions.size()];
	G4String* outputAlpha = new G4String[fReactions.size()];
	G4String* outputReactionRadius = new G4String[fReactions.size()];
	G4String* outputReactionRadiusEff = new G4String[fReactions.size()];
	G4String* outputReactionRadiusEff1 = new G4String[fReactions.size()];
	G4String* outputProbability = new G4String[fReactions.size()];
	G4String* outputModel = new G4String[fReactions.size()];
	
	G4cout << G4endl;
	G4int n = 0;
	for ( size_t i = 1; i <= temporal.size(); i++ ) {
		for ( auto& reactions : temporal[i] ) {
			G4int type = reactions.reactionType;
			G4int molA = reactions.reactorA;
			G4int molB = reactions.reactorB;
			std::vector<G4int> products = reactions.products;
			
			outputReactionParentA[n] = fMoleculesName[molA];
			outputReactionParentB[n] = fMoleculesName[molB];
			
			if ( products.size() > 0 ) {
				outputReactionProducts[n] = fMoleculesName[products[0]];
				for ( size_t u = 1; u < products.size(); u++ ) {
					outputReactionProducts[n] += " + " + fMoleculesName[products[u]];
				}
			} else {
				outputReactionProducts[n] = "None";
			}
			
			outputKobs[n] = G4UIcommand::ConvertToString(reactions.kobs/fPm->GetUnitValue("/M/s"));
			outputKdif[n] = G4UIcommand::ConvertToString(reactions.kdif/fPm->GetUnitValue("/M/s"));
			outputKact[n] = G4UIcommand::ConvertToString(reactions.kact/fPm->GetUnitValue("/M/s"));
			outputAlpha[n] = G4UIcommand::ConvertToString(reactions.alpha * nm);
			outputReactionRadius[n] = G4UIcommand::ConvertToString(reactions.reactionRadius/nm);
			outputReactionRadiusEff[n] = G4UIcommand::ConvertToString(reactions.effectiveReactionRadius/nm);
			outputReactionRadiusEff1[n] = G4UIcommand::ConvertToString(reactions.effectiveTildeReactionRadius/nm);
			outputProbability[n] = G4UIcommand::ConvertToString(reactions.probabilityOfReaction);
			if (type == 6 ) {
				outputScav[n] = G4UIcommand::ConvertToString(reactions.scavengingCapacity * s);
				outputModel[n] = reactions.sampleExponential ? "true" : "false";
			} else {
				outputScav[n] = "0.0";
				outputModel[n] = "none";
			}
			n++;
		}
	}
	
	G4int maxLengthOutputReactionParentA = -1;
	G4int maxLengthOutputReactionParentB = -1;
	G4int maxLengthOutputReactionProducts = -1;
	G4int maxLengthOutputKobs = -1;
	G4int maxLengthOutputKdif = -1;
	G4int maxLengthOutputKact = -1;
	G4int maxLengthOutputProbability = -1;
	G4int maxLengthOutputReactionRadius = -1;
	G4int maxLengthOutputReactionRadiusEff = -1;
	G4int maxLengthOutputReactionRadiusEff1 = -1;
	G4int maxLengthOutputAlpha = -1;
	
	for ( int i = 0; i < n; i++ ) {
		if ( maxLengthOutputReactionParentA < (G4int)outputReactionParentA[i].length())
			maxLengthOutputReactionParentA = outputReactionParentA[i].length();
		
		if ( maxLengthOutputReactionParentB < (G4int)outputReactionParentB[i].length())
			maxLengthOutputReactionParentB = outputReactionParentB[i].length();
		
		if ( maxLengthOutputReactionProducts < (G4int)outputReactionProducts[i].length())
			maxLengthOutputReactionProducts = outputReactionProducts[i].length();
		
		if ( maxLengthOutputKobs < (G4int)outputKobs[i].length())
			maxLengthOutputKobs = outputKobs[i].length();
		
		if ( maxLengthOutputKdif < (G4int)outputKdif[i].length())
			maxLengthOutputKdif = outputKdif[i].length();
		
		if ( maxLengthOutputKact < (G4int)outputKact[i].length())
			maxLengthOutputKact = outputKact[i].length();
		
		if ( maxLengthOutputAlpha < (G4int)outputAlpha[i].length())
			maxLengthOutputAlpha = outputAlpha[i].length();
		
		if ( maxLengthOutputProbability < (G4int) outputProbability[i].length())
			maxLengthOutputProbability = outputProbability[i].length();
		
		if ( maxLengthOutputReactionRadius < (G4int) outputReactionRadius[i].length())
			maxLengthOutputReactionRadius = outputReactionRadius[i].length();
		
		if ( maxLengthOutputReactionRadiusEff < (G4int) outputReactionRadiusEff[i].length())
			maxLengthOutputReactionRadiusEff = outputReactionRadiusEff[i].length();
		
		if ( maxLengthOutputReactionRadiusEff1 < (G4int) outputReactionRadiusEff1[i].length())
			maxLengthOutputReactionRadiusEff1 = outputReactionRadiusEff1[i].length();
	}
	
	/*maxLengthOutputReactionParentA += 2;
	 maxLengthOutputReactionParentB += 2;*/
	maxLengthOutputReactionProducts += 2;
	maxLengthOutputKobs += 4;
	maxLengthOutputKdif += 4;
	maxLengthOutputKact += 4;
	maxLengthOutputAlpha += 2;
	 maxLengthOutputProbability += 2;
	 maxLengthOutputReactionRadius += 2;
	 maxLengthOutputReactionRadiusEff += 2;
	 maxLengthOutputReactionRadiusEff1 += 2;
	
	G4String* title = new G4String[11];
	title[0] = "Type";
	title[1] = "Molecules";
	title[2] = "Products";
	title[3] = "Kobs (/M/s)";
	title[4] = "Kdif (/M/s)";
	title[5] = "Kact (/M/s)";
	title[6] = "r (nm)";
	title[7] = "Preact";
	title[8] = "alpha (/nm)";
	title[9] = "reff (nm)";
	title[10] = "reff1 (nm)";
	n = 0;
	for ( size_t i = 1; i <= temporal.size(); i++ ) {
		for (int k = 0; k < 120; k++ )
			G4cout << "=";
		G4cout << G4endl;
		if ( i == 1 )
			G4cout << std::setfill(' ') << std::setw(7) << title[0]
			<< std::setw(maxLengthOutputReactionParentA+maxLengthOutputReactionParentB+3) << std::left << title[1]
			<< std::setw(maxLengthOutputReactionProducts+7) << std::left << title[2]
			<< std::setw(maxLengthOutputKobs) << std::left << title[3]
			<< std::setw(maxLengthOutputReactionRadius) << std::left << title[6]
			<< std::setw(maxLengthOutputProbability) << std::left << title[7] << G4endl;
		else if ( i == 2 )
			G4cout << std::setfill(' ') << std::setw(7) << title[0]
			<< std::setw(maxLengthOutputReactionParentA+maxLengthOutputReactionParentB+3) << std::left << title[1]
			<< std::setw(maxLengthOutputReactionProducts+7) << std::left << title[2]
			<< std::setw(maxLengthOutputKobs) << std::left << title[3]
			<< std::setw(maxLengthOutputKdif) << std::left << title[4]
			<< std::setw(maxLengthOutputKact) << std::left << title[5]
			<< std::setw(maxLengthOutputReactionRadius) << std::left << title[6]
			<< std::setw(maxLengthOutputProbability) << std::left << title[7]
			<< std::setw(maxLengthOutputAlpha) << std::left << title[8] << G4endl;
		else if ( i == 3 )
			G4cout << std::setfill(' ') << std::setw(7) << title[0]
			<< std::setw(maxLengthOutputReactionParentA+maxLengthOutputReactionParentB+3) << std::left << title[1]
			<< std::setw(maxLengthOutputReactionProducts+7) << std::left << title[2]
			<< std::setw(maxLengthOutputKobs) << std::left << title[3]
			<< std::setw(maxLengthOutputReactionRadius) << std::left << title[6]
			<< std::setw(maxLengthOutputReactionRadiusEff) << std::left << title[9]
			<< std::setw(maxLengthOutputProbability) << std::left << title[7] << G4endl;
		else if ( i == 4 )
			G4cout << std::setfill(' ') << std::setw(7) << title[0]
			<< std::setw(maxLengthOutputReactionParentA+maxLengthOutputReactionParentB+3) << std::left << title[1]
			<< std::setw(maxLengthOutputReactionProducts+7) << std::left << title[2]
			<< std::setw(maxLengthOutputKobs) << std::left << title[3]
			<< std::setw(maxLengthOutputKdif) << std::left << title[4]
			<< std::setw(maxLengthOutputKact) << std::left << title[5]
			<< std::setw(maxLengthOutputReactionRadius) << std::left << title[6]
			<< std::setw(maxLengthOutputReactionRadiusEff) << std::left << title[9]
			<< std::setw(maxLengthOutputReactionRadiusEff1) << std::left << title[10]
			<< std::setw(maxLengthOutputProbability) << std::left << title[7] << G4endl;
		else if ( i == 5 )
			G4cout << std::setfill(' ') << std::setw(7) << title[0]
			<< std::setw(maxLengthOutputReactionParentA+maxLengthOutputReactionParentB+3) << std::left << title[1]
			<< std::setw(maxLengthOutputReactionProducts+7) << std::left << title[2]
			<< std::setw(maxLengthOutputKobs) << std::left << title[3]
			<< std::setw(maxLengthOutputReactionRadius) << std::left << title[6]
			<< std::setw(maxLengthOutputReactionRadiusEff) << std::left << title[9]
			<< std::setw(maxLengthOutputProbability) << std::left << title[7] << G4endl;
		else if ( i == 6 )
			G4cout << std::setfill(' ') << std::setw(7) << title[0]
			<< std::setw(maxLengthOutputReactionParentA+maxLengthOutputReactionParentB+3) << std::left << title[1]
			<< std::setw(maxLengthOutputReactionProducts+7) << std::left << title[2]
			<< std::setw(maxLengthOutputKobs) << std::left << "Scavenging capacity /s" << G4endl;
		
		for (int k = 0; k < 120; k++ )
			G4cout << "=";
		G4cout << G4endl;
		
		for ( auto& reactions : temporal[i] ) {
			if ( i == 1 ) {
				G4cout << std::setw(7) << " I "
				<< std::setw(maxLengthOutputReactionParentA) << std::left << outputReactionParentA[n] << " + "
				<< std::setw(maxLengthOutputReactionParentB) << std::left << outputReactionParentB[n] << " -> "
				<< std::setw(maxLengthOutputReactionProducts) << std::left << outputReactionProducts[n]
				<< std::setw(maxLengthOutputKobs) << std::left << outputKobs[n]
				<< std::setw(maxLengthOutputReactionRadius) << std::left << outputReactionRadius[n]
				<< std::setw(maxLengthOutputProbability) << std::left << outputProbability[n] <<
				G4endl;
			} else if ( i == 2 ) {
				G4cout << std::setw(7) << " II "
				<< std::setw(maxLengthOutputReactionParentA) << std::left << outputReactionParentA[n] << " + "
				<< std::setw(maxLengthOutputReactionParentB) << std::left << outputReactionParentB[n] << " -> "
				<< std::setw(maxLengthOutputReactionProducts) << std::left << outputReactionProducts[n]
				<< std::setw(maxLengthOutputKobs) << std::left << outputKobs[n]
				<< std::setw(maxLengthOutputKdif) << std::left << outputKdif[n]
				<< std::setw(maxLengthOutputKact) << std::left << outputKact[n]
				<< std::setw(maxLengthOutputReactionRadius) << std::left << outputReactionRadius[n]
				<< std::setw(maxLengthOutputProbability) << std::left << outputProbability[n]
				<< std::setw(maxLengthOutputAlpha) << std::left << outputAlpha[n] <<
				G4endl;
			} else if ( i == 3 ) {
				G4cout << std::setw(7) << " III "
				<< std::setw(maxLengthOutputReactionParentA) << std::left << outputReactionParentA[n] << " + "
				<< std::setw(maxLengthOutputReactionParentB) << std::left << outputReactionParentB[n] << " -> "
				<< std::setw(maxLengthOutputReactionProducts) << std::left << outputReactionProducts[n]
				<< std::setw(maxLengthOutputKobs) << std::left << outputKobs[n]
				<< std::setw(maxLengthOutputReactionRadius+3) << std::left << outputReactionRadius[n]
				<< std::setw(maxLengthOutputReactionRadiusEff+3) << std::left << outputReactionRadiusEff[n]
				<< std::setw(maxLengthOutputProbability+3) << std::left << outputProbability[n] <<
				G4endl;
			} else if ( i == 4 ) {
				G4cout << std::setw(7) << " IV "
				<< std::setw(maxLengthOutputReactionParentA) << std::left << outputReactionParentA[n] << " + "
				<< std::setw(maxLengthOutputReactionParentB) << std::left << outputReactionParentB[n] << " -> "
				<< std::setw(maxLengthOutputReactionProducts) << std::left << outputReactionProducts[n]
				<< std::setw(maxLengthOutputKobs) << std::left << outputKobs[n]
				<< std::setw(maxLengthOutputKdif) << std::left << outputKdif[n]
				<< std::setw(maxLengthOutputKact) << std::left << outputKact[n]
				<< std::setw(maxLengthOutputReactionRadius) << std::left << outputReactionRadius[n]
				<< std::setw(maxLengthOutputReactionRadiusEff) << std::left << outputReactionRadiusEff[n]
				<< std::setw(maxLengthOutputReactionRadiusEff1) << std::left << outputReactionRadiusEff1[n]
				<< std::setw(maxLengthOutputProbability) << std::left << outputProbability[n] <<
				G4endl;
			} else if ( i == 5 ) {
				G4cout << std::setw(7) << " V "
				<< std::setw(maxLengthOutputReactionParentA) << std::left << outputReactionParentA[n] << " + "
				<< std::setw(maxLengthOutputReactionParentB) << std::left << outputReactionParentB[n] << " -> "
				<< std::setw(maxLengthOutputReactionProducts) << std::left << outputReactionProducts[n]
				<< std::setw(maxLengthOutputKobs) << std::left << outputKobs[n]
				<< std::setw(maxLengthOutputReactionRadius) << std::left << outputReactionRadius[n]
				<< std::setw(maxLengthOutputReactionRadiusEff) << std::left << outputReactionRadiusEff[n]
				<< std::setw(maxLengthOutputProbability) << std::left << outputProbability[n] <<
				G4endl;
			} else if ( i == 6 ) {
				G4String scavModel = outputModel[n] == "false" ? " Single factor " : " Double factor ";
				G4cout << std::setw(7) << " VI "
				<< std::setw(maxLengthOutputReactionParentA) << std::left << outputReactionParentA[n] << " + "
				<< std::setw(maxLengthOutputReactionParentB) << std::left << outputReactionParentB[n] << " -> "
				<< std::setw(maxLengthOutputReactionProducts) << std::left << outputReactionProducts[n]
				<< std::setw(maxLengthOutputKobs) << std::left << outputScav[n] << std::setw(15) << scavModel <<
				G4endl;
			}
			n++;
		}
		G4cout << G4endl;
	}
	
	delete[] title;
	delete[] outputProbability;
	delete[] outputReactionRadius;
	delete[] outputReactionRadiusEff;
	delete[] outputReactionRadiusEff1;
	delete[] outputKobs;
	delete[] outputKdif;
	delete[] outputKact;
	delete[] outputScav;
	delete[] outputAlpha;
	delete[] outputReactionProducts;
	delete[] outputReactionParentA;
	delete[] outputReactionParentB;
	delete[] outputModel;
	temporal.clear();
}


G4int TsIRTConfiguration::GetReactionIndex(G4int pdgA, G4int pdgB) {
	for ( size_t u = 0; u < fMoleculeCanReactWith[pdgA].size(); u++ ) {
		if ( pdgB == fMoleculeCanReactWith[pdgA][u].first )
				return fMoleculeCanReactWith[pdgA][u].second;
	}
	return -1;
}


G4double TsIRTConfiguration::GetOnsagerRadius(G4int molA, G4int molB) {
	if ( !fScaleForTemperature )
		return (fMoleculesDefinition[molA].charge * fMoleculesDefinition[molB].charge)/
		(4*CLHEP::pi*epsilon0*k_Boltzmann) / (293.15 * 80.1);
	else {
		G4double temperatureInKelvin = fTemperature + 273.15;
		G4double epsilon =((5321.0 * std::pow(temperatureInKelvin,-1)) + (233.76 * std::pow(temperatureInKelvin,+0))
						   - (0.9297 * std::pow(temperatureInKelvin,+1)) + (0.001417*std::pow(temperatureInKelvin,+2))
						   - (8.292E-7*std::pow(temperatureInKelvin,+3))) *  8.85418781762037E-12;
		
		G4double electronCharge = 1.60217662e-19 ;
		G4double rc = std::pow(electronCharge,2) / (4 * CLHEP::pi * epsilon * 1.3806488E-23 * temperatureInKelvin)
		* fMoleculesDefinition[molA].charge * fMoleculesDefinition[molB].charge;
		return rc*1E9*nm;
	}
}


void TsIRTConfiguration::ResolveRemainerReactionParameters() {
	for ( int i = 0; i < (int)fReactions.size(); i++ ) {
		G4int molA = fReactions[i].reactorA;
		G4int molB = fReactions[i].reactorB;
		G4int reactionType = fReactions[i].reactionType;
		
		G4double observedReactionRate = fReactions[i].kobs;
		G4double diffusionReactionRate = fReactions[i].kdif;
		G4double activationReactionRate = fReactions[i].kact;
		G4double probability = fReactions[i].probabilityOfReaction;
		G4double rc = fReactions[i].OnsagerRadius;
		
		G4double reactionRadius = 0;
		G4double effectiveReactionRadius = 0;
		G4double effectiveTildeReactionRadius = 0;
		G4double alpha = 0;
		
		if (reactionType == 1 || reactionType == 3 || reactionType == 5) {
			G4double sumDiffCoeff = 0;
			if (molA == molB) {
				sumDiffCoeff = fMoleculesDefinition[molA].diffusionCoefficient;
				effectiveReactionRadius = observedReactionRate / (4. * CLHEP::pi * sumDiffCoeff * CLHEP::Avogadro);
			} else {
				sumDiffCoeff = fMoleculesDefinition[molA].diffusionCoefficient
				+ fMoleculesDefinition[molB].diffusionCoefficient;
				effectiveReactionRadius = observedReactionRate / (4. * CLHEP::pi * sumDiffCoeff * CLHEP::Avogadro);
			}
			
			effectiveReactionRadius /= probability;
			
			if ( rc == 0 ) {
				reactionRadius = effectiveReactionRadius;
			} else {
				reactionRadius = rc/std::log(1 + rc/effectiveReactionRadius);
			}
			
		} else if ( reactionType == 2 || reactionType == 4 ) {
			// R = RA + RB, see Plante 2011 after Eq.17
			reactionRadius = fMoleculesDefinition[molA].radius + fMoleculesDefinition[molB].radius;
			G4double sumDiffCoeff = fMoleculesDefinition[molA].diffusionCoefficient + fMoleculesDefinition[molB].diffusionCoefficient;
			
			if ( reactionType == 2 ) {
				effectiveReactionRadius = reactionRadius;
				effectiveTildeReactionRadius = effectiveReactionRadius;
				alpha = (activationReactionRate + diffusionReactionRate)/(diffusionReactionRate*reactionRadius);
				
			} else {
				effectiveReactionRadius = -rc / (1-exp(rc / reactionRadius));
				rc /= nm;
				G4double kact = activationReactionRate/fPm->GetUnitValue("/M/s");
				G4double r = reactionRadius/nm;
				G4double reff = effectiveReactionRadius/nm;
				G4double v = kact * exp(rc/r) / (4 * CLHEP::pi * std::pow(r,2) * 6.022140857e-1);

				sumDiffCoeff /= nm*nm/s;
				effectiveTildeReactionRadius =  rc/(std::exp(rc/r) * (1 + sumDiffCoeff*rc/(std::pow(r,2)*v))-1) * nm;
				
				G4double nm3persToPerMPers = 6.022140857e-1;
				kact /= nm3persToPerMPers;
				kact *= nm*nm*nm/s;
				rc *= nm;
				sumDiffCoeff *= nm*nm/s;
				alpha = 0.25*kact/(CLHEP::pi*reactionRadius*reactionRadius) +
				rc*sumDiffCoeff/(reactionRadius*reactionRadius*(1.0 - std::exp(-rc/reactionRadius)));
			}
		} else {
			continue; // TypeVI
		}
		
		if ( fQualityAssurance ) {
			reactionRadius = 1.0 * nm;
			effectiveReactionRadius = 1.0 * nm;
			probability = 1;
			fReactions[i].reactionRadius = reactionRadius;
			fReactions[i].effectiveReactionRadius = effectiveReactionRadius;
			fReactions[i].probabilityOfReaction = probability;
		}
		
		fReactions[i].reactionRadius = reactionRadius;
		fReactions[i].effectiveReactionRadius = effectiveReactionRadius;
		fReactions[i].effectiveTildeReactionRadius = effectiveTildeReactionRadius;
		fReactions[i].alpha = alpha;
	}
}


void TsIRTConfiguration::CalculateContactProbabilities() {
	for ( int i = 0; i < (int)fReactions.size(); i++ ) {
		G4double probability = 1;
		G4int reactionType = fReactions[i].reactionType;
		G4int molA = fReactions[i].reactorA;
		G4int molB = fReactions[i].reactorB;
		if ( reactionType == 1 || reactionType == 3 || reactionType == 5 ) {
			if ( reactionType == 5 )
				probability = 0.25;
			else
				probability = 1.0;
			
		} else if ( reactionType == 2 || reactionType == 4 ){
			G4double Rs = 0.29 * nm;
			G4double diffusionReactionRate = fReactions[i].kdif;
			G4double observedReactionRate = fReactions[i].kobs;
			G4double rc = GetOnsagerRadius(molA, molB);
			
			G4double sigmaEff = fMoleculesDefinition[molA].radius + fMoleculesDefinition[molB].radius;
			G4double sigmaEffRs = fMoleculesDefinition[molA].radius + fMoleculesDefinition[molB].radius + Rs;

			if(reactionType == 4){
				sigmaEff = -rc / (1-exp(rc /sigmaEff));
				sigmaEffRs = -rc / (1-exp(rc /sigmaEffRs));
			}

			probability = observedReactionRate / diffusionReactionRate * ((1 - (sigmaEff/sigmaEffRs))
					/(1 - observedReactionRate / diffusionReactionRate * sigmaEff/sigmaEffRs));	

		} else {
			continue;
		}
		
		fReactions[i].probabilityOfReaction = probability;
	}
}


void TsIRTConfiguration::ResolveReactionRateCoefficients() {
	for ( int i = 0; i < (int)fReactions.size(); i++ ) {
		G4int molA = fReactions[i].reactorA;
		G4int molB = fReactions[i].reactorB;
		G4int reactionType = fReactions[i].reactionType;
		G4double kobs = fReactions[i].kobs;
		
		G4double rc = GetOnsagerRadius(molA, molB);
		G4double diffusionReactionRate = 0;
		G4double activationReactionRate = 0;
		fReactions[i].OnsagerRadius = rc;

		if (fReactions[i].kdif != -1 && fReactions[i].kact != -1) {
			continue;
		}
		
		if (reactionType == 1 || reactionType == 3 || reactionType == 5) {
			diffusionReactionRate = kobs;
			
		} else if ( reactionType == 2 || reactionType == 4 ) {
			
			// R = RA + RB, see Plante 2011 after Eq.17
			G4double ReactionRadius = fMoleculesDefinition[molA].radius + fMoleculesDefinition[molB].radius;
			G4double sumDiffCoeff = fMoleculesDefinition[molA].diffusionCoefficient + fMoleculesDefinition[molB].diffusionCoefficient;
			
			if ( reactionType == 2 ) {
				diffusionReactionRate = 4 * pi * sumDiffCoeff * ReactionRadius * Avogadro;
				if (molA == molB)
					diffusionReactionRate/=2;
				
				activationReactionRate = diffusionReactionRate * kobs / (diffusionReactionRate - kobs);
				
			} else {
				G4double effectiveReactionRadius = -rc / (1-exp(rc / ReactionRadius));
				diffusionReactionRate = 4 * pi * sumDiffCoeff * effectiveReactionRadius * Avogadro;
				
				if (molA == molB) diffusionReactionRate/=2;
				
				activationReactionRate = diffusionReactionRate * kobs / (diffusionReactionRate - kobs);
			}
		} else {
			continue; //Type VI;
		}
		
		fReactions[i].kdif = diffusionReactionRate;
		fReactions[i].kact = activationReactionRate;
	}
}


void TsIRTConfiguration::InsertReaction(G4int molA, G4int molB, std::vector<G4int> products,
										G4double kobs, G4int reactionType)
{
	G4int index = fReactionID;
	fReactionID++;
	fTotalBinaryReaction++;
	
	TsMolecularReaction aMolecularReaction;
	aMolecularReaction.reactorA = molA;
	aMolecularReaction.reactorB = molB;
	aMolecularReaction.products = products;
	aMolecularReaction.index = index;
	aMolecularReaction.kact  = -1;
	aMolecularReaction.kdif  = -1;
	
	if (fAllTotallyDiffusionControlled) {
		if ( reactionType == 2 ) {
			reactionType = 1;
		} else if ( reactionType == 4 ) {
			reactionType = 3;
		} else if ( reactionType == 5 ) {
			if ( fMoleculesDefinition[molA].charge * fMoleculesDefinition[molB].charge != 0 )
				reactionType = 3;
			else
				reactionType = 1;
		}
	}
	
	aMolecularReaction.reactionType = reactionType;
	aMolecularReaction.kobs = kobs;
	aMolecularReaction.scavengingCapacity = 0.0;
	
	G4int pdgA = molA;
	G4int pdgB = molB;
	if ( fMoleculeCanReactWith.find(pdgA) == fMoleculeCanReactWith.end() ) { // key not found
		fMoleculeCanReactWith[pdgA].push_back(std::make_pair(pdgB, index));
	} else {
		G4bool found = false;
		for ( size_t u = 0; u < fMoleculeCanReactWith[pdgA].size(); u++ ) {
			if ( pdgB == fMoleculeCanReactWith[pdgA][u].first) {
				found = true;
				break;
			}
		}
		if ( !found ) {
			fMoleculeCanReactWith[pdgA].push_back(std::make_pair(pdgB,index));
		}
	}
	
	if ( fMoleculeCanReactWith.find(pdgB) == fMoleculeCanReactWith.end() ) { // key not found
		fMoleculeCanReactWith[pdgB].push_back(std::make_pair(pdgA,index));
	} else {
		G4bool found = false;
		for ( size_t u = 0; u < fMoleculeCanReactWith[pdgB].size(); u++ ) {
			if ( pdgA == fMoleculeCanReactWith[pdgB][u].first ) {
				found = true;
				break;
			}
		}
		if ( !found ) {
			fMoleculeCanReactWith[pdgB].push_back(std::make_pair(pdgA,index));
		}
	}
	
	fReactions[index] = aMolecularReaction;
}


void TsIRTConfiguration::InsertReaction(G4String A, G4String B, G4String p1, G4String p2, G4String p3,
										G4double kobs, G4int reactionType)
{
	QuitIfMoleculeNotFound(A);
	QuitIfMoleculeNotFound(B);
	G4String molNameA  = fExistingMolecules[A];
	G4String molNameB  = fExistingMolecules[B];

	std::vector<G4int> products;
	G4String molNameP1 = p1;
	if (p1 != "None" && p1 != "none") {
		QuitIfMoleculeNotFound(p1);
		molNameP1 = fExistingMolecules[p1];
		products.push_back(fMoleculesID[molNameP1]);
	}

	G4String molNameP2 = p2;
	if (p2 != "None" && p2 != "none") {
		QuitIfMoleculeNotFound(p2);
		molNameP2 = fExistingMolecules[p2];
		products.push_back(fMoleculesID[molNameP2]);
	}

	G4String molNameP3 = p3;
	if (p3 != "None" && p3 != "none") {
		QuitIfMoleculeNotFound(p3);
		molNameP3 = fExistingMolecules[p3];
		products.push_back(fMoleculesID[molNameP3]);
	}

	G4int molA = fMoleculesID[molNameA];
	G4int molB = fMoleculesID[molNameB];
	
	if ( GetReactionIndex(molA, molB) > -1 ) {
		G4cerr << "TOPAS is exiting due to a fatal error in IRT Chemistry setup!" << G4endl;
		G4cerr << "--- Reaction: " << molNameA << " + " << molNameB << " -> " 
		                           << molNameP1 << " + " << molNameP2 << " + " << molNameP3 << " already exists." << G4endl;
		fPm->AbortSession(1);
		return;
	}
	
	InsertReaction(molA, molB, products, kobs, reactionType);
}


void TsIRTConfiguration::InsertReaction(G4String A, G4String B, std::vector<G4String> prods,
										G4double kobs, G4int reactionType)
{
	QuitIfMoleculeNotFound(A);
	QuitIfMoleculeNotFound(B);
	G4String molNameA  = fExistingMolecules[A];
	G4String molNameB  = fExistingMolecules[B];
	std::vector<G4int> products;

	for (size_t i = 0; i < prods.size(); i++) {
		G4String molName = prods[i];
		if (molName != "None" && molName != "none") {
			QuitIfMoleculeNotFound(molName);
			molName = fExistingMolecules[molName];
			products.push_back(fMoleculesID[molName]);
		}
	}

	G4int molA = fMoleculesID[molNameA];
	G4int molB = fMoleculesID[molNameB];
	
	if ( GetReactionIndex(molA, molB) > -1 ) {
		G4cerr << "TOPAS is exiting due to a fatal error in IRT Chemistry setup!" << G4endl;
		G4cerr << "--- Reaction: " << molNameA << " + " << molNameB << " -> ";
		for (size_t i = 0; i < prods.size(); i++) {
			if (i != prods.size() -1) 
				G4cerr << prods[i] << " + ";
			else
				G4cerr << prods[i];
		}
		G4cerr << " already exists." << G4endl;
		fPm->AbortSession(1);
		return;
	}
	
	InsertReaction(molA, molB, products, kobs, reactionType);
}


void TsIRTConfiguration::InsertBackgroundReaction(G4String A, G4String B, G4String p1,
												  G4String p2, G4String p3,
												  G4double scavengingCapacity, G4bool sampleExponential)
{
	QuitIfMoleculeNotFound(A);
	QuitIfMoleculeNotFound(B);
	G4String molNameA = fExistingMolecules[A];
	G4String molNameB = fExistingMolecules[B];
	
	G4int index = fReactionID;
	fReactionID++;
	G4int molA = fMoleculesID[molNameA];
	G4int molB = fMoleculesID[molNameB];
	// No need to test to confirm if this reaction already exists, because it could exists as a second order reaction

	std::vector<G4int> products;
	G4String molNameP1 = p1;
	if (p1 != "None" && p1 != "none") {
		QuitIfMoleculeNotFound(p1);
		molNameP1 = fExistingMolecules[p1];
		products.push_back(fMoleculesID[molNameP1]);
	}

	G4String molNameP2 = p2;
	if (p2 != "None" && p2 != "none") {
		QuitIfMoleculeNotFound(p2);
		molNameP2 = fExistingMolecules[p2];
		products.push_back(fMoleculesID[molNameP2]);
	}

	G4String molNameP3 = p3;
	if (p3 != "None" && p3 != "none") {
		QuitIfMoleculeNotFound(p3);
		molNameP3 = fExistingMolecules[p3];
		products.push_back(fMoleculesID[molNameP3]);
	}
	
	TsMolecularReaction aMolecularReaction;
	aMolecularReaction.reactorA = molA;
	aMolecularReaction.reactorB = molB;
	aMolecularReaction.products = products;
	aMolecularReaction.scavengingCapacity = scavengingCapacity;
	aMolecularReaction.kobs = 0.0;
	aMolecularReaction.index = index;
	aMolecularReaction.reactionType = 6;
	aMolecularReaction.sampleExponential = sampleExponential;

	fReactions[index] = aMolecularReaction;
}


void TsIRTConfiguration::InsertBackgroundReaction(G4String A, G4String B, G4String p1,
												  G4String p2, G4String p3, G4double kobs,
												  G4double concentration, G4bool sampleExponential)
{
	QuitIfMoleculeNotFound(A);
	QuitIfMoleculeNotFound(B);
	G4String molNameA = fExistingMolecules[A];
	G4String molNameB = fExistingMolecules[B];
	
	G4int index = fReactionID;
	fReactionID++;
	G4int molA = fMoleculesID[molNameA];
	G4int molB = fMoleculesID[molNameB];
	// No test to confirm if this reaction already exists, it could be a first order reaction
	std::vector<G4int> products;
	G4String molNameP1 = p1;
	if (p1 != "None" && p1 != "none") {
		QuitIfMoleculeNotFound(p1);
		molNameP1 = fExistingMolecules[p1];
		products.push_back(fMoleculesID[molNameP1]);
	}

	G4String molNameP2 = p2;
	if (p2 != "None" && p2 != "none") {
		QuitIfMoleculeNotFound(p2);
		molNameP2 = fExistingMolecules[p2];
		products.push_back(fMoleculesID[molNameP2]);
	}

	G4String molNameP3 = p3;
	if (p3 != "None" && p3 != "none") {
		QuitIfMoleculeNotFound(p3);
		molNameP3 = fExistingMolecules[p3];
		products.push_back(fMoleculesID[molNameP3]);
	}
	
	TsMolecularReaction aMolecularReaction;
	aMolecularReaction.reactorA = molA;
	aMolecularReaction.reactorB = molB;
	aMolecularReaction.products = products;
	aMolecularReaction.kobs = kobs;
	G4double sumDiffCoeff = fMoleculesDefinition[molA].diffusionCoefficient
	+ fMoleculesDefinition[molB].diffusionCoefficient;
	aMolecularReaction.effectiveReactionRadius = kobs / (4. * CLHEP::pi * sumDiffCoeff * CLHEP::Avogadro);
	aMolecularReaction.concentration = concentration;
	aMolecularReaction.scavengingCapacity = kobs * concentration;
	aMolecularReaction.index = index;
	aMolecularReaction.reactionType = 6;
	aMolecularReaction.sampleExponential = sampleExponential;// true;
	fReactions[index] = aMolecularReaction;
}


void TsIRTConfiguration::InsertBackgroundReaction(G4String A, G4String B, std::vector<G4String> p,
												  G4double scavengingCapacity, G4bool sampleExponential)
{
	QuitIfMoleculeNotFound(A);
	QuitIfMoleculeNotFound(B);
	G4String molNameA = fExistingMolecules[A];
	G4String molNameB = fExistingMolecules[B];
	
	G4int index = fReactionID;
	fReactionID++;
	G4int molA = fMoleculesID[molNameA];
	G4int molB = fMoleculesID[molNameB];
	// No need to test to confirm if this reaction already exists, because it could exists as a second order reaction

	std::vector<G4int> products;
	for (size_t i = 0; i < p.size(); i++) {
		G4String molNameP = p[i];
		if (p[i] != "None" && p[i] != "none") {
			QuitIfMoleculeNotFound(p[i]);
			molNameP = fExistingMolecules[p[i]];
			products.push_back(fMoleculesID[molNameP]);
		}
	}
	
	TsMolecularReaction aMolecularReaction;
	aMolecularReaction.reactorA = molA;
	aMolecularReaction.reactorB = molB;
	aMolecularReaction.products = products;
	aMolecularReaction.scavengingCapacity = scavengingCapacity;
	aMolecularReaction.kobs = 0.0;
	aMolecularReaction.index = index;
	aMolecularReaction.reactionType = 6;
	aMolecularReaction.sampleExponential = sampleExponential;
	
	fReactions[index] = aMolecularReaction;
}


void TsIRTConfiguration::InsertBackgroundReaction(G4String A, G4String B, std::vector<G4String> p,
												  G4double kobs, G4double concentration, G4bool sampleExponential)
{
	QuitIfMoleculeNotFound(A);
	QuitIfMoleculeNotFound(B);
	G4String molNameA = fExistingMolecules[A];
	G4String molNameB = fExistingMolecules[B];
	
	G4int index = fReactionID;
	fReactionID++;
	G4int molA = fMoleculesID[molNameA];
	G4int molB = fMoleculesID[molNameB];
	// No test to confirm if this reaction already exists, it could be a first order reaction
	
	std::vector<G4int> products;
	for (size_t i = 0; i < p.size(); i++) {
		G4String molNameP = p[i];
		if (p[i] != "None" && p[i] != "none") {
			QuitIfMoleculeNotFound(p[i]);
			molNameP = fExistingMolecules[p[i]];
			products.push_back(fMoleculesID[molNameP]);
		}
	}
	
	TsMolecularReaction aMolecularReaction;
	aMolecularReaction.reactorA = molA;
	aMolecularReaction.reactorB = molB;
	aMolecularReaction.products = products;
	aMolecularReaction.kobs = kobs;
	G4double sumDiffCoeff = fMoleculesDefinition[molA].diffusionCoefficient
	+ fMoleculesDefinition[molB].diffusionCoefficient;
	aMolecularReaction.effectiveReactionRadius = kobs / (4. * CLHEP::pi * sumDiffCoeff * CLHEP::Avogadro);
	aMolecularReaction.concentration = concentration;
	aMolecularReaction.scavengingCapacity = kobs * concentration;
	aMolecularReaction.index = index;
	aMolecularReaction.reactionType = 6;
	aMolecularReaction.sampleExponential = sampleExponential;// true;
	fReactions[index] = aMolecularReaction;
}


void TsIRTConfiguration::AdjustReactionAndDiffusionRateForTemperature() {
	std::map<G4String, std::pair<G4String,G4String>> ReactionLabels;
	std::map<G4String, G4double> DefaultValues;

	ReactionLabels["R1"]  = std::make_pair("e_aq^-1","e_aq^-1");
	DefaultValues["R1"]   = 5.5e9;
	ReactionLabels["R2"]  = std::make_pair("e_aq^-1","H3O^1");
    DefaultValues["R2"]   = 2.3e10;
	ReactionLabels["R3"]  = std::make_pair("e_aq^-1","H^0");
    DefaultValues["R3"]   = 2.5e10;
	ReactionLabels["R4"]  = std::make_pair("e_aq^-1","OH^0");
    DefaultValues["R4"]   = 3.0e10;
	ReactionLabels["R5"]  = std::make_pair("e_aq^-1","H2O2^0");
    DefaultValues["R5"]   = 1.1e10;
	ReactionLabels["R6"]  = std::make_pair("H3O^1",  "OH^-1");
    DefaultValues["R6"]   = 14.3e10;
	ReactionLabels["R7"]  = std::make_pair("H^0",    "H^0");
	DefaultValues["R7"]   = 7.8e9;
	ReactionLabels["R8"]  = std::make_pair("OH^0",   "H^0");
	DefaultValues["R8"]   = 1.55e10;
	ReactionLabels["R9"]  = std::make_pair("H^0",    "H2O2^0");
	DefaultValues["R9"]   = 9.0e7;
	ReactionLabels["R10"] = std::make_pair("OH^0",   "OH^0");
	DefaultValues["R10"]  = 5.5e9;

	for(size_t i = 0; i < fReactions.size(); i++) {
		G4double radiusA = fMoleculesDefinition[fReactions[i].reactorA].radius / m;
		G4double radiusB = fMoleculesDefinition[fReactions[i].reactorB].radius / m;
		G4String ReactA = fMoleculesName[fReactions[i].reactorA];
		G4String ReactB = fMoleculesName[fReactions[i].reactorB];

		G4double newKobs = -1;
		G4double newKdif = -1;
		G4double newKact = -1;

		G4double refKobs = -1;
		G4double refKdif = -1;
		G4double refKact = -1;

		G4double factor  = -1;

		G4int type = fReactions[i].reactionType;

		G4bool Found = false;
		if ((ReactA == ReactionLabels["R1"].first  && ReactB == ReactionLabels["R1"].second) ||
			(ReactA == ReactionLabels["R1"].second && ReactB == ReactionLabels["R1"].first)) {
			newKobs = 2*fUtils->ArrheniusFunction(2.33E13, fTemperature+273.15, 20.3);
			refKobs = 2*fUtils->ArrheniusFunction(2.33E13, 298.15, 20.3);
			factor  = 2*DefaultValues["R1"] / (refKobs);
			Found = true;
		}

		else if ((ReactA == ReactionLabels["R2"].first  && ReactB == ReactionLabels["R2"].second) ||
			     (ReactA == ReactionLabels["R2"].second && ReactB == ReactionLabels["R2"].first)) {
			G4double R = radiusA + radiusB;
			G4double Da = fUtils->ElectronDiffusionRate(fTemperature);
			G4double Db = fUtils->H3ODiffusionRate(fTemperature+273.15);
			newKobs = fUtils->ArrheniusFunction(1.24E12,fTemperature+273.15,10.1);
			newKdif = fUtils->DebyeFunction(1, R, Da+Db, fTemperature, -1, 1);
			newKact = fUtils->NoyesRelationship(newKobs, 0, newKdif);
			refKobs = fUtils->ArrheniusFunction(1.24E12,298.15,10.1);
			factor  = DefaultValues["R2"] / refKobs;
			Found = true;
		}

		else if ((ReactA == ReactionLabels["R3"].first  && ReactB == ReactionLabels["R3"].second) ||
			     (ReactA == ReactionLabels["R3"].second && ReactB == ReactionLabels["R3"].first)) {
			newKobs = fUtils->ArrheniusFunction(7.52E12, fTemperature+273.15, 14.0);
			refKobs = fUtils->ArrheniusFunction(7.52E12, 298.15, 14.0);
			factor  = DefaultValues["R3"] / refKobs;
			Found = true;
		}

		else if ((ReactA == ReactionLabels["R4"].first  && ReactB == ReactionLabels["R4"].second) ||
			     (ReactA == ReactionLabels["R4"].second && ReactB == ReactionLabels["R4"].first)) {
			G4double R = radiusA + radiusB;
			G4double D = fUtils->ElectronDiffusionRate(25.) + 2.3E-9;
			G4double Da = fUtils->ElectronDiffusionRate(fTemperature);
			G4double Db = fUtils->WaterDiffusionRate(fTemperature+273.15,2.3E-9);
			newKact = fUtils->ArrheniusFunction(3.04E10, fTemperature+273.15, -3.5);
			newKdif = fUtils->SmoluchowskiFunction(1, R, Da+Db);
			newKobs = fUtils->NoyesRelationship(0, newKact, newKdif);
			refKact = fUtils->ArrheniusFunction(3.04E10, 298.15, -3.5);
			refKdif = fUtils->SmoluchowskiFunction(1, R, D);
			refKobs = fUtils->NoyesRelationship(0,refKact,refKdif);
			factor  = DefaultValues["R4"] / refKobs;
			Found = true;
		}

		else if ((ReactA == ReactionLabels["R5"].first  && ReactB == ReactionLabels["R5"].second) ||
			     (ReactA == ReactionLabels["R5"].second && ReactB == ReactionLabels["R5"].first)) {
			G4double R = radiusA + radiusB;
			G4double D = fUtils->ElectronDiffusionRate(25.) + 2.3E-9;
			G4double Da = fUtils->ElectronDiffusionRate(fTemperature);
			G4double Db = fUtils->WaterDiffusionRate(fTemperature+273.15,2.3E-9);
			newKact = fUtils->ArrheniusFunction(6.26E12, fTemperature+273.15, 14.0);
			newKdif = fUtils->SmoluchowskiFunction(1, R, Da+Db);
			newKobs = fUtils->NoyesRelationship(0, newKact, newKdif);
			refKact = fUtils->ArrheniusFunction(6.26E12, 298.15, 14.0);
			refKdif = fUtils->SmoluchowskiFunction(1, R, D);
			refKobs = fUtils->NoyesRelationship(0,refKact,refKdif);
			factor  = DefaultValues["R5"] / refKobs;
			Found = true;
		}

		else if ((ReactA == ReactionLabels["R6"].first  && ReactB == ReactionLabels["R6"].second) ||
			     (ReactA == ReactionLabels["R6"].second && ReactB == ReactionLabels["R6"].first)) {
			G4double R = radiusA + radiusB;
			G4double D = fUtils->H3ODiffusionRate(298.15) + fUtils->OHmDiffusionRate(298.15);//diffusionA + diffusionB;
			G4double Da = fUtils->H3ODiffusionRate(fTemperature+273.15);
			G4double Db = fUtils->OHmDiffusionRate(fTemperature+273.15);
			newKobs = fUtils->DebyeFunction(1, R, Da+Db, fTemperature, 1, -1);
			refKobs = fUtils->DebyeFunction(1, R, D, 25.,-1,1);
			factor  = DefaultValues["R6"] / refKobs;
			Found = true;
		}

		else if ((ReactA == ReactionLabels["R7"].first  && ReactB == ReactionLabels["R7"].second) ||
			     (ReactA == ReactionLabels["R7"].second && ReactB == ReactionLabels["R7"].first)) {
			G4double R = radiusA + radiusB;
			G4double Da = fUtils->WaterDiffusionRate(fTemperature+273.15,8.0E-9);
			newKobs = 2*fUtils->SmoluchowskiFunction(1, R, 2*Da);
			refKobs = 2*fUtils->SmoluchowskiFunction(1, R, 1.6E-8);
			factor  = 2*DefaultValues["R7"] / (refKobs);
			Found = true;
		}

		else if ((ReactA == ReactionLabels["R8"].first  && ReactB == ReactionLabels["R8"].second) ||
			     (ReactA == ReactionLabels["R8"].second && ReactB == ReactionLabels["R8"].first)) {
			G4double R = radiusA + radiusB;
			G4double Da = fUtils->WaterDiffusionRate(fTemperature+273.15,8.0E-9);
			G4double Db = fUtils->WaterDiffusionRate(fTemperature+273.15,2.3E-9);
			newKact = fUtils->ArrheniusFunction(0.178E12, fTemperature+273.15, 4.5);
			newKdif = fUtils->SmoluchowskiFunction(1, R, Da+Db);
			newKobs = fUtils->NoyesRelationship(0,newKact,newKdif);
			refKact = fUtils->ArrheniusFunction(0.178E12, 298.15, 4.5);
			refKdif = fUtils->SmoluchowskiFunction(1, R, 10.3E-9);
			refKobs = fUtils->NoyesRelationship(0,refKact,refKdif);
			factor  = DefaultValues["R8"] / refKobs;
			Found = true;
		}

		else if ((ReactA == ReactionLabels["R9"].first  && ReactB == ReactionLabels["R9"].second) ||
			     (ReactA == ReactionLabels["R9"].second && ReactB == ReactionLabels["R9"].first)) {
			G4double R = radiusA + radiusB;
			G4double Da = fUtils->WaterDiffusionRate(fTemperature+273.15,8.0E-9);
			G4double Db = fUtils->WaterDiffusionRate(fTemperature+273.15,2.3E-9);
			newKobs = fUtils->ArrheniusFunction(3.21E10, fTemperature+273.15, 15.94);
			newKdif = fUtils->SmoluchowskiFunction(1, R, Da+Db);
			newKact = fUtils->NoyesRelationship(newKobs, 0, newKdif);
			refKobs = fUtils->ArrheniusFunction(3.21E10, 298.15, 15.94);
			factor  = DefaultValues["R9"] / refKobs;
			Found = true;
		}

		else if ((ReactA == ReactionLabels["R10"].first  && ReactB == ReactionLabels["R10"].second) ||
			     (ReactA == ReactionLabels["R10"].second && ReactB == ReactionLabels["R10"].first)) {
			G4double R = radiusA + radiusB;
			G4double Da = fUtils->WaterDiffusionRate(fTemperature+273.15,2.3E-9);
			newKact = fUtils->ArrheniusFunction(3.69E10, fTemperature+273.15, 3);
			newKdif = fUtils->SmoluchowskiFunction(1, R, 2*Da)/2;
			newKobs = fUtils->NoyesRelationship(0,newKact,newKdif);
			refKact = fUtils->ArrheniusFunction(3.69E10, 298.15, 3);
			refKdif = fUtils->SmoluchowskiFunction(1, R, 4.6E-9)/2;
			refKobs = fUtils->NoyesRelationship(0,refKact,refKdif);
			newKdif = -1;
			newKact = -1;
			factor  = DefaultValues["R10"] / refKobs;
			Found = true;
		}

		if (fKick == true)
			newKobs = newKobs * factor;

		if (Found) {
			if (type != 6) {
				if (newKobs > 0)
					fReactions[i].kobs = newKobs * fPm->GetUnitValue("/M/s");
				if (newKact > 0)
					fReactions[i].kact = newKact * fPm->GetUnitValue("/M/s");
				if (newKdif > 0)
					fReactions[i].kdif = newKdif * fPm->GetUnitValue("/M/s");
			}
		}
		else {
			fReactions[i].kobs = newKobs * fPm->GetUnitValue("/M/s");
			fReactions[i].scavengingCapacity = newKobs * (fReactions[i].scavengingCapacity / refKobs);
		}
	}

	for (auto& IndexAndMolDef:fMoleculesDefinition) {
		G4int Index                 = IndexAndMolDef.first;
		G4String MolName            = fMoleculesName[Index];

		if (MolName == "e_aq^-1") {
			G4double Diff = fUtils->ElectronDiffusionRate(fTemperature);
			fMoleculesDefinition[Index].diffusionCoefficient = Diff * (m*m/s);
		} else if (MolName == "H3O^1") {
			G4double Diff = fUtils->H3ODiffusionRate(fTemperature+273.15);
			fMoleculesDefinition[Index].diffusionCoefficient = Diff * (m*m/s);
		} else if (MolName == "OH^-1") {
			G4double Diff = fUtils->OHmDiffusionRate(fTemperature+273.15);
			fMoleculesDefinition[Index].diffusionCoefficient = Diff * (m*m/s);
		} else if (MolName == "OH^0") {
			G4double Diff = fUtils->WaterDiffusionRate(fTemperature+273.15,2.3E-9);
			fMoleculesDefinition[Index].diffusionCoefficient = Diff * (m*m/s);
		} else if (MolName == "H^0") {
			G4double Diff = fUtils->WaterDiffusionRate(fTemperature+273.15,8.0E-9);
			fMoleculesDefinition[Index].diffusionCoefficient = Diff * (m*m/s);
		} else if (MolName == "H_2^0") {
			G4double Diff = fUtils->WaterDiffusionRate(fTemperature+273.15,4.8E-9);
			fMoleculesDefinition[Index].diffusionCoefficient = Diff * (m*m/s);
		} else if (MolName == "H2O2^0") {
			G4double Diff = fUtils->WaterDiffusionRate(fTemperature+273.15,2.3E-9);
			fMoleculesDefinition[Index].diffusionCoefficient = Diff * (m*m/s);
		}
	}
}

void TsIRTConfiguration::AdjustReactionRateForPH(G4String pHOrConcentration) {
	std::vector<G4double> AcidComponents;
	G4double Ionic   = 0.0;
	G4double HCon    = 0.0;
	G4double HCon25  = 1.00E-7;
	G4double OHCon   = 0.0;
	G4double OHCon25 = 1.00E-7;
	G4double HSO4Con = 0.0;
	
	if (pHOrConcentration == "PH" && fpHSolvent == "h2so4") {
		AcidComponents = GetH2SO4ComponentsConcentrationPH(fpHValue);
		Ionic          = GetIonicStrength(AcidComponents);
		HCon           = AcidComponents[0];
		HSO4Con        = AcidComponents[1];
		OHCon          = AcidComponents[3];
		G4cout << "-- Adjust for PH of H2SO4 " << G4endl;
		
	}
	
	else if (pHOrConcentration == "Concentration" && fpHSolvent == "h2so4") {
		AcidComponents = GetH2SO4ComponentsConcentrationP(fpHSolventConcentration/fPm->GetUnitValue("M"));
		Ionic          = GetIonicStrength(AcidComponents);
		HCon           = AcidComponents[0];
		HSO4Con        = AcidComponents[1];
		OHCon          = AcidComponents[3];
		G4cout << "-- Adjust for Concentration of H2SO4 " << G4endl;
	}
	
	else if (fpHSolvent == "generic") {
		HCon  = pow(10,-fpHValue);
		OHCon = 1E-14 / HCon;
		AcidComponents = {HCon, 0.0, 0.0, OHCon, 0.0, 0.0};
		Ionic          = GetIonicStrength(AcidComponents);
		G4cout << "-- Adjust for a generic substance " << G4endl;
		
	} else {
		G4cout << "-- Is doing nothing " << pHOrConcentration << " " << fpHSolvent << G4endl;
	}
	/*if ( fpHSolvent == "h2so4" ) {
	 if ( pHOrConcentration == "PH" ) {
	 AcidComponents = GetH2SO4ComponentsConcentrationPH(fpHValue);
	 } else {
	 AcidComponents = GetH2SO4ComponentsConcentrationP(fpHSolventConcentration);
	 }
	 
	 Ionic          = GetIonicStrength(AcidComponents);
	 HCon           = AcidComponents[0];
	 HSO4Con        = AcidComponents[1];
	 OHCon          = AcidComponents[3];
	 } else {
	 HCon  = pow(10,-fpHValue);
	 OHCon = 1E-14 / HCon;
	 AcidComponents = {HCon, 0.0, 0.0, OHCon, 0.0, 0.0};
	 Ionic          = GetIonicStrength(AcidComponents);
	 }*/
	
	G4cout << G4endl;
	G4cout << " ###-------- pH Scaling Starts ---------###" << G4endl;
	
	for(size_t i = 0; i < fReactions.size(); i++) {
		G4int chargeA   = fMoleculesDefinition[fReactions[i].reactorA].charge;
		G4int chargeB   = fMoleculesDefinition[fReactions[i].reactorB].charge;
		G4String ReactA = fMoleculesName[fReactions[i].reactorA];
		G4String ReactB = fMoleculesName[fReactions[i].reactorB];
		std::vector<G4String> products;
		G4double k_Before = 0;
		
		for (size_t vsize = 0; vsize < fReactions[i].products.size(); vsize++){
			products.push_back(fMoleculesName[fReactions[i].products[vsize]]);
		}
		
		fReactions[i].kobs /= fPm->GetUnitValue("/M/s");
		fReactions[i].scavengingCapacity /= 1/s;
		
		if ((chargeA != 0 && chargeB != 0) && (!(ReactA == "e_aq^-1" && ReactB == "e_aq^-1"))) {
			if (fReactions[i].reactionType != 6) {
				k_Before = fReactions[i].kobs;
				fReactions[i].kobs = IonicRate(Ionic, fReactions[i]);
			}
			
			else if (fReactions[i].reactionType == 6 && ReactB == "H3O^1") {
				k_Before = fReactions[i].scavengingCapacity;
				fReactions[i].scavengingCapacity = IonicRate(Ionic, fReactions[i]) * HCon;
			}
			
			else if (fReactions[i].reactionType == 6 && ReactB == "OH^-1") {
				k_Before = fReactions[i].scavengingCapacity;
				fReactions[i].scavengingCapacity = IonicRate(Ionic, fReactions[i]) * OHCon;
			}
			
			else {
				G4cout << "========= "<< fReactions[i].reactionType << G4endl;
				k_Before = fReactions[i].kobs;
				fReactions[i].kobs = fReactions[i].kobs * 1E-7;
				fReactions[i].kobs = IonicRate(Ionic, fReactions[i]);
			}
		}
		
		else if ((fReactions[i].reactionType == 6) && (ReactB == "H3O^1")) {
			k_Before = fReactions[i].scavengingCapacity;
			fReactions[i].scavengingCapacity = fReactions[i].scavengingCapacity * HCon / HCon25;
		}
		
		else if ((fReactions[i].reactionType == 6) && (ReactB == "OH^-1")) {
			k_Before = fReactions[i].scavengingCapacity;
			fReactions[i].scavengingCapacity = fReactions[i].scavengingCapacity * OHCon / OHCon25;
		}
		
		else if ((fReactions[i].reactionType == 6) && (ReactB == "HSO4^-1")) {
			k_Before = fReactions[i].scavengingCapacity;
			fReactions[i].scavengingCapacity = fReactions[i].scavengingCapacity * HSO4Con;
		}
		
		if (k_Before > 0) {
			G4cout << " Reaction Type: " << fReactions[i].reactionType << " | Reaction: " << ReactA << " + " << ReactB;
			for (int prod = 0; prod < (int)products.size(); prod++) {
				if (prod == 0)
					G4cout << " -> ";
				G4cout << products[prod];
				if (prod < (int)products.size() - 1) {
					G4cout << " + ";
				}
			}
			if ( fReactions[i].reactionType == 6 )
				G4cout << G4endl << " ---- scav: " << k_Before << " ---> "  << fReactions[i].scavengingCapacity << G4endl << G4endl;
			else
				G4cout << G4endl << " ---- kobs: " << k_Before << " ---> "  << fReactions[i].kobs << G4endl << G4endl;
		}
		fReactions[i].kobs *= fPm->GetUnitValue("/M/s");
		fReactions[i].scavengingCapacity *= 1/s;
	}
	
	G4cout << " ###-------- pH Scaling Ends ---------###" << G4endl;
	G4cout << G4endl;
}


void TsIRTConfiguration::ResampleReactantsPosition(TsMolecule& molA, TsMolecule& molB, G4int index, G4double time) {
	// Position approach
	G4double D1 = fMoleculesDefinition[molA.id].diffusionCoefficient;
	G4double D2 = fMoleculesDefinition[molB.id].diffusionCoefficient;
	
	G4ThreeVector r1 = molA.position;
	G4ThreeVector r2 = molB.position;
	G4double dtA = std::fabs(time-molA.time);
	G4double dtB = std::fabs(time-molB.time);
	G4double dt  = dtA > dtB ? dtA : dtB;
	
	if ( D1 == 0 ) {
		molB.position = r1;
		return;
	} else if ( D2 == 0 ) {
		molA.position = r2;
		return;
	}
	if (dt == 0) {
		dt = 1*ps;
	}
	
	G4ThreeVector S1 = r1-r2;
	G4double r0 = S1.mag();
	G4double effectiveReactionRadius = fReactions[index].effectiveReactionRadius;
	if ( fReactions[index].reactionType == 4 )
		effectiveReactionRadius = fReactions[index].effectiveTildeReactionRadius;
	
	if (S1 != G4ThreeVector())
		S1.setMag(effectiveReactionRadius);

	G4double s12 = 2.0 * D1 * dt;
	G4double s22 = 2.0 * D2 * dt;
	G4double alpha = effectiveReactionRadius * r0/(2*(D1+D2)*dt);
	
	G4ThreeVector S2 = (r1 + (s12/s22)*r2) + G4ThreeVector(G4RandGauss::shoot(0.0, s12 + s22*s22/s12),
														   G4RandGauss::shoot(0.0, s12 + s22*s22/s12),
														   G4RandGauss::shoot(0.0, s12 + s22*s22/s12));
	
	if (S1 != G4ThreeVector()) {
		S1.setPhi(rad*G4UniformRand()*2.0*CLHEP::pi);
		S1.setTheta(rad*std::acos(1.0 + 1./alpha * std::log(1.0 - G4UniformRand()*(1.-std::exp(-2.0*alpha)))));
	}

	G4ThreeVector R1 = (D1*S1 + D2*S2)/(D1+D2);
	G4ThreeVector R2 = D2*(S2-S1)/(D1+D2);

	molA.position = R1;
	molA.time = time;
	molB.position = R2;
	molB.time = time;
	
}


std::vector<G4ThreeVector> TsIRTConfiguration::GetPositionOfProducts(TsMolecule molA, TsMolecule molB, G4int index) {
	G4double D1 = fMoleculesDefinition[molA.id].diffusionCoefficient;
	G4double D2 = fMoleculesDefinition[molB.id].diffusionCoefficient;
	G4ThreeVector r1 = molA.position;
	G4ThreeVector r2 = molB.position;
	
	std::vector<G4ThreeVector> result;
	// weighted-position
	G4ThreeVector position = r1*std::sqrt(D2)/(std::sqrt(D1) + std::sqrt(D2)) + r2*std::sqrt(D1)/(std::sqrt(D1) + std::sqrt(D2));
	if ( fReactions[index].products.size() == 1 ) {
		// At weighted position.
		//result.push_back(position);
		if ( G4UniformRand() > 0.5 )
			result.push_back(r1);
		else
			result.push_back(r2);
		
	} else if ( fReactions[index].products.size() == 3 ) {
		// at weighted and at parents position
		result.push_back(position);
		result.push_back(r1);
		result.push_back(r2);
	} else if ( fReactions[index].products.size() == 2 ) {
		// at parents position
		result.push_back(r1);
		result.push_back(r2);
	}
	else if ( fReactions[index].products.size() == 4) {
		result.push_back(position);
		result.push_back(r1);
		result.push_back(r2);
		result.push_back(r1 + 0.5*(r1+r2));
	}
	
	return result;
	
}

G4double TsIRTConfiguration::GetRCutOff(G4double tCutOff) {
	G4double probabilityOfReaction = 0.01;
	G4double maximumReactionRadius = 1.45*nm;
	G4double maximumRelativeDiffusionCoefficient = 2.0*9.46e9 *nm*nm/s;
	G4double erfcInv = fUtils->erfcInv(probabilityOfReaction);
	return maximumReactionRadius + 2.0 * std::sqrt(maximumRelativeDiffusionCoefficient * tCutOff) * erfcInv;
}


G4double TsIRTConfiguration::GetIndependentReactionTime(TsMolecule molA, TsMolecule molB, G4int indexOfReaction) {
	G4int typeOfReaction = fReactions[indexOfReaction].reactionType;
	G4double result = -1.0*ps;
	
	if ( typeOfReaction == 1 || typeOfReaction == 3 || typeOfReaction == 5) {
		result = SampleIRTTotallyDiffusionControlled(molA, molB, indexOfReaction);
	} else if ( typeOfReaction == 2 || typeOfReaction == 4 ) {
		result = SampleIRTPartiallyDiffusionControlled(molA, molB, indexOfReaction);
	}
	
	return result;
}


G4double TsIRTConfiguration::SampleIRTTotallyDiffusionControlled(TsMolecule molA, TsMolecule molB, G4int indexOfReaction) {
	G4int typeOfReaction = fReactions[indexOfReaction].reactionType;
	G4double probFactor = 1.0;
	G4double reff = (molA.position - molB.position).mag();
	
	G4double Reff = fReactions[indexOfReaction].effectiveReactionRadius;
	if ( typeOfReaction == 4 ) // Only if considering  all reactions totally diffusion controlled (for testing purposes)
		Reff = fReactions[indexOfReaction].effectiveTildeReactionRadius;
	
	if ( fReactions[indexOfReaction].OnsagerRadius != 0 ) {
		G4double rc = fReactions[indexOfReaction].OnsagerRadius;
		reff = -rc/(1 - std::exp(rc/reff));
	}
	
	G4double Winf = probFactor * Reff/reff;
	G4double D = fMoleculesDefinition[molA.id].diffusionCoefficient + fMoleculesDefinition[molB.id].diffusionCoefficient;
	
	G4double W = G4UniformRand();
	G4double irt = -1.0 * ps;
	
	if ( W < Winf ) {
		irt = (0.25/D) * std::pow( (reff-Reff)/fUtils->erfcInv(reff*W/Reff), 2 );
	}

	return irt;
}


void TsIRTConfiguration::TestSampling(G4int indexOfReaction, G4int nHistories) {
	G4double sigma = 1.0 * nm;
	G4int molA = -1, molB = -1;
	for ( auto& reactions : fReactions ) {
		if ( reactions.second.index == indexOfReaction ) {
			molA = reactions.second.reactorA;
			molB = reactions.second.reactorA;
			break;
		}
	}
	G4ThreeVector posA;
	G4ThreeVector posB;
	std::ofstream out("TestSampling.csv");
	out << "# Reaction ID " << indexOfReaction << std::endl;
	for ( int i = 0; i < nHistories; i++ ) {
		TsMolecule A;
		A.id = molA;
		A.time = 0.01*ps;
		A.reacted = false;
		A.trackID = 0;
		A.spin = -1;
		posA = G4ThreeVector(
							 G4RandGauss::shoot(0., 1.0)*sigma,
							 G4RandGauss::shoot(0., 1.0)*sigma,
							 G4RandGauss::shoot(0., 1.0)*sigma);
		A.position = posA;
		
		TsMolecule B;
		B.id = molB;
		B.time = 0.01*ps;
		B.reacted = false;
		B.trackID = 0;
		B.spin = -1;
		posB = G4ThreeVector(
							 G4RandGauss::shoot(0., 1.0)*sigma,
							 G4RandGauss::shoot(0., 1.0)*sigma,
							 G4RandGauss::shoot(0., 1.0)*sigma);
		B.position = posB;
		
		G4double irt = GetIndependentReactionTime(A, B, indexOfReaction);
		if ( irt > 0.0 )
			out << (irt)/ps << std::endl;
	}
	out.close();
}


G4double TsIRTConfiguration::SampleIRTPartiallyDiffusionControlled(TsMolecule molA, TsMolecule molB, G4int indexOfReaction) {

	G4double r0 = (molA.position-molB.position).mag();

	G4double sigma = 0;
    G4double kdif = fReactions[indexOfReaction].kdif;
    G4double kobs = fReactions[indexOfReaction].kobs;

	if ( fReactions[indexOfReaction].reactionType == 2 ) {
	    sigma = fReactions[indexOfReaction].reactionRadius;
	} else {
	    G4double rc = fReactions[indexOfReaction].OnsagerRadius;
		sigma = fReactions[indexOfReaction].effectiveReactionRadius;
		r0 = -rc/(1-std::exp(rc/r0));
	}

    if(sigma/r0 * kobs / kdif > G4UniformRand()) {
        G4double kact = fReactions[indexOfReaction].kact;
    	G4double D = fMoleculesDefinition[molA.id].diffusionCoefficient + fMoleculesDefinition[molB.id].diffusionCoefficient;
    	G4double a = kact / kobs / sigma;
    	G4double b = (r0 - sigma) /2;
    	return fUtils->SamplePDC(a, b)/D;
    }

	return -1;
}


G4double TsIRTConfiguration::SolveTime(TsMolecule, G4int indexOfReaction, G4double offset) {
	G4double Reff = (fReactions[indexOfReaction]).effectiveReactionRadius;
	G4double CsB = (fReactions[indexOfReaction]).probabilityOfReaction/s;
	G4double nm3persToPerMPers = 6.022140857e-1;
	G4double k = (fReactions[indexOfReaction]).kobs/nm3persToPerMPers * (nm*nm*nm/s);
	
	G4double sqrtArgument = std::pow(4.0*Reff*std::sqrt(Reff/k),2) - 4.0/CsB * std::log(1.0 - offset);
	if (sqrtArgument < 0 )
		return -1.0*ps;
	G4double t1 = std::pow(-2.0*Reff*sqrt(Reff/k) + 0.5 * sqrtArgument, 2);
	G4double t2 = std::pow(-2.0*Reff*sqrt(Reff/k) - 0.5 * sqrtArgument, 2);
	
	if ( t1 <= t2 )
		return t1;
	else
		return t2;
}


G4int  TsIRTConfiguration::ContactFirstOrderAndBackgroundReactions(TsMolecule molA) {
	G4int pdgA = molA.id;
	std::vector<size_t> index;
	for ( size_t u = 0; u < fReactions.size(); u++ ) {
		if (fReactions[u].reactorA == pdgA && fReactions[u].reactionType == 6 ) {
			index.push_back(u);
		}
	}
	
	size_t sizeIndex = index.size();
	if ( 0 < sizeIndex ) {
		std::random_device rd;
		std::mt19937 RandomGenerator(rd());
		std::shuffle(std::begin(index),std::end(index),RandomGenerator);
	}
	//Pimblott S M, 1991 Stochastic models of spur kinetics in water. International Journal of Radiation Applications and Instrumentation. Part 37 377–88
	for ( size_t v = 0; v < sizeIndex; v++ ) {
		size_t u = index[v];
		G4double R  = fMoleculesDefinition[pdgA].radius/nm;
		G4double R3 = R*R*R;
		G4double Cs = fReactions[u].concentration/fPm->GetUnitValue("M");
		if (Cs > 50) {continue;}          // No Contact Reactions with water molecules
		Cs *= 6.022140857e-1;  // nm3 to M multiply by 10^-24 Nav = 6.02214076×10^23x10^-24 = 0.602214076
		G4double prob1 = std::exp(-1.33333333*CLHEP::pi*R3*Cs);

		if ( G4UniformRand() < 1. - prob1 ) {
			return (int)u;
		}
	}
	
	return -1;
}


G4double TsIRTConfiguration::SampleExponentialTime(G4int pdgA, G4int pdgB, G4int indexOfReaction) {
	G4double D = fMoleculesDefinition[pdgA].diffusionCoefficient + fMoleculesDefinition[pdgB].diffusionCoefficient;
	D /= nm*nm/s;
	G4double Rreact = fReactions[indexOfReaction].effectiveReactionRadius;
	Rreact /= nm;
	G4double Cs = fReactions[indexOfReaction].concentration/0.60221/fPm->GetUnitValue("M");
	G4double A = -Rreact/std::sqrt(CLHEP::pi * D);
	G4double B = 0.5/std::sqrt(CLHEP::pi * D);
	G4double C = 4.0 * Rreact - std::log(1.0 - G4UniformRand())/(Rreact*Cs);
	G4double timeA = A + B * std::sqrt(C);
	G4double timeB = A - B * std::sqrt(C);
	if ( timeA < timeB )
		return timeA*timeA*s;
	else
		return timeB*timeB*s;
}


std::vector<std::pair<G4int, G4double>> TsIRTConfiguration::SampleAllIRTFirstOrderAndBackgroundReactions(TsMolecule molA ) {
	std::vector<std::pair<G4int, G4double>> result;
	G4int pdgA = molA.id;
	G4double scavengingCapacity, prob, time;
	
	for ( size_t u = 0; u < fReactions.size(); u++ ) { // TODO-> Review this algorithm: use t or dt???
		if ( fReactions[u].reactionType == 6 && pdgA == fReactions[u].reactorA) {
			scavengingCapacity = fReactions[u].scavengingCapacity;
			prob = G4UniformRand();
			if ( fReactions[u].sampleExponential )
				time = SampleExponentialTime(pdgA,fReactions[u].reactorB, G4int(u));
			else
				time = -(std::log(1.0 - prob)/scavengingCapacity);
			
			if ( time > 0.0 )
				result.push_back(std::make_pair((int)u, time));
		}
	}
	return result;
}


std::pair<G4int, G4double> TsIRTConfiguration::SampleIRTFirstOrderAndBackgroundReactions(TsMolecule molA ) {
	G4int pdgA = molA.id;
	// Search for the IRT with the minium value
	G4double irt = 1000.0*s;
	
	G4int index = -1;
	G4bool found = false;
	
	G4double scavengingCapacity, prob, time;
	
	for ( size_t u = 0; u < fReactions.size(); u++ ) { // TODO-> Review this algorithm: use t or dt???
		if ( fReactions[u].reactionType == 6 && pdgA == fReactions[u].reactorA) {
			scavengingCapacity = fReactions[u].scavengingCapacity;
			prob = G4UniformRand();
			if ( fReactions[u].sampleExponential )
				time = SampleExponentialTime(pdgA,fReactions[u].reactorB, G4int(u));
			else
				if (scavengingCapacity <= 0)
					time = -1;
				else
					time = -(std::log(1.0 - prob)/scavengingCapacity);
			
			if ( time < irt && time > 0.0) {
				irt = time;
				index = (int)u;
			}
			found = true;
		}
	}
	
	if ( !found ) {
		return std::make_pair(-1,-1*ps);
	}
	
	return std::make_pair(index, irt);
}


G4double TsIRTConfiguration::CalculateProbabilityPartiallyDiffusionControlled(TsMolecule molA, TsMolecule molB, G4int indexOfReaction, G4double t) {
	G4double DA = fMoleculesDefinition[molA.id].diffusionCoefficient;
	G4double DB = fMoleculesDefinition[molB.id].diffusionCoefficient;
	G4double D = DA + DB;
	G4double r0 = (molA.position-molB.position).mag();
	
	if ( fReactions[indexOfReaction].reactionType == 2 ) { // between neutral particles
		// From paper A Monte Carlo step-by-step
		G4double sigma = (fReactions[indexOfReaction]).effectiveReactionRadius;
		G4double kact = (fReactions[indexOfReaction]).kact;
		
		G4double alpha =(fReactions[indexOfReaction]).alpha;
		G4double factor = kact/(4.0*CLHEP::pi*sigma*D*r0*alpha);
		G4double x = (r0-sigma)/std::sqrt(4.0*D*t);
		G4double y = alpha*std::sqrt(D*t);
		
		return factor * (fUtils->erfc(x) - std::exp(-x*x) * fUtils->erfcx(x+y));
		
	} else {
		G4double rc = (fReactions[indexOfReaction]).OnsagerRadius;
		G4double reff = -rc/(1 - std::exp(rc/r0));
		G4double sigma = (fReactions[indexOfReaction]).reactionRadius;
		G4double sigmaEffEff = (fReactions[indexOfReaction]).effectiveTildeReactionRadius;
		
		G4double alpha = (fReactions[indexOfReaction]).alpha;
		G4double a = (4.0*sigma*sigma*alpha/(rc*rc))*std::sqrt(t/D)*std::sinh(0.5*rc/sigma)*std::sinh(0.5*rc/sigma);
		G4double b = (0.25*rc/std::sqrt(D*t))*(std::cosh(0.5*rc/r0)/std::sinh(0.5*rc/r0) - std::cosh(0.5*rc/sigma)/std::sinh(0.5*rc/sigma));
		
		G4double factor = sigmaEffEff/reff;
		
		return factor * (fUtils->erfc(b) - std::exp(-b*b) * fUtils->erfcx(a+b));
	}
}


void TsIRTConfiguration::Quit(const G4String& name, G4String message) {
	G4cerr << G4endl;
	G4cerr << "Topas is exiting due to a serious error in Chemistry IRT setup." << G4endl;
	G4cerr << "--- Parameter name: " << name << G4endl;
	G4cerr << "--- " << message << G4endl;
	G4cerr << G4endl;
	fPm->AbortSession(1);
}


void TsIRTConfiguration::SetTimeLimits(G4double lower, G4double upper) {
	fUpperTime = upper;
	fLowerTime = lower;
}


G4int TsIRTConfiguration::GetNumberOfReactions() {
	return (G4int)fReactions.size();
}


std::pair<G4String, G4String> TsIRTConfiguration::GetReactants(G4int ReactionIndex) {
	std::pair<G4String, G4String> Reactants;
	
	Reactants.first  = fMoleculesName[fReactions[ReactionIndex].reactorA];
	Reactants.second = fMoleculesName[fReactions[ReactionIndex].reactorB];
	
	return Reactants;
}


std::vector<G4String> TsIRTConfiguration::GetProducts(G4int ReactionIndex) {
	std::vector<G4String> Products;
	
	for ( size_t i = 0; i < (fReactions[ReactionIndex].products).size(); i++) {
		Products.push_back(fMoleculesName[(fReactions[ReactionIndex].products)[i]]);
	}
	
	return Products;
}


G4double TsIRTConfiguration::brents_fun(TsMolecule molA, TsMolecule molB, G4int indexOfReaction, G4double offset) {
	G4double lower = fLowerTime;
	G4double upper = fUpperTime;
	G4double tol = 0.001*ps;
	unsigned int max_iter = 1000;
	
	G4double a = lower;
	G4double b = upper;
	G4double fa = CalculateProbabilityPartiallyDiffusionControlled(molA, molB, indexOfReaction, a)+offset;
	G4double fb = CalculateProbabilityPartiallyDiffusionControlled(molA, molB, indexOfReaction, b)+offset;
	G4double fs = 0;
	
	if (!(fa * fb < 0)) {
		return -11;
	}
	
	if (std::abs(fa) < std::abs(b)) {
		std::swap(a,b);
		std::swap(fa,fb);
	}
	
	G4double c = a;
	G4double fc = fa;
	G4bool mflag = true;
	G4double ss = 0;
	G4double d = 0;
	
	for (unsigned int iter = 1; iter < max_iter; ++iter) {
		// stop if converged on root or error is less than tolerance
		if (std::abs(b-a) < tol) {
			return ss;
		}
		
		if (fa != fc && fb != fc) { // use inverse quadratic interopolation
			ss =      ( a * fb * fc / ((fa - fb) * (fa - fc)) )
			+ ( b * fa * fc / ((fb - fa) * (fb - fc)) )
			+ ( c * fa * fb / ((fc - fa) * (fc - fb)) );
		}
		else{ // secant method
			ss = b - fb * (b - a) / (fb - fa);
		}
		
		// checks to see whether we can use the faster converging quadratic && secant methods or if we need to use bisection
		if (    ( (ss < (3 * a + b) * 0.25) || (ss > b) ) ||
			( mflag && (std::abs(ss-b) >= (std::abs(b-c) * 0.5)) ) ||
			( !mflag && (std::abs(ss-b) >= (std::abs(c-d) * 0.5)) ) ||
			( mflag && (std::abs(b-c) < tol) ) ||
			( !mflag && (std::abs(c-d) < tol))    ) {
			// bisection method
			ss = (a+b)*0.5;
			mflag = true;
		}
		else {
			mflag = false;
		}
		
		fs = CalculateProbabilityPartiallyDiffusionControlled(molA, molB, indexOfReaction, ss)+offset;
		d = c;
		c = b;
		fc = fb;
		
		if ( fa * fs < 0) {
			b = ss;
			fb = fs;
		} else {
			a = ss;
			fa = fs;
		}
		
		if (std::abs(fa) < std::abs(fb)) {
			std::swap(a,b);
			std::swap(fa,fb);
		}
	}
	return -1.0*ps;
}


G4double TsIRTConfiguration::CalculateProbabilityOfScavenger(TsMolecule, G4int indexOfReaction, G4double t) {
	// Vars: Reff, k, [B]=Cs
	G4double Reff = (fReactions[indexOfReaction]).effectiveReactionRadius;
	G4double CsB = (fReactions[indexOfReaction]).probabilityOfReaction/s;
	G4double nm3persToPerMPers = 6.022140857e-1;
	G4double k = (fReactions[indexOfReaction]).kobs/nm3persToPerMPers * (nm*nm*nm/s);;
	
	G4double Wscav = 1.0 - std::exp(-CsB * (t + 4.0 * Reff * std::sqrt(t * Reff / k )));
	return Wscav;
}


G4double TsIRTConfiguration::brents_fun_scav(TsMolecule molA, G4int indexOfReaction, G4double offset) {
	G4double lower = fLowerTime;
	G4double upper = fUpperTime;
	G4double tol = 0.001*ps;
	unsigned int max_iter = 100000;
	
	G4double a = lower;
	G4double b = upper;
	G4double fa = CalculateProbabilityOfScavenger(molA, indexOfReaction, a)+offset;
	G4double fb = CalculateProbabilityOfScavenger(molA, indexOfReaction, b)+offset;
	G4double fs = 0;
	
	if (!(fa * fb < 0)) {
		return -11;
	}
	
	if (std::abs(fa) < std::abs(b)) {
		std::swap(a,b);
		std::swap(fa,fb);
	}
	
	G4double c = a;
	G4double fc = fa;
	G4bool mflag = true;
	G4double ss = 0;
	G4double d = 0;
	
	for (unsigned int iter = 1; iter < max_iter; ++iter) {
		// stop if converged on root or error is less than tolerance
		if (std::abs(b-a) < tol) {
			return ss;
		}
		
		if (fa != fc && fb != fc) { // use inverse quadratic interopolation
			ss =      ( a * fb * fc / ((fa - fb) * (fa - fc)) )
			+ ( b * fa * fc / ((fb - fa) * (fb - fc)) )
			+ ( c * fa * fb / ((fc - fa) * (fc - fb)) );
		}
		else{ // secant method
			ss = b - fb * (b - a) / (fb - fa);
		}
		
		// checks to see whether we can use the faster converging quadratic && secant methods or if we need to use bisection
		if (    ( (ss < (3 * a + b) * 0.25) || (ss > b) ) ||
			( mflag && (std::abs(ss-b) >= (std::abs(b-c) * 0.5)) ) ||
			( !mflag && (std::abs(ss-b) >= (std::abs(c-d) * 0.5)) ) ||
			( mflag && (std::abs(b-c) < tol) ) ||
			( !mflag && (std::abs(c-d) < tol))    ) {
			// bisection method
			ss = (a+b)*0.5;
			mflag = true;
		}
		else {
			mflag = false;
		}
		
		fs = CalculateProbabilityOfScavenger(molA, indexOfReaction, ss)+offset;
		d = c;
		c = b;
		fc = fb;
		
		if ( fa * fs < 0) {
			b = ss;
			fb = fs;
		} else {
			a = ss;
			fa = fs;
		}
		
		if (std::abs(fa) < std::abs(fb)) {
			std::swap(a,b);
			std::swap(fa,fb);
		}
	}
	return -1.0*ps;
}


G4bool TsIRTConfiguration::Inside(G4ThreeVector p ) {
	G4double delta = 0.5 * 1e-3 * nm;
	G4double fDx = 0.5*um;
	G4double fDy = 0.5*um;
	G4double fDz = 0.5*um;
	G4double dist = std::max(std::max(std::abs(p.x())-fDx,
									  std::abs(p.y())-fDy),
							          std::abs(p.z())-fDz);
	if (dist > delta) return false;
	return true;
}


G4bool TsIRTConfiguration::MakeReaction(std::unordered_map<G4int,TsMolecule> &initialSpecies, G4int& speciesIndex,
										std::unordered_map<G4int, std::unordered_map<G4int, std::unordered_map<G4int, std::unordered_map<G4int,G4bool>>>> &spaceBinned,
										G4int NX, G4int NY, G4int NZ, G4double XMin, G4double XMax, G4double YMin,
										G4double YMax, G4double ZMin, G4double ZMax,
										std::map<G4int, std::map<G4int, G4int>> &theGvalue,
										std::vector<G4double> timeSteps,
										G4int iM, G4int indexOfReaction, G4double irt,
										std::unordered_map<G4int,G4bool> &used, std::vector<G4int>& prods) {
	
	G4ThreeVector positions = initialSpecies[iM].position;
	std::vector<G4int> products = (GetReaction(indexOfReaction)).products;
	G4int tBin = fUtils->FindBin(irt, timeSteps);
	if ( tBin < 0 ) return false;
	
	for ( size_t u = 0; u < products.size(); u++ ) {
		TsMolecule aProd;
		aProd.id = products[u];
		aProd.position = positions;
		aProd.time = irt;
		aProd.trackID = 0;
		aProd.reacted = false;
		aProd.isDNA = false;
		aProd.isNew = true;
		if ( products[u] == 1 || products[u] == 5)
			aProd.spin = G4UniformRand() > 0.5 ? 1 : 0;
		else
			aProd.spin = -1;
		
		G4int i = fUtils->FindBin(NX, XMin, XMax, positions.x());
		G4int j = fUtils->FindBin(NY, YMin, YMax, positions.y());
		G4int k = fUtils->FindBin(NZ, ZMin, ZMax, positions.z());
		
		spaceBinned[i][j][k][speciesIndex] = true;
		initialSpecies[speciesIndex] = aProd;
		used[speciesIndex] = false;
		prods.push_back(speciesIndex);
		speciesIndex++;
		
		for ( int ti = tBin; ti < (int)timeSteps.size(); ti++ ) {
			theGvalue[aProd.id][ti]++;
		}
	}
	for ( int ti = tBin; ti < (int)timeSteps.size(); ti++ )
		theGvalue[initialSpecies[iM].id][ti]--;
	
	initialSpecies[iM].reacted = true;
	return true;
}


G4bool TsIRTConfiguration::MakeReaction(std::unordered_map<G4int,TsMolecule> &initialSpecies, G4int& speciesIndex,
										std::unordered_map<G4int, std::unordered_map<G4int, std::unordered_map<G4int, std::unordered_map<G4int,G4bool>>>> &spaceBinned,
										G4int NX, G4int NY, G4int NZ, G4double XMin, G4double XMax, G4double YMin,
										G4double YMax, G4double ZMin, G4double ZMax,
										std::map<G4int, std::map<G4int, G4int>> &theGvalue,
										std::vector<G4double> timeSteps,
										G4int iM, G4int jM, G4int indexOfReaction, G4double irt,
										G4double probabilityOfReaction, std::unordered_map<G4int,G4bool> &used, std::vector<G4int>& prods) {
	
	if ( G4UniformRand() < probabilityOfReaction ) {
		std::vector<G4ThreeVector> positions = GetPositionOfProducts(initialSpecies[iM], initialSpecies[jM], indexOfReaction);
		std::vector<G4int> products = (GetReaction(indexOfReaction)).products;
		G4int tBin = fUtils->FindBin(irt, timeSteps);
		if ( tBin < 0 ) return false;
		
		for ( size_t u = 0; u < products.size(); u++ ) {
			TsMolecule aProd;
			aProd.id = products[u];
			aProd.position = positions[u];
			aProd.time = irt;
			aProd.trackID = 0;
			aProd.reacted = false;
			aProd.isDNA = false;
			aProd.isNew = true;
			if ( products[u] == 1 || products[u] == 5)
				aProd.spin = G4UniformRand() > 0.5 ? 1 : 0;
			else
				aProd.spin = -1;
			
			G4int i = fUtils->FindBin(NX, XMin, XMax, positions[u].x());
			G4int j = fUtils->FindBin(NY, YMin, YMax, positions[u].y());
			G4int k = fUtils->FindBin(NZ, ZMin, ZMax, positions[u].z());
			
			spaceBinned[i][j][k][speciesIndex] = true;
			initialSpecies[speciesIndex] = aProd;
			used[speciesIndex] = false;
			prods.push_back(speciesIndex);
			speciesIndex++;
			
			for ( int ti = tBin; ti < (int)timeSteps.size(); ti++ ) {
				theGvalue[aProd.id][ti]++;
			}
		}
		for ( int ti = tBin; ti < (int)timeSteps.size(); ti++ ) {
			theGvalue[initialSpecies[iM].id][ti]--;
			theGvalue[initialSpecies[jM].id][ti]--;
		}
		initialSpecies[iM].reacted = true;
		initialSpecies[jM].reacted = true;
		return true;
	} else {
		return false;
	}
}


G4bool TsIRTConfiguration::MakeReaction(std::unordered_map<G4int,TsMolecule> &initialSpecies, G4int& speciesIndex,
										std::unordered_map<G4int, std::unordered_map<G4int, std::unordered_map<G4int, std::unordered_map<G4int,G4bool>>>> &spaceBinned,
										G4int NX, G4int NY, G4int NZ, G4double XMin, G4double XMax, G4double YMin,
										G4double YMax, G4double ZMin, G4double ZMax,
										std::map<G4int, std::map<G4int, G4int>> &theGvalue,
										std::map<G4int, std::map<G4int, G4int>> &theGvalueInVolume,
										std::vector<G4double> timeSteps,
										G4int iM, G4int jM, G4int indexOfReaction, G4double irt,
										G4double probabilityOfReaction, std::unordered_map<G4int,G4bool> &used, std::vector<G4int>& prods) {
	
	if ( G4UniformRand() < probabilityOfReaction ) {
		std::vector<G4ThreeVector> positions = GetPositionOfProducts(initialSpecies[iM], initialSpecies[jM], indexOfReaction);
		std::vector<G4int> products = (GetReaction(indexOfReaction)).products;
		G4int tBin = fUtils->FindBin(irt, timeSteps);
		if ( tBin < 0 ) return false;
		
		G4bool inVolume = false;
		if ( Inside(initialSpecies[iM].position) || Inside(initialSpecies[jM].position))  // at least one specie is at scoring region
			inVolume = true;
		
		for ( size_t u = 0; u < products.size(); u++ ) {
			TsMolecule aProd;
			aProd.id = products[u];
			aProd.position = positions[u];
			aProd.time = irt;
			aProd.trackID = 0;
			aProd.reacted = false;
			aProd.isDNA = false;
			aProd.isNew = true;
			if ( products[u] == 1 || products[u] == 5)
				aProd.spin = G4UniformRand() > 0.5 ? 1 : 0;
			else
				aProd.spin = -1;
			
			G4int i = fUtils->FindBin(NX, XMin, XMax, positions[u].x());
			G4int j = fUtils->FindBin(NY, YMin, YMax, positions[u].y());
			G4int k = fUtils->FindBin(NZ, ZMin, ZMax, positions[u].z());
			
			spaceBinned[i][j][k][speciesIndex] = true;;
			initialSpecies[speciesIndex] = aProd;
			used[speciesIndex] = false;
			prods.push_back(speciesIndex);
			speciesIndex++;
			
			for ( int ti = tBin; ti < (int)timeSteps.size(); ti++ ) {
				theGvalue[aProd.id][ti]++;
				if ( inVolume )
					theGvalueInVolume[aProd.id][ti]++;
			}
		}
		for ( int ti = tBin; ti < (int)timeSteps.size(); ti++ ) {
			theGvalue[initialSpecies[iM].id][ti]--;
			theGvalue[initialSpecies[jM].id][ti]--;
			if ( inVolume ) {
				theGvalueInVolume[initialSpecies[iM].id][ti]--;
				theGvalueInVolume[initialSpecies[jM].id][ti]--;
			}
		}
		initialSpecies[iM].reacted = true;
		initialSpecies[jM].reacted = true;
		return true;
	} else {
		return false;
	}
}


void TsIRTConfiguration::ScoreGvalue(std::vector<TsMolecule> &initialSpecies,
									 std::map<G4int, std::map<G4int, G4int>> &theGvalueInVolume,
									 std::vector<G4double> timeSteps,
									 G4int iM, G4int jM, G4int indexOfReaction, G4double irt) {
	
	std::vector<G4int> products = (GetReaction(indexOfReaction)).products;
	G4int tBin = fUtils->FindBin(irt, timeSteps);
	if ( tBin < 0 ) return;
	
	for ( size_t u = 0; u < products.size(); u++ ) {
		for ( int ti = tBin; ti < (int)timeSteps.size(); ti++ ) {
			theGvalueInVolume[products[u]][ti]++;
		}
	}
	
	for ( int ti = tBin; ti < (int)timeSteps.size(); ti++ ) {
		theGvalueInVolume[initialSpecies[iM].id][ti]--;
		theGvalueInVolume[initialSpecies[jM].id][ti]--;
	}
	
	return;
}


std::vector<G4double> TsIRTConfiguration::GetH2SO4ComponentsConcentrationPH(G4double pH) {
	G4double Ka1 = pow(10,3);
	G4double Ka2 = pow(10,-1.987);
	G4double _Kw = 1E-14;
	
	G4double _H_pos = pow(10,-pH);
	G4double _SO_4    = (_H_pos - (_Kw / _H_pos)) / ((_H_pos / Ka2) + 2);
	G4double _HSO_4   = _SO_4 * _H_pos / Ka2;
	G4double _OH_me   = _Kw / _H_pos;
	G4double _H_2SO_4 = _HSO_4 * _H_pos / Ka1;
	G4double H_Poly = 0;
	
	std::vector<G4double> Results = {_H_pos, _HSO_4, _SO_4, _OH_me, _H_2SO_4, H_Poly};
	
	return Results;
}


std::vector<G4double> TsIRTConfiguration::GetH2SO4ComponentsConcentrationP(G4double Concentration) {
	G4double Ka1 = pow(10,3);
	G4double Ka2 = pow(10,-1.987);
	G4double _C = Concentration;
	G4double _Kw = 1E-14;
	
	G4double _p4 = 1;
	G4double _p3 = Ka1;
	G4double _p2 = (Ka1 * Ka2) - (_C * Ka1) - _Kw;
	G4double _p1 = - ((Ka1 * _Kw) + (2 * _C * Ka1 * Ka2));
	G4double _p0 = - (_Kw * Ka1 * Ka2);
	
	std::vector<G4double> P = {_p4, _p3, _p2, _p1, _p0};
	
	std::vector<G4double> Result = fUtils->GetRoots(4, P);
	
	G4double _H_pos = -10;
	
	for (size_t i = 0; i < Result.size(); i++) {
		if ((Result[i] > 0 and Result[i] < 3.0*_C) and i % 2 != 1) {
			if (Result[i + 1] == 0) {
				_H_pos = Result[i];
			}
		}
	}
	
	if (_H_pos == -10)
		return {0, 0, 0, 0, 0, 0};
	
	G4double _SO_4    = (_H_pos - (_Kw / _H_pos)) / ((_H_pos / Ka2) + 2);
	G4double _HSO_4   = _SO_4 * _H_pos / Ka2;
	G4double _OH_me   = _Kw / _H_pos;
	G4double _H_2SO_4 = _HSO_4 * _H_pos / Ka1;
	
	G4double H_Poly = (P[0] * pow(_H_pos,4)) + (P[1] * pow(_H_pos,3)) +
	(P[2] * pow(_H_pos,2)) + (P[3] * _H_pos) + (P[4]);
	
	std::vector<G4double> Results = {_H_pos, _HSO_4, _SO_4, _OH_me, _H_2SO_4, H_Poly};
	
	return Results;
}


G4double TsIRTConfiguration::GetIonicStrength(std::vector<G4double> Components) {
	G4double I = .5 * ((Components[0]) + (Components[1]) + (Components[2] * 4) + (Components[3]));
	return I;
}


G4double TsIRTConfiguration::IonicRate(G4double IonicStrength, G4double Rate, G4int Charge1, G4int Charge2) {
	G4double a = 0.15;
	G4double b = a * Charge1 * Charge2;
	G4double scaleFactor = pow(10,((1.02 * Charge1 * Charge2 * ((sqrt(IonicStrength)) / (1 + sqrt(IonicStrength))))  -  (2 * b * IonicStrength)));
	G4double K = Rate * scaleFactor;
	return K;
}


G4double TsIRTConfiguration::IonicRate(G4double IonicStrength, TsMolecularReaction Reaction) {
	G4double a = 0.15;
	G4int Charge1 = fMoleculesDefinition[Reaction.reactorA].charge;
	G4int Charge2 = fMoleculesDefinition[Reaction.reactorB].charge;
	G4int BackGround = Reaction.reactionType;
	G4double Rate = 0.0;
	if ( BackGround == 6 ) {
		Rate = Reaction.scavengingCapacity;
		Rate /= 1E-7;
	} else {
		Rate = Reaction.kobs;
	}
	
	G4double b = a * Charge1 * Charge2;
	G4double scaleFactor = pow(10,((1.02 * Charge1 * Charge2 * ((sqrt(IonicStrength)) / (1 + sqrt(IonicStrength))))  -  (2 * b * IonicStrength)));
	G4double K = Rate * scaleFactor;
	return K;
}
