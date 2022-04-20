// Scorer for DNADamageStepByStep
// ********************************************************************
// *																  *
// * This file is part of the TOPAS-nBio extensions to the			  *
// *   TOPAS Simulation Toolkit.									  *
// * The TOPAS-nBio extensions are freely available under the license *
// *   agreement set forth at: https://topas-nbio.readthedocs.io/	  *
// *																  *
// ********************************************************************
//
// Authors: Alejandro Bertolet, Jan Schuemann

#include "TsScoreDNADamageSBS.hh"

#include "G4SystemOfUnits.hh"
#include "G4VTouchable.hh"

#include "G4VProcess.hh"
#include "G4Molecule.hh"

TsScoreDNADamageSBS::TsScoreDNADamageSBS(TsParameterManager* pM, TsMaterialManager* mM, TsGeometryManager* gM, TsScoringManager* scM, TsExtensionManager* eM,
										G4String scorerName, G4String quantity, G4String outFileName, G4bool isSubScorer)
							: TsVNtupleScorer(pM, mM, gM, scM, eM, scorerName, quantity, outFileName, isSubScorer)
{
	// Initialize physical quantities
	fEdep = 0;
	fTrackAveragedLET = 0;
	fDoseInThisExposure = 0;
	fExposureID = 0;

	// Initialize quantification of damage
	fNumSB = 0; fNumSBDirect = 0; fNumSBQuasiDirect = 0; fNumSBIndirect = 0;
	fNumSSB = 0; fNumSSBDirect = 0; fNumSSBQuasiDirect = 0; fNumSSBIndirect = 0;
	fNumDSB = 0; fNumDSBDirect = 0; fNumDSBIndirect = 0; fNumDSBDirectIndirect = 0; fNumDSBDirectQuasiDirect = 0; fNumDSBQuasiDirectQuasiDirect = 0; fNumDSBIndirectQuasiDirect = 0;
	fNumBaseDamage = 0; fNumBaseDamageDirect = 0; fNumBaseDamageQuasiDirect = 0; fNumBaseDamageIndirect = 0;
	fYSB = 0; fYSSB = 0; fYDSB = 0;
	fYBaseDam = 0;

	//---------------
	// Get parameters
	//---------------
	fNumberOfHistoriesInRun = fPm->GetIntegerParameter(GetFullParmName("NumberOfHistoriesInRun"));

	// Output filename
	fOutFileName = fPm->GetStringParameter(GetFullParmName("OutputFile"));

	// Parameters for material filter
	fBasePairDepth = 0;
	if (fPm->ParameterExists(GetFullParmName("BasePairPositionAtGeometricHierarchy")))
		fBasePairDepth = fPm->GetIntegerParameter(GetFullParmName("BasePairPositionAtGeometricHierarchy"));
	G4String* strand1Materials;
	G4String* strand2Materials;
	G4int strand1Length, strand2Length;
	if (fPm->ParameterExists(GetFullParmName("Strand1MaterialNames")))
	{
		strand1Materials = fPm->GetStringVector(GetFullParmName("Strand1MaterialNames"));
		strand1Length = fPm->GetVectorLength(GetFullParmName("Strand1MaterialNames"));
	}
	else
		strand1Materials[0] = "G4_WATER";
	if (fPm->ParameterExists(GetFullParmName("Strand2MaterialNames")))
	{
		strand2Materials = fPm->GetStringVector(GetFullParmName("Strand2MaterialNames"));
		strand2Length = fPm->GetVectorLength(GetFullParmName("Strand2MaterialNames"));
	}
	else
		strand2Materials[0] = "G4_WATER";
	for (G4int i = 0; i < strand1Length; i++)
		fStrand1Materials.push_back(GetMaterial(strand1Materials[i]));
	for (G4int i = 0; i < strand2Length; i++)
		fStrand2Materials.push_back(GetMaterial(strand2Materials[i]));

	// Options for direct damage
	fDirectDamageThreshold = 11.75 * eV; // 17,5 eV for half-cylinder
	if (fPm->ParameterExists(GetFullParmName("DirectDamageThreshold")))
		fDirectDamageThreshold = fPm->GetBooleanParameter(GetFullParmName("DirectDamageThreshold"));
	fUseLinearProbabilityForDirectDamage = false;
	if ( fPm->ParameterExists(GetFullParmName("UseLinearProbabilityForDirectDamage")) )
		fUseLinearProbabilityForDirectDamage = fPm->GetBooleanParameter(GetFullParmName("UseLinearProbabilityForDirectDamage"));
	if (!fUseLinearProbabilityForDirectDamage && (fPm->ParameterExists(GetFullParmName("LowerLimitForLinearProbabilityFunction")) || fPm->ParameterExists(GetFullParmName("UpperLimitForLinearProbabilityFunction"))))
	{
		G4cerr << "TOPAS is exiting due to a serious error in scoring setup." << G4endl;
		G4cerr << GetFullParmName("UseLinearProbabilityForDirectDamage") << " is set to False and limits for linear probability function in DNA components were given." << G4endl;
		fPm->AbortSession(1);
	}
	else if (fUseLinearProbabilityForDirectDamage)
	{
		fLowerLimitLinearProbability = 5 * eV;
		if (fPm->ParameterExists(GetFullParmName("LowerLimitForLinearProbabilityFunction")))
			fLowerLimitLinearProbability = fPm->GetDoubleParameter(GetFullParmName("LowerLimitForLinearProbabilityFunction"), "Energy");
		fUpperLimitLinearProbability = 37.5 * eV;
		if (fPm->ParameterExists(GetFullParmName("UpperLimitForLinearProbabilityFunction")))
			fUpperLimitLinearProbability = fPm->GetDoubleParameter(GetFullParmName("UpperLimitForLinearProbabilityFunction"), "Energy");
	}
	// Options for quasi-direct damage
	fProbabilityOfChargeTransferFromHydrationShellToBackbone = 0.33333;
	if (fPm->ParameterExists(GetFullParmName("ProbabilityOfChargeTransferFromHydrationShellToBackbone")))
		fProbabilityOfChargeTransferFromHydrationShellToBackbone = fPm->GetUnitlessParameter(GetFullParmName("ProbabilityOfChargeTransferFromHydrationShellToBackbone"));

	// Options for indirect damage
	fAlwaysScavengeSpeciesInDNAComponents = false;
	if (fPm->ParameterExists(GetFullParmName("AlwaysScavengeSpeciesInDNAComponents")))
		fAlwaysScavengeSpeciesInDNAComponents = fPm->GetBooleanParameter(GetFullParmName("AlwaysScavengeSpeciesInDNAComponents"));

	if (fAlwaysScavengeSpeciesInDNAComponents && (fPm->ParameterExists(GetFullParmName("ProbabilityOfScavengingInBackbone")) || fPm->ParameterExists(GetFullParmName("ProbabilityOfScavengingInBase"))))
	{
		G4cerr << "TOPAS is exiting due to a serious error in scoring setup." << G4endl;
		G4cerr << GetFullParmName("AlwaysScavengeSpeciesInDNAComponents") << " is set to True and probabilities for scavenging in DNA components were given." << G4endl;
		fPm->AbortSession(1);
	}
	if (!fAlwaysScavengeSpeciesInDNAComponents)
	{
		fProbabilityOfScavengingInBackbone = 0.25; // 0.0585 for half-cylinder
		if (fPm->ParameterExists(GetFullParmName("ProbabilityOfScavengingInBackbone")))
			fProbabilityOfScavengingInBackbone = fPm->GetUnitlessParameter(GetFullParmName("ProbabilityOfScavengingInBackbone"));
		fProbabilityOfScavengingInBase = 1.0;
		if (fPm->ParameterExists(GetFullParmName("ProbabilityOfScavengingInBase")))
			fProbabilityOfScavengingInBase = fPm->GetUnitlessParameter(GetFullParmName("ProbabilityOfScavengingInBase"));
		fProbabilityOfDamageInBackbone = 0.55;
		if (fPm->ParameterExists(GetFullParmName("ProbabilityOfIndirectDamageToBackbone")))
			fProbabilityOfDamageInBackbone = fPm->GetUnitlessParameter(GetFullParmName("ProbabilityOfIndirectDamageToBackbone"));
		fProbabilityOfDamageInBase = 1.0;
		if (fPm->ParameterExists(GetFullParmName("ProbabilityOfIndirectDamageToBase")))
			fProbabilityOfDamageInBase = fPm->GetUnitlessParameter(GetFullParmName("ProbabilityOfIndirectDamageToBase"));
	}
	else
	{
		fProbabilityOfScavengingInBackbone = 1.0;
		fProbabilityOfScavengingInBase = 1.0;
		fProbabilityOfDamageInBackbone = 0.4;
		if (fPm->ParameterExists(GetFullParmName("ProbabilityOfIndirectDamageToBackbone")))
			fProbabilityOfDamageInBackbone = fPm->GetUnitlessParameter(GetFullParmName("ProbabilityOfIndirectDamageToBackbone"));
		fProbabilityOfDamageInBase = 0.4;
		if (fPm->ParameterExists(GetFullParmName("ProbabilityOfIndirectDamageToBase")))
			fProbabilityOfDamageInBase = fPm->GetUnitlessParameter(GetFullParmName("ProbabilityOfIndirectDamageToBase"));
	}

	fScavengeInHistones = true;
	if (fPm->ParameterExists(GetFullParmName("ScavengeInHistones")))
		fScavengeInHistones = fPm->GetBooleanParameter(GetFullParmName("ScavengeInHistones"));

	// Classify damage as SSBs and DSBs
	fScoreDSB = true;
	if (fPm->ParameterExists(GetFullParmName("ScoreDSBs")))
		fScoreDSB = fPm->GetBooleanParameter(GetFullParmName("ScoreDSBs"));
	fNumberOfBasePairsForDSB = 10;
	if (fPm->ParameterExists(GetFullParmName("MaximumBasePairDistanceToConsiderDSB")))
		fNumberOfBasePairsForDSB = fPm->GetIntegerParameter(GetFullParmName("MaximumBasePairDistanceToConsiderDSB"));

	// Classify damage as foci
	fScoreFoci = true;
	if (fPm->ParameterExists(GetFullParmName("ScoreNumberOfFoci")))
		fScoreFoci = fPm->GetBooleanParameter(GetFullParmName("ScoreNumberOfFoci"));
	fFociSize = 500 * nm;
	if (fPm->ParameterExists(GetFullParmName("FociSize")))
		fFociSize = fPm->GetDoubleParameter(GetFullParmName("FociSize"), "Length");

	// Considering fragments
	fExcludeShortFragments = false;
	if (fPm->ParameterExists(GetFullParmName("ExcludeShortFragments")))
		fExcludeShortFragments = fPm->GetBooleanParameter(GetFullParmName("ExcludeShortFragments"));
	if (!fExcludeShortFragments && (fPm->ParameterExists(GetFullParmName("LowerThresholdForFragmentDetection")) || fPm->ParameterExists(GetFullParmName("UpperThresholdForFragmentDetection"))))
	{
		G4cerr << "TOPAS is exiting due to a serious error in scoring setup." << G4endl;
		G4cerr << GetFullParmName("ExcludeShortFragments") << " is set to False and limits for fragment detection were given." << G4endl;
		fPm->AbortSession(1);
	}
	else if (fExcludeShortFragments)
	{
		fLowerThresholdForFragmentDetection = 0.0;
		if (fPm->ParameterExists(GetFullParmName("LowerThresholdForFragmentDetection")))
			fLowerThresholdForFragmentDetection = fPm->GetIntegerParameter(GetFullParmName("LowerThresholdForFragmentDetection"));
		fUpperThresholdForFragmentDetection = 3E8;
		if (fPm->ParameterExists(GetFullParmName("UpperThresholdForFragmentDetection")))
			fUpperThresholdForFragmentDetection = fPm->GetIntegerParameter(GetFullParmName("UpperThresholdForFragmentDetection"));
	}

	// Options for the output
	fWriteCSVWithExtensiveDamage = false;
	if (fPm->ParameterExists(GetFullParmName("WriteCSVOutputWithAllDamageSpecification")))
		fWriteCSVWithExtensiveDamage = fPm->GetBooleanParameter(GetFullParmName("WriteCSVOutputWithAllDamageSpecification"));
	fScoreDirectDamage = true;
	if (fPm->ParameterExists(GetFullParmName("ScoreDirectDamage")))
		fScoreDirectDamage = fPm->GetBooleanParameter(GetFullParmName("ScoreDirectDamage"));
	fScoreIndirectDamage = true;
	if (fPm->ParameterExists(GetFullParmName("ScoreIndirectDamage")))
		fScoreIndirectDamage = fPm->GetBooleanParameter(GetFullParmName("ScoreIndirectDamage"));
	fScoreQuasiDirectDamage = true;
	if (fPm->ParameterExists(GetFullParmName("ScoreQuasiDirectDamage")))
		fScoreQuasiDirectDamage = fPm->GetBooleanParameter(GetFullParmName("ScoreQuasiDirectDamage"));
	fScoreOnBases = true;
	if (fPm->ParameterExists(GetFullParmName("ScoreBaseDamages")))
		fScoreOnBases = fPm->GetBooleanParameter(GetFullParmName("ScoreBaseDamages"));
	fScoreOnBackbones = true;
	if (fPm->ParameterExists(GetFullParmName("ScoreBackboneDamages")))
		fScoreOnBackbones = fPm->GetBooleanParameter(GetFullParmName("ScoreBackboneDamages"));
	fBreakDownPerDamageOrigin = true;
	if (fPm->ParameterExists(GetFullParmName("BreakDownOutputPerDamageOrigin")))
		fBreakDownPerDamageOrigin = fPm->GetBooleanParameter(GetFullParmName("BreakDownOutputPerDamageOrigin"));

	// For SDD specification
	fDosePerExposure = 1 * gray;
	if (fPm->ParameterExists(GetFullParmName("DosePerExposure")))
		fDosePerExposure = fPm->GetDoubleParameter(GetFullParmName("DosePerExposure"), "Dose");
	fOnlyIncludeDSBinSDD = false;
	if (fPm->ParameterExists(GetFullParmName("IncludeDSBDamageSitesOnlyinSDD")))
		fOnlyIncludeDSBinSDD = fPm->GetBooleanParameter(GetFullParmName("IncludeDSBDamageSitesOnlyinSDD"));
	fWriteMinimalSDDOutput = false;
	if (fPm->ParameterExists(GetFullParmName("WriteMinimalSDDOutput")))
		fWriteMinimalSDDOutput = fPm->GetBooleanParameter(GetFullParmName("WriteMinimalSDDOutput"));
	fPrimaryParticle = "proton";
	if (fPm->ParameterExists(GetFullParmName("PrimaryParticle")))
		fWriteMinimalSDDOutput = fPm->GetStringParameter(GetFullParmName("PrimaryParticle"));

	// Parameters for the SDD header
	G4String author = "@";
	if ( fPm->ParameterExists(GetFullParmName("AuthorName")) )
		author = fPm->GetStringParameter(GetFullParmName("AuthorName"));
	G4String simulationDetails = "Sim details";
	if ( fPm->ParameterExists(GetFullParmName("SimulationDetails")) )
		simulationDetails = fPm->GetStringParameter(GetFullParmName("SimulationDetails"));
	G4String sourceDetails = "Source details";
	if ( fPm->ParameterExists(GetFullParmName("SourceDetails")) )
		sourceDetails = fPm->GetStringParameter(GetFullParmName("SourceDetails"));
	G4int sourceType = 1;
	if ( fPm->ParameterExists(GetFullParmName("SourceType")) )
		sourceType = fPm->GetIntegerParameter(GetFullParmName("SourceType"));
	G4double meanEnergy = 0.0;
	if ( fPm->ParameterExists(GetFullParmName("MeanEnergy")))
		meanEnergy = fPm->GetDoubleParameter(GetFullParmName("MeanEnergy"), "Energy");
	G4String energyDist = "M, 0";
	if ( fPm->ParameterExists(GetFullParmName("EnergyDistribution")) )
		energyDist = fPm->GetStringParameter(GetFullParmName("EnergyDistribution"));
	G4String irrTarget = "";
	if ( fPm->ParameterExists(GetFullParmName("IrradiationTarget")) )
		irrTarget = fPm->GetStringParameter(GetFullParmName("IrradiationTarget"));
	G4String cellCycle = "0";
	if ( fPm->ParameterExists(GetFullParmName("CellCycleStage")))
		cellCycle = fPm->GetStringParameter(GetFullParmName("CellCycleStage"));
	G4String DNAStructure = "0, 1";
	if ( fPm->ParameterExists(GetFullParmName("DNAStructure")))
		cellCycle = fPm->GetStringParameter(GetFullParmName("DNAStructure"));
	G4int inVitroOrInVivo = 0;
	if ( fPm->ParameterExists(GetFullParmName("InVitroOrInVivo")) )
		inVitroOrInVivo = fPm->GetIntegerParameter(GetFullParmName("InVitroOrInVivo"));
	G4String proliferationStatus = "1";
	if ( fPm->ParameterExists(GetFullParmName("ProliferationStatus")))
		proliferationStatus = fPm->GetStringParameter(GetFullParmName("ProliferationStatus"));
	G4String microenvironment = "20, 0.01";
	if ( fPm->ParameterExists(GetFullParmName("Microenvironment")))
		microenvironment = fPm->GetStringParameter(GetFullParmName("Microenvironment"));
	G4double time = 0;
	if ( fPm->ParameterExists(GetFullParmName("Time")))
		time = fPm->GetDoubleParameter(GetFullParmName("Time"), "Time");
	G4String addInfo = "";
	if ( fPm->ParameterExists(GetFullParmName("AdditionalInfo")))
		addInfo = fPm->GetStringParameter(GetFullParmName("AdditionalInfo"));

	// =============================
	//       PRINT PARAMETERS
	// =============================
	G4cout << "*********************************************************************************" << G4endl;
	if (fScoreDirectDamage)
	{
		G4cout << "DIRECT DAMAGE" << G4endl;
		G4cout << "--------------" << G4endl;
		if (!fUseLinearProbabilityForDirectDamage)
			G4cout << "Single energy threshold for damage = " << fDirectDamageThreshold/eV << " eV" << G4endl;
		else
			G4cout << "Linearly increasing probability of damage from 0 to 1 from " << fLowerLimitLinearProbability/eV << " eV to " << fUpperLimitLinearProbability/eV << " eV" << G4endl;
	}
	if (fScoreQuasiDirectDamage)
	{
		G4cout << "QUASI-DIRECT DAMAGE" << G4endl;
		G4cout << "-------------------" << G4endl;
		G4cout << "Probability for charge transfer from hydration shell to backbone to produce strand break = " << fProbabilityOfChargeTransferFromHydrationShellToBackbone << G4endl;
	}
	if (fScoreIndirectDamage)
	{
		G4cout << "INDIRECT DAMAGE" << G4endl;
		G4cout << "-------------------" << G4endl;
		if (fAlwaysScavengeSpeciesInDNAComponents)
			G4cout << "Probability of scavenging in bases = 1.0; Probability of scavenging in backbones = 1.0" << G4endl;
		else
			G4cout << "Probability of scavenging in bases = " << fProbabilityOfScavengingInBase << "; Probability of scavenging in backbones = " << fProbabilityOfScavengingInBackbone << G4endl;
		G4cout << "Probability of base damage after scavenging = " << fProbabilityOfDamageInBase << "; Probability of strand break after scavenging = " << fProbabilityOfDamageInBackbone << G4endl;
	}
	G4cout << "*********************************************************************************" << G4endl;

	// Register variables in nTuple
	fNtuple->RegisterColumnD(&fEdep, "Energy_imparted_per_event", "keV");
	fNtuple->RegisterColumnD(&fDoseInThisExposure, "Dose_per_event", "Gy");
	fNtuple->RegisterColumnD(&fTrackAveragedLET, "LET_kev/um", "");
	if (fScoreDSB)
	{
		if (fScoreOnBackbones)
		{
			fNtuple->RegisterColumnD(&fYDSB, "DSB/Gy/Gbp", "");
			fNtuple->RegisterColumnD(&fYSSB, "SSB/Gy/Gbp", "");
			fNtuple->RegisterColumnD(&fYSB, "SB/Gy/Gbp", "");
		}
		if (fScoreOnBases)
			fNtuple->RegisterColumnD(&fYBaseDam, "BD/Gy/Gbp", "");
		if (fBreakDownPerDamageOrigin)
		{
			if (fScoreOnBackbones)
			{
				fNtuple->RegisterColumnI(&fNumDSB, "DSBs");
				if (fScoreDirectDamage) fNtuple->RegisterColumnI(&fNumDSBDirect, "DSBs_Direct");
				if (fScoreIndirectDamage) fNtuple->RegisterColumnI(&fNumDSBIndirect, "DSBs_Indirect");
				if (fScoreDirectDamage && fScoreIndirectDamage) fNtuple->RegisterColumnI(&fNumDSBDirectIndirect, "DSBs_Hybrid");
				if (fScoreDirectDamage && fScoreQuasiDirectDamage) fNtuple->RegisterColumnI(&fNumDSBDirectQuasiDirect, "DSBs_Direct_WithOneQuasiDirect");
				if (fScoreQuasiDirectDamage) fNtuple->RegisterColumnI(&fNumDSBQuasiDirectQuasiDirect, "DSBs_Direct_WithBothQuasiDirect");
				if (fScoreIndirectDamage && fScoreQuasiDirectDamage) fNtuple->RegisterColumnI(&fNumDSBIndirectQuasiDirect, "DSBs_Hybrid_WithOneQuasiDirect");
				fNtuple->RegisterColumnI(&fNumSSB, "SSBs");
				if (fScoreDirectDamage) fNtuple->RegisterColumnI(&fNumSSBDirect, "SSBs_Direct");
				if (fScoreQuasiDirectDamage) fNtuple->RegisterColumnI(&fNumSSBQuasiDirect, "SSBs_QuasiDirect");
				if (fScoreIndirectDamage) fNtuple->RegisterColumnI(&fNumSSBIndirect, "SSBs_Indirect");
				fNtuple->RegisterColumnI(&fNumSB, "SBs");
				if (fScoreDirectDamage) fNtuple->RegisterColumnI(&fNumSSBDirect, "SBs_Direct");
				if (fScoreQuasiDirectDamage) fNtuple->RegisterColumnI(&fNumSSBQuasiDirect, "SBs_QuasiDirect");
				if (fScoreIndirectDamage) fNtuple->RegisterColumnI(&fNumSSBIndirect, "SBs_Indirect");
			}
			if (fScoreOnBases)
			{
				fNtuple->RegisterColumnI(&fNumBaseDamage, "BDs");
				if (fScoreDirectDamage) fNtuple->RegisterColumnI(&fNumBaseDamageDirect, "BDs_Direct");
				if (fScoreQuasiDirectDamage) fNtuple->RegisterColumnI(&fNumBaseDamageQuasiDirect, "BDs_QuasiDirect");
				if (fScoreIndirectDamage) fNtuple->RegisterColumnI(&fNumBaseDamageIndirect, "BDs_Indirect");
			}
		}
	}
	if (fScoreFoci)
		fNtuple->RegisterColumnI(&fNumFoci, "Foci");
}

TsScoreDNADamageSBS::~TsScoreDNADamageSBS() {}

G4bool TsScoreDNADamageSBS::ProcessHits(G4Step* aStep, G4TouchableHistory*)
{
	if (!fIsActive)
	{
		fSkippedWhileInactive++;
		return false;
	}

	// Gets position
	G4ThreeVector pos = aStep->GetPreStepPoint()->GetPosition();

	// Accumulates energy for this event
	G4double edep = aStep->GetTotalEnergyDeposit();
	fEdep += edep;

	// Increments the number of steps this track has made
	G4int trackID = aStep->GetTrack()->GetTrackID();
	fTrackSteps[trackID] += 1;

	// Gets current volume and last volume for this track. Updates new volume for this track
	G4TouchableHistory* touchable = (G4TouchableHistory*)(aStep->GetPreStepPoint()->GetTouchable());
	G4String volumeName = touchable->GetVolume(fBasePairDepth)->GetName();
	G4bool enteringInNewVolume = false;
	if (volumeName != fTrackLastVolume[trackID])
		enteringInNewVolume = true;
	fTrackLastVolume[trackID] = volumeName;

	// Checks materials for determining if step is happening in DNA componets
	G4Material* material = aStep->GetPreStepPoint()->GetMaterial();
	G4bool materialMatched = false;
	for (G4int i = 0; i < fStrand1Materials.size(); i++)
	{
		if (material == fStrand1Materials[i])
		{
			materialMatched = true;
			break;
		}
	}
	if (materialMatched == 0)
	{
		for (G4int i = 0; i < fStrand2Materials.size(); i++)
		{
			if (material == fStrand2Materials[i])
			{
				materialMatched = true;
				break;
			}
		}
	}
	// Goes on only if DNA materials has been matched
	if (materialMatched)
	{
		// Gets bp ID
		G4int bpID = touchable->GetVolume(fBasePairDepth)->GetCopyNo();
		// Sets bpID to -1 if histone is touched
		if (strstr(volumeName, "Histone") != NULL)
			bpID = -1;

		// Gets strand number and DNA component ID (see header file for component IDs)
		G4int componentID = -1; G4int strandID = -1;
		if (strstr(volumeName, "Base1") != NULL) { componentID = base; strandID = 1; }
		else if (strstr(volumeName, "Base2") != NULL) { componentID = base; strandID = 2; }
		else if (strstr(volumeName, "Backbone1") != NULL) { componentID = backbone; strandID = 1; }
		else if (strstr(volumeName, "Backbone2") != NULL) { componentID = backbone; strandID = 2; }
		else if (strstr(volumeName, "HydrationShell1") != NULL) { componentID = hydrationshell; strandID = 1; }
		else if (strstr(volumeName, "HydrationShell2") != NULL) { componentID = hydrationshell; strandID = 2; }
		else if (strstr(volumeName, "Histone") != NULL) { componentID = histone; }

		// Gets particle and process info
		G4ParticleDefinition* particle = aStep->GetTrack()->GetParticleDefinition();
		G4String particleName = particle->GetParticleName();
		G4String processName = (G4String)aStep->GetPostStepPoint()->GetProcessDefinedStep()->GetProcessName();

		// Sets hit information
		TsHitInDNA* hit = new TsHitInDNA();
		hit->SetEventID(GetEventID() + fNumberOfHistoriesInRun * GetRunID());
		hit->SetEdep(edep);
		hit->SetParticleName(particleName);
		hit->SetBasePairID(bpID);
		hit->SetStrandNumber(strandID);
		hit->SetDNAComponentID(componentID);
		hit->SetPosition(pos);

		// Adds direct damage
		if (trackID >= 0 && edep > 0 && fScoreDirectDamage)
		{
			fHits.push_back(hit);
			return true;
		}
		// Adds quasi-direct damage
		if (componentID == hydrationshell && strstr(processName, "Ionisation") != NULL)
		{
			fHits.push_back(hit);
			return true;
		}
		// Adds indirect damage
		if (trackID < 0 && fScoreIndirectDamage)
		{
			G4String speciesName = GetMolecule(aStep->GetTrack())->GetName();
			G4bool isSpeciesToKill = (speciesName == "OH^0" || speciesName == "e_aq^-1" || speciesName == "H^0");
			G4bool isHydroxil = (speciesName == "OH^0");
			G4bool isHydElectron =  (speciesName == "e^aq^-1");
			// Kills all species generated inside DNA volumes except for the hydration shell
			if (fTrackSteps[trackID] == 1 && componentID != hydrationshell)
			{
				aStep->GetTrack()->SetTrackStatus(fStopAndKill);
				delete hit;
				return false;
			}
			// Makes damage to bases. OH and e_aq induce damage to bases. Only one step per species is considered (otherwise all would end up reacting), so we use the enteringInNewVolume flag
			else if ((isHydroxil || isHydElectron) && componentID == base && enteringInNewVolume)
			{
				hit->SetEdep(-0.001 * eV);
				G4bool scavenged = false;
				if (fAlwaysScavengeSpeciesInDNAComponents)
					scavenged = true;
				else
					if (G4UniformRand() < fProbabilityOfScavengingInBase) scavenged = true;
				if (scavenged)
				{
					aStep->GetTrack()->SetTrackStatus(fStopAndKill);
					if (G4UniformRand() < fProbabilityOfDamageInBase)
					{
						fHits.push_back(hit);
						return true;
					}
					else
					{
						delete hit;
						return false;
					}
				}
			}
			// Makes damage to backbones. Only OH induces damage to backbones. Only one step per species is considered (otherwise all would end up reacting), so we use the enteringInNewVolume flag
			else if (isHydroxil && componentID == backbone && enteringInNewVolume)
			{
				hit->SetEdep(-0.001 * eV);
				G4bool scavenged = false;
				if (fAlwaysScavengeSpeciesInDNAComponents)
					scavenged = true;
				else
					if (G4UniformRand() < fProbabilityOfScavengingInBackbone) scavenged = true;
				if (scavenged)
				{
					aStep->GetTrack()->SetTrackStatus(fStopAndKill);
					if (G4UniformRand() < fProbabilityOfDamageInBackbone)
					{
						fHits.push_back(hit);
						return true;
					}
					else
					{
						delete hit;
						return false;
					}
				}
			}
			// Scavenge species by histones
			else if (isSpeciesToKill && fScavengeInHistones && componentID == histone)
			{
				aStep->GetTrack()->SetTrackStatus(fStopAndKill);
				delete hit;
				return false;
			}
		}
		delete hit;
	}
	return false;
}

void TsScoreDNADamageSBS::AccumulateEvent()
{
	fCollectionsOfHits.push_back(fHits);
	fHits.clear();
	G4double edep = fEdep;
	fEventsEdep.push_back(edep);
	fEdep = 0.;
}

void TsScoreDNADamageSBS::UserHookForEndOfRun()
{
	// Print info
	G4cout << "--------------------" << G4endl;
	G4cout << "Number of events comprised in Run: " << fCollectionsOfHits.size() << G4endl;
	G4int numberOfLesions = 0;
	for (G4int i = 0; i < fCollectionsOfHits.size(); i++)
	{
		numberOfLesions += Analyze(fCollectionsOfHits[i], i);
		fNtuple->Fill();
	}
	fCollectionsOfHits.clear();
	fEventsEdep.clear();
}

void TsScoreDNADamageSBS::AbsorbResultsFromWorkerScorer(TsVScorer* workerScorer)
{
	TsVNtupleScorer::AbsorbResultsFromWorkerScorer(workerScorer);
	TsScoreDNADamageSBS* workerMTScorer = dynamic_cast<TsScoreDNADamageSBS*>(workerScorer);

	for(G4int i=0; i < workerMTScorer->fCollectionsOfHits.size(); i++)
		fCollectionsOfHits.push_back(workerMTScorer->fCollectionsOfHits[i]);
	workerMTScorer->fCollectionsOfHits.clear();

	for(G4int i=0; i < workerMTScorer->fEventsEdep.size(); i++)
		fEventsEdep.push_back(workerMTScorer->fEventsEdep[i]);
	workerMTScorer->fEventsEdep.clear();
}
