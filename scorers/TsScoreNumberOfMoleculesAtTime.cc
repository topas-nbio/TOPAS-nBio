// Scorer for TsSBSMoleculesAtTime
//
// ********************************************************************
// *                                                                  *
// * This file is part of the TOPAS-nBio extensions to the            *
// *   TOPAS Simulation Toolkit.                                      *
// * The TOPAS-nBio extensions are freely available under the license *
// *   agreement set forth at: https://topas-nbio.readthedocs.io/     *
// *                                                                  *
// ********************************************************************
//

#include "TsScoreNumberOfMoleculesAtTime.hh"

#include "TsChemistryManager.hh"

#include "G4SystemOfUnits.hh"
#include "G4DNAChemistryManager.hh"
#include "G4MoleculeCounter.hh"
#include "G4MolecularConfiguration.hh"
#include "G4Molecule.hh"
#include "G4RunManager.hh"

TsScoreNumberOfMoleculesAtTime::TsScoreNumberOfMoleculesAtTime(TsParameterManager *pM, TsMaterialManager *mM, TsGeometryManager *gM, TsScoringManager *scM, TsExtensionManager *eM,
                                                               G4String scorerName, G4String quantity, G4String outFileName, G4bool isSubScorer)
    : TsVNtupleScorer(pM, mM, gM, scM, eM, scorerName, quantity, outFileName, isSubScorer),
      fPm(pM)
{
    fNtuple->RegisterColumnI(&fNumberOfMolecules, "Total number of molecules at time");
    fNtuple->RegisterColumnF(&fMeanNumberOfMolecules, "Mean number of molecules per event at time", "");
    fNtuple->RegisterColumnF(&fTime, "Time", "picosecond");
    fNtuple->RegisterColumnS(&fMoleculeName, "MoleculeName");

    G4int tBins = 10;
    if (fPm->ParameterExists(GetFullParmName("TimeBins")))
        tBins = fPm->GetIntegerParameter(GetFullParmName("TimeBins"));

    G4double tBinMin = 1.0 * ps;
    if (fPm->ParameterExists(GetFullParmName("TimeBinMin")))
        tBinMin = fPm->GetDoubleParameter(GetFullParmName("TimeBinMin"), "Time");

    G4double tBinMax = 1.0 * us;
    if (fPm->ParameterExists(GetFullParmName("TimeBinMax")))
        tBinMax = fPm->GetDoubleParameter(GetFullParmName("TimeBinMax"), "Time");

    G4bool tBinLog = true;
    if (fPm->ParameterExists(GetFullParmName("TimeBinLog")))
        tBinLog = fPm->GetBooleanParameter(GetFullParmName("TimeBinLog"));

    if (!tBinLog)
    {
        for (int i = 0; i < tBins; i++)
            fTimesToRecord.push_back(tBinMin + i * (tBinMax - tBinMin) / tBins);
    }
    else
    {
        if (tBinMin < 1.0 * ps)
        {
            G4cerr << "Topas is exiting due to a serious error in scoring." << G4endl;
            G4cerr << "The scorer named: " << GetName() << " has TimeBinMin lower than 1 ps" << G4endl;
            fPm->AbortSession(1);
        }

        G4double logXMin = std::log10(tBinMin);
        G4double logXMax = std::log10(tBinMax);
        for (int i = 0; i < tBins; i++)
            fTimesToRecord.push_back(std::pow(10, logXMin + i * (logXMax - logXMin) / (tBins - 1)));
    }

    fEnergyLossKill = 0.0;
    if (fPm->ParameterExists(GetFullParmName("KillPrimaryIfEnergyLossExceeds")))
        fEnergyLossKill = fPm->GetDoubleParameter(GetFullParmName("KillPrimaryIfEnergyLossExceeds"), "Energy");

    fEnergyLossAbort = 0.0;
    if (fPm->ParameterExists(GetFullParmName("AbortEventIfPrimaryEnergyLossExceeds")))
        fEnergyLossAbort = fPm->GetDoubleParameter(GetFullParmName("AbortEventIfPrimaryEnergyLossExceeds"), "Energy");

    fMaximumTrackLength = 0.0;
    if (fPm->ParameterExists(GetFullParmName("KillPrimaryBasedOnTrackLength")))
        fMaximumTrackLength = fPm->GetDoubleParameter(GetFullParmName("KillPrimaryBasedOnTrackLength"), "Length");

    fNbOfScoredEvents = 0;

    fMoleculeCounterName = TsChemistryManager::GetDefaultFixedPrecisionCounterName();

    // Ensure the MoleculeCounterManager is configured appropriately for this counter
    if (G4MoleculeCounterManager::Instance()->GetResetCountersBeforeEvent() == false)
    {
        // must reset before each event, otherwise G-values will be off
        G4cerr << "TOPAS is exiting due to a serious error." << G4endl;
        G4cerr << "TsSBSScoreGValue requires the molecule counters to be reset before each event!" << G4endl;
        fPm->AbortSession(1);
    }
}

G4bool TsScoreNumberOfMoleculesAtTime::ProcessHits(G4Step *aStep, G4TouchableHistory *)
{
    if (!fIsActive)
    {
        fSkippedWhileInactive++;
        return false;
    }
    if (-1 < aStep->GetTrack()->GetTrackID())
    { // physical tracks
        if (aStep->GetTrack()->GetParentID() == 0)
        {
            G4double kineticE = aStep->GetPostStepPoint()->GetKineticEnergy();
            G4double eLoss = aStep->GetPreStepPoint()->GetKineticEnergy() - kineticE;
            fEnergyLoss += eLoss;
            fTotalTrackLength += aStep->GetStepLength();

            if (0 < fEnergyLossAbort && fEnergyLoss > fEnergyLossAbort)
            {
                std::cout << " -- Aborting event " << GetEventID() << std::endl;
                G4RunManager::GetRunManager()->AbortEvent();
                return true;
            }

            if (0 < fEnergyLossKill && fEnergyLoss >= fEnergyLossKill)
            {
                std::cout << " -- Killing primary based on maximum accumulated energ loss." << std::endl;
                aStep->GetTrack()->SetTrackStatus(fStopAndKill);
                return true;
            }

            if (0 < fMaximumTrackLength && fTotalTrackLength >= fMaximumTrackLength)
            {
                aStep->GetTrack()->SetTrackStatus(fStopAndKill);
                std::cout << " -- Killing primary based on maximum accumulated track-length." << std::endl;
                return true;
            }
        }

        if (aStep->GetTotalEnergyDeposit() > 0)
            return true;
    }

    return false;
}

void TsScoreNumberOfMoleculesAtTime::AccumulateEvent()
{

    if (!G4EventManager::GetEventManager()->GetConstCurrentEvent()->IsAborted())
    {
        auto counters = G4MoleculeCounterManager::Instance()->GetMoleculeCounters(fMoleculeCounterName);

        G4MoleculeCounter *counter = nullptr;
        if (counters.rbegin() != counters.rend())
        {
            counter = const_cast<G4MoleculeCounter *>(dynamic_cast<const G4MoleculeCounter *>(*counters.begin()));
            if (counter == nullptr)
                G4Exception("ScoreMoleculeCounter::EndOfEvent", "SCOREMOLCOUNT", FatalException, "Molecule counter has wrong type!");
        }
        else
        {
            G4Exception("ScoreMoleculeCounter::EndOfEvent", "SCOREMOLCOUNT", FatalException, "No molecule counter with given name found!");
        }

        for (auto const &entry : counter->GetCounterMap())
        {
            auto key = entry.first.GetMolecule()->GetName();
            for (auto const &time : fTimesToRecord)
            {
                G4int nbOfMoleculesAtTime = counter->GetNbMoleculesAtTime(entry.first, time);
                fMoleculeCountPerIndexPerTime[key][time] += nbOfMoleculesAtTime;
            }
        }

        fNbOfScoredEvents++;
    }

    fEnergyLoss = 0.0;
    fTotalTrackLength = 0.0;
}

void TsScoreNumberOfMoleculesAtTime::AbsorbResultsFromWorkerScorer(TsVScorer *workerScorer)
{
    TsVNtupleScorer::AbsorbResultsFromWorkerScorer(workerScorer);

    TsScoreNumberOfMoleculesAtTime *worker = dynamic_cast<TsScoreNumberOfMoleculesAtTime *>(workerScorer);

    fNbOfScoredEvents += worker->fNbOfScoredEvents;

    for (auto const &workerEntry : worker->fMoleculeCountPerIndexPerTime)
    {
        auto emplacedPair = fMoleculeCountPerIndexPerTime.emplace(workerEntry.first, workerEntry.second);
        if (!emplacedPair.second)
        { // workerEntry.first already exists, emplacedPair.first == it
            auto it = emplacedPair.first;

            for (auto const &workerMoleculeCountWithTimeIter : workerEntry.second)
            {
                auto masterMoleculeCountWithTimeIter = it->second.find(workerMoleculeCountWithTimeIter.first);
                masterMoleculeCountWithTimeIter->second += workerMoleculeCountWithTimeIter.second;
            }
        }
    }

    worker->fMoleculeCountPerIndexPerTime.clear();
    worker->fNbOfScoredEvents = 0;
}

void TsScoreNumberOfMoleculesAtTime::Output()
{

    if (fNbOfScoredEvents == 0)
        return;

    for (auto const &wIter : fMoleculeCountPerIndexPerTime)
    {
        for (auto const &iter : wIter.second)
        {
            fNumberOfMolecules = iter.second;
            fMeanNumberOfMolecules = fNumberOfMolecules / fNbOfScoredEvents;
            fTime = iter.first;
            fMoleculeName = wIter.first;
            fNtuple->Fill();
        }
    }

    fNtuple->Write();
}

void TsScoreNumberOfMoleculesAtTime::Clear()
{
    fMoleculeCountPerIndexPerTime.clear();
    fNbOfScoredEvents = 0;
    UpdateFileNameForUpcomingRun();
}
