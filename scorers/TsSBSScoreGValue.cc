// Scorer for TsSBSGValue
//
// ********************************************************************
// *                                                                  *
// * This file is part of the TOPAS-nBio extensions to the            *
// *   TOPAS Simulation Toolkit.                                      *
// * The TOPAS-nBio extensions are freely available under the license *
// *   agreement set forth at: https://topas-nbio.readthedocs.io/     *
// *                                                                  *
// ********************************************************************

#include "TsSBSScoreGValue.hh"

#include "TsChemistryManager.hh"

#include "G4RunManager.hh"
#include "G4SystemOfUnits.hh"
#include "G4DNAChemistryManager.hh"
#include "G4Scheduler.hh"
#include "G4MoleculeCounter.hh"
#include "G4MolecularConfiguration.hh"
#include "G4Molecule.hh"
#include "G4H2O.hh"

TsSBSScoreGValue::TsSBSScoreGValue(TsParameterManager *pM, TsMaterialManager *mM, TsGeometryManager *gM, TsScoringManager *scM, TsExtensionManager *eM,
                                   G4String scorerName, G4String quantity, G4String outFileName, G4bool isSubScorer)
    : TsVNtupleScorer(pM, mM, gM, scM, eM, scorerName, quantity, outFileName, isSubScorer),
      fPm(pM), fEnergyDepositPerEvent(0)
{
    fNtuple->RegisterColumnD(&fGValue, "GValue: number of molecules per energy deposit per eV", "");
    fNtuple->RegisterColumnD(&fGValueError, "GValue statistical error", "");
    fNtuple->RegisterColumnD(&fTime, "Time", "picosecond");
    fNtuple->RegisterColumnS(&fMoleculeName, "MoleculeName");

    fNbTimeToRecord = fPm->GetVectorLength(GetFullParmName("TimeToRecord"));
    {
        auto times = std::unique_ptr<G4double[]>(fPm->GetDoubleVector(GetFullParmName("TimeToRecord"), "Time"));
        for (auto i = 0; i < fNbTimeToRecord; ++i)
            fTimesToRecord.push_back(times[i]);
    }

    fEnergyLossKill = fPm->GetDoubleParameter(GetFullParmName("KillPrimaryIfEnergyLossExceeds"), "Energy");
    fEnergyLossAbort = fPm->GetDoubleParameter(GetFullParmName("AbortEventIfPrimaryEnergyLossExceeds"), "Energy");
    fMaximumTrackLength = fPm->GetDoubleParameter(GetFullParmName("KillPrimaryBasedOnTrackLength"), "Length");

    fNbOfMoleculesToScavenge = -1;
    if (fPm->ParameterExists(GetFullParmName("ScavengeTheseMolecules")))
    {
        fNbOfMoleculesToScavenge = fPm->GetVectorLength(GetFullParmName("ScavengeTheseMolecules"));
        fMoleculeIDToScavenge = fPm->GetIntegerVector(GetFullParmName("ScavengeTheseMolecules"));
        fScavengingCapacity = fPm->GetDoubleVector(GetFullParmName("ScavengingCapacities"), "perTime");
    }

    fNbOfScoredEvents = 0;
    fEnergyLoss = 0.0;
    fTotalTrackLength = 0.0;

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

G4bool TsSBSScoreGValue::ProcessHits(G4Step *aStep, G4TouchableHistory *)
{
    if (!fIsActive)
    {
        fSkippedWhileInactive++;
        return false;
    }

    G4Track *aTrack = aStep->GetTrack();
    if (-1 < aTrack->GetTrackID())
    {
        if (aTrack->GetParentID() == 0)
        {
            G4double kineticE = aStep->GetPostStepPoint()->GetKineticEnergy();
            G4double eLoss = aStep->GetPreStepPoint()->GetKineticEnergy() - kineticE;
            fEnergyLoss += eLoss;
            fTotalTrackLength += aStep->GetStepLength();

            if (fEnergyLoss > fEnergyLossAbort)
            {
                std::cout << " -- Aborting event " << GetEventID() << std::endl;
                G4RunManager::GetRunManager()->AbortEvent();
                return true;
            }
            if (fEnergyLoss >= fEnergyLossKill)
            {
                std::cout << " -- Killing primary track of event " << GetEventID() << std::endl;
                aStep->GetTrack()->SetTrackStatus(fStopAndKill);
                return true;
            }
            if (fTotalTrackLength >= fMaximumTrackLength)
            {
                aStep->GetTrack()->SetTrackStatus(fStopAndKill);
                std::cout << " Track killed with track: " << fTotalTrackLength / nm << " nm with energy lost " << fEnergyLoss / keV << " keV " << std::endl;
                return true;
            }
        }

        G4double edep = aStep->GetTotalEnergyDeposit();
        if (edep > 0)
        {
            fEnergyDepositPerEvent += edep * aStep->GetPreStepPoint()->GetWeight();
            return true;
        }
    }
    else
    {
        if (0 < fNbOfMoleculesToScavenge)
        {
            G4Molecule *molecule = GetMolecule(aTrack);
            G4int moleculeID = molecule->GetMoleculeID();
            G4double t = aStep->GetPostStepPoint()->GetLocalTime() - aStep->GetPreStepPoint()->GetLocalTime();
            // this should only use the time diff since the last step (i.e. probability to have scavenged since then)
            for (int i = 0; i < fNbOfMoleculesToScavenge; i++)
            {
                if (moleculeID == fMoleculeIDToScavenge[i])
                {
                    G4double probability = 1. - std::exp(-fScavengingCapacity[i] * t);
                    if (G4UniformRand() < probability)
                    {
                        std::cout << " scavenged " << molecule->GetName() << std::endl;
                        aStep->GetTrack()->SetTrackStatus(fStopAndKill);
                        return true;
                    }
                }
            }
        }
    }

    return false;
}

void TsSBSScoreGValue::AccumulateEvent()
{
    if (!G4EventManager::GetEventManager()->GetConstCurrentEvent()->IsAborted() && fEnergyDepositPerEvent > 0)
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
                G4double gvalue = 100.0 * nbOfMoleculesAtTime / (fEnergyDepositPerEvent / eV);

                fGValuePerSpeciePerTime[key][time] += gvalue;
                fGValuePerSpeciePerTime2[key][time] += gvalue * gvalue;
            }
        }

        fNbOfScoredEvents++;
    }

    fTotalTrackLength = 0.0;
    fEnergyDepositPerEvent = 0.0;
    fEnergyLoss = 0.0;
}

void TsSBSScoreGValue::AbsorbResultsFromWorkerScorer(TsVScorer *workerScorer)
{
    TsVNtupleScorer::AbsorbResultsFromWorkerScorer(workerScorer);

    TsSBSScoreGValue *worker = dynamic_cast<TsSBSScoreGValue *>(workerScorer);

    fNbOfScoredEvents += worker->fNbOfScoredEvents;

    for (auto const &workerEntry : worker->fGValuePerSpeciePerTime)
    {
        auto emplacedPair = fGValuePerSpeciePerTime.emplace(workerEntry.first, workerEntry.second);
        if (!emplacedPair.second)
        { // workerEntry.first already exists, emplacedPair.first == it
            auto it = emplacedPair.first;

            for (auto const &workerGValueWithTimeIter : workerEntry.second)
            {
                auto masterGValueWithTimeIter = it->second.find(workerGValueWithTimeIter.first);
                masterGValueWithTimeIter->second += workerGValueWithTimeIter.second;
            }
        }
    }

    for (auto const &workerEntry : worker->fGValuePerSpeciePerTime2)
    {
        auto emplacedPair = fGValuePerSpeciePerTime2.emplace(workerEntry.first, workerEntry.second);
        if (!emplacedPair.second)
        { // workerEntry.first already exists, emplacedPair.first == it
            auto it = emplacedPair.first;

            for (auto const &workerGValueWithTimeIter : workerEntry.second)
            {
                auto masterGValueWithTimeIter = it->second.find(workerGValueWithTimeIter.first);
                masterGValueWithTimeIter->second += workerGValueWithTimeIter.second;
            }
        }
    }

    worker->fGValuePerSpeciePerTime.clear();
    worker->fGValuePerSpeciePerTime2.clear();
    worker->fNbOfScoredEvents = 0;
    worker->fEnergyDepositPerEvent = 0.0;
}

void TsSBSScoreGValue::Output()
{
    if (fNbOfScoredEvents == 0)
        return;

    std::map<G4String, std::map<G4double, G4double>>::iterator wIter;
    std::map<G4String, std::map<G4double, G4double>>::iterator wIter2;
    std::map<G4double, G4double>::iterator iter;
    std::map<G4double, G4double>::iterator iter2;

    for (wIter = fGValuePerSpeciePerTime.begin(); wIter != fGValuePerSpeciePerTime.end(); wIter++)
    {
        wIter2 = fGValuePerSpeciePerTime2.find(wIter->first);

        for (iter = (wIter->second).begin(); iter != (wIter->second).end(); iter++)
        {
            iter2 = (wIter2->second).find(iter->first);

            fGValue = iter->second / fNbOfScoredEvents;
            if (fNbOfScoredEvents > 1)
            {
                fGValueError = sqrt((1.0 / (fNbOfScoredEvents - 1)) * ((iter2->second) / fNbOfScoredEvents - fGValue * fGValue));
            }
            else
            {
                fGValueError = 1.0;
            }
            fTime = iter->first;
            fMoleculeName = wIter->first;
            fNtuple->Fill();
        }
    }

    fNtuple->Write();
}

void TsSBSScoreGValue::Clear()
{
    fGValuePerSpeciePerTime.clear();
    fNbOfScoredEvents = 0;
    UpdateFileNameForUpcomingRun();
}
