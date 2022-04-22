// Extra Class for TsScoreDNADamageSBS
//
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

#include "TsFociAnalysis.hh"

TsFociAnalysis::TsFociAnalysis()
{
	fFociSize = -1;
	fGet3DFociImage = false;
	fGet2DFociImage = false;
	fMicroscopePSFShape = "none";
	fMicroscopePSFWidth = -1;
}

TsFociAnalysis::~TsFociAnalysis() { }

G4int TsFociAnalysis::GetNumberOfFoci(std::vector<G4ThreeVector> dsbPositions)
{
	std::vector<G4bool> indexIsAvailable;
	for (G4int i = 0; i < dsbPositions.size(); i++)
		indexIsAvailable.push_back(true);

	std::vector<std::vector<G4int>> vectorOfDSBsInEachFocus;
	std::vector<G4int> dsbIdsInThisFocus;
	for (G4int i = 0; i < dsbPositions.size(); i++)
	{
		if (indexIsAvailable[i])
		{
			indexIsAvailable[i] = false;
			dsbIdsInThisFocus.push_back(i);
			for (G4int j = 0; j < dsbPositions.size(); j++)
			{
				if (indexIsAvailable[j] && GetDistance(dsbPositions[i], dsbPositions[j]) < fFociSize / 2)
				{
					indexIsAvailable[j] = false;
					dsbIdsInThisFocus.push_back(j);
				}
			}
			vectorOfDSBsInEachFocus.push_back(dsbIdsInThisFocus);
		}
		dsbIdsInThisFocus.clear();
	}
	for (G4int i = 0; i < vectorOfDSBsInEachFocus.size(); i++)
	{
		for (G4int j = 0; j < vectorOfDSBsInEachFocus[i].size(); j++)
			G4cout << vectorOfDSBsInEachFocus[i][j] << " ";
		G4cout << G4endl;
	}
	return vectorOfDSBsInEachFocus.size();
}

G4double TsFociAnalysis::GetDistance(G4ThreeVector a, G4ThreeVector b)
{
	return std::sqrt(std::pow(a.x()-b.x(), 2) + std::pow(a.y()-b.y(), 2) + std::pow(a.z()-b.z(), 2));
}
