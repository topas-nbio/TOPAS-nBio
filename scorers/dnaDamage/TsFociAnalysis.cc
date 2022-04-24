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
#include "G4VisExtent.hh"

TsFociAnalysis::TsFociAnalysis(TsVGeometryComponent* component)
{
	fComponent = component;

	fMicroscopePSFShape = "none";
	fMicroscopePSFWidth = -1;

	fxmin = -5 * um;
	fxmax = 5 * um;
	fymin = -5 * um;
	fymax = 5 * um;
	fzmin = -5 * um;
	fzmax = 5 * um;
}

TsFociAnalysis::~TsFociAnalysis() { }

std::vector<G4int> TsFociAnalysis::GetNumberOfFoci(std::vector<G4ThreeVector> dsbPositions)
{
	std::vector<G4int> numFoci;
	for (G4int iSize = 0; iSize < fFociSizes.size(); iSize++)
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
					if (indexIsAvailable[j] && GetDistance(dsbPositions[i], dsbPositions[j]) < fFociSizes[iSize] / 2)
					{
						indexIsAvailable[j] = false;
						dsbIdsInThisFocus.push_back(j);
					}
				}
				vectorOfDSBsInEachFocus.push_back(dsbIdsInThisFocus);
			}
			dsbIdsInThisFocus.clear();
		}
		numFoci.push_back(vectorOfDSBsInEachFocus.size());
	}
	return numFoci;
}

void TsFociAnalysis::Produce3DImage(std::vector<G4ThreeVector> dsbPositions)
{
	for (G4int iRes = 0; iRes < fResolutions.size(); iRes++)
	{
		// Creates 3D matrix
		G4float dx = (G4float)fResolutions[iRes];
		G4float dy = (G4float)fResolutions[iRes];
		G4float dz = (G4float)fResolutions[iRes];
		G4int nx = std::floor((fxmax-fxmin) / dx + 1e-14);
		G4int ny = std::floor((fymax-fymin) / dy + 1e-14);
		G4int nz = std::floor((fzmax-fzmin) / dz + 1e-14);
		G4float image3d[nx][ny][nz];

		// Initialization
		for (G4int i = 0; i < nx; i++)
		{
			for (G4int j = 0; j < ny; j++)
			{
				for (G4int k = 0; k < nz; k++)
					image3d[i][j][k] = 0.0;
			}
		}
		// Gets the PSF (Gaussian)
		if (strstr(fMicroscopePSFShape, "Gaussian") != NULL)
		{
			G4int twoSigmas = std::floor(2 * fMicroscopePSFWidth / (G4float)fResolutions[iRes]);
			G4int nkernel = 2 * twoSigmas + 1;
			G4float psf[nkernel][nkernel][nkernel];
			G4float minkernel = -twoSigmas * (G4float)fResolutions[iRes];
			G4float sum = 0;
			for (G4int ix = 0; ix < nkernel; ix++)
			{
				G4float xpos = minkernel + ix * (G4float)fResolutions[iRes];
				for (G4int iy = 0; iy < nkernel; iy++)
				{
					G4float ypos = minkernel + iy * (G4float)fResolutions[iRes];
					for (G4int iz = 0; iz < nkernel; iz++)
					{
						G4float zpos = minkernel + iz * (G4float)fResolutions[iRes];
						G4float v = Gaussian3D(xpos, ypos, zpos, fMicroscopePSFWidth);
						if (v > 1e-9)
						{
							psf[ix][iy][iz] = v;
							sum += v;
						}
						else psf[ix][iy][iz] = 0.0;
					}
				}
			}
			// Normalize
			for (G4int ix = 0; ix < nkernel; ix++)
			{
				for (G4int iy = 0; iy < nkernel; iy++)
				{
					for (G4int iz = 0; iz < nkernel; iz++)
						psf[ix][iy][iz] /= sum;
				}
			}
			// Convolves PSF for each DSB position
			for (G4int iDsb = 0; iDsb < dsbPositions.size(); iDsb++)
			{
				G4int idx = std::floor((dsbPositions[iDsb].x() - fxmin) / dx + 1e-14);
				G4int idy = std::floor((dsbPositions[iDsb].y() - fymin) / dy + 1e-14);
				G4int idz = std::floor((dsbPositions[iDsb].z() - fzmin) / dz + 1e-14);
				G4int startx = (idx - twoSigmas) < 0 ? 0 : idx - twoSigmas;
				G4int starty = (idy - twoSigmas) < 0 ? 0 : idy - twoSigmas;
				G4int startz = (idz - twoSigmas) < 0 ? 0 : idz - twoSigmas;
				G4int endx = (idx + twoSigmas) > nx ? nx : idx + twoSigmas;
				G4int endy = (idy + twoSigmas) > ny ? ny : idy + twoSigmas;
				G4int endz = (idz + twoSigmas) > nz ? nz : idz + twoSigmas;
				for (G4int i = startx; i < endx; i++)
				{
					for (G4int j = starty; j < endy; j++)
					{
						for (G4int k = startz; k < endz; k++)
							image3d[i][j][k] += psf[i - idx][j - idy][k - idz];
					}
				}
			}
		}

		// Writes csv file
		G4String filename = "Foci3D_" + std::to_string(int(fResolutions[iRes]*1e6)) + "nm.csv";
		std::fstream out;
		// Creates a new file
		out.open(filename, std::ios::out | std::ios::trunc);
		// Inserts headers
		out << "ix,iy,iz,v\n";
		for (G4int i = 0; i < nx; i++)
		{
			for (G4int j = 0; j < ny; j++)
			{
				for (G4int k = 0; k < nz; k++)
					for (G4int j = 0; j < nz; j++)
					{
						if (image3d[i][j][k] > 1e-9) out << i << "," << j << "," << k << "," << image3d[i][j][k] << "\n";
						else out << i << "," << j << "," << k << "," << 0.0 << "\n";
					}
			}
		}
		out.close();
	}
}

void TsFociAnalysis::Produce2DImages(std::vector<G4ThreeVector> dsbPositions)
{
	for (G4int iRes = 0; iRes < fResolutions.size(); iRes++)
	{
		for (G4int iPlane = 0; iPlane < f2DPlanesForFociImage.size(); iPlane++)
		{
			// Gets the PSF (Gaussian)
			G4int twoSigmas = std::floor(2 * fMicroscopePSFWidth / fResolutions[iRes]);
			G4int nkernel = 2 * twoSigmas + 1;
			G4float psf[nkernel][nkernel];

			if (strstr(fMicroscopePSFShape, "Gaussian") != NULL)
			{
				G4float minkernel = -twoSigmas * (G4float)fResolutions[iRes];
				G4float sum = 0;
				for (G4int ix = 0; ix < nkernel; ix++)
				{
					G4float xpos = minkernel + ix * (G4float)fResolutions[iRes];
					for (G4int iy = 0; iy < nkernel; iy++)
					{
						G4float ypos = minkernel + iy * (G4float)fResolutions[iRes];
						G4float v = Gaussian2D(xpos, ypos, fMicroscopePSFWidth);
						if (v >= 1e-9)
						{
							psf[ix][iy] = v;
							sum += v;
						}
						else
							psf[ix][iy] = 0.0;
					}
				}
				// Normalize
				for (G4int ix = 0; ix < nkernel; ix++)
				{
					for (G4int iy = 0; iy < nkernel; iy++)
						psf[ix][iy] /= sum;
				}
			}
			if (f2DPlanesForFociImage[iPlane] == "Z")
			{
				// Creates 2D matrix
				G4float dx = (G4float)fResolutions[iRes];
				G4float dy = (G4float)fResolutions[iRes];

				G4int nx = std::floor((fxmax-fxmin) / dx + 1e-14);
				G4int ny = std::floor((fymax-fymin) / dy + 1e-14);

				G4float image2d[nx][ny];
				// Initialization
				for (G4int i = 0; i < nx; i++)
				{
					for (G4int j = 0; j < ny; j++)
							image2d[i][j] = 0.0;
				}

				// Convolves PSF for each DSB position
				for (G4int iDsb = 0; iDsb < dsbPositions.size(); iDsb++)
				{
					G4int idx = std::floor((dsbPositions[iDsb].x() - fxmin) / dx + 1e-14);
					G4int idy = std::floor((dsbPositions[iDsb].y() - fymin) / dy + 1e-14);
					G4int startx = (idx - twoSigmas) < 0 ? 0 : idx - twoSigmas;
					G4int starty = (idy - twoSigmas) < 0 ? 0 : idy - twoSigmas;
					G4int endx = (idx + twoSigmas) > nx ? nx : idx + twoSigmas;
					G4int endy = (idy + twoSigmas) > ny ? ny : idy + twoSigmas;
					for (G4int i = startx; i < endx; i++)
					{
						for (G4int j = starty; j < endy; j++)
							image2d[i][j] += psf[i - idx][j - idy];
					}
				}
				// Writes csv file
				G4String filename = "Foci2D_ZPlane_" + std::to_string(int(fResolutions[iRes]*1e6)) + "nm.csv";
				std::fstream out;
				// Creates a new file
				out.open(filename, std::ios::out | std::ios::trunc);
				// Inserts headers
				out << "ix,iy,v\n";
				for (G4int i = 0; i < nx; i++)
				{
					for (G4int j = 0; j < ny; j++)
					{
						if (image2d[i][j] > 1e-9) out << i << "," << j << "," << image2d[i][j] << "\n";
						else out << i << "," << j << "," << 0.0 << "\n";
					}
				}
				out.close();
			}
			if (f2DPlanesForFociImage[iPlane] == "Y")
			{
				// Creates 2D matrix
				G4float dx = (G4float)fResolutions[iRes];
				G4float dz = (G4float)fResolutions[iRes];

				G4int nx = std::floor((fxmax-fxmin) / dx + 1e-14);
				G4int nz = std::floor((fzmax-fzmin) / dz + 1e-14);

				G4float image2d[nx][nz];
				// Initialization
				for (G4int i = 0; i < nx; i++)
				{
					for (G4int j = 0; j < nz; j++)
						image2d[i][j] = 0;
				}

				// Convolves PSF for each DSB position
				for (G4int iDsb = 0; iDsb < dsbPositions.size(); iDsb++)
				{
					G4int idx = std::floor((dsbPositions[iDsb].x() - fxmin) / dx + 1e-14);
					G4int idz = std::floor((dsbPositions[iDsb].z() - fzmin) / dz + 1e-14);
					G4int startx = (idx - twoSigmas) < 0 ? 0 : idx - twoSigmas;
					G4int startz = (idz - twoSigmas) < 0 ? 0 : idz - twoSigmas;
					G4int endx = (idx + twoSigmas) > nx ? nx : idx + twoSigmas;
					G4int endz = (idz + twoSigmas) > nz ? nz : idz + twoSigmas;
					for (G4int i = startx; i < endx; i++)
					{
						for (G4int j = startz; j < endz; j++)
							image2d[i][j] += psf[i - idx][j - idz];
					}
				}
				// Writes csv file
				G4String filename = "Foci2D_YPlane_" + std::to_string(int(fResolutions[iRes]*1e6)) + "nm.csv";
				std::fstream out;
				// Creates a new file
				out.open(filename, std::ios::out | std::ios::trunc);
				// Inserts headers
				out << "ix,iz,v\n";
				for (G4int i = 0; i < nx; i++)
				{
					for (G4int j = 0; j < nz; j++)
					{
						if (image2d[i][j] > 1e-9) out << i << "," << j << "," << image2d[i][j] << "\n";
						else out << i << "," << j << "," << 0.0 << "\n";
					}
				}
				out.close();
			}
			if (f2DPlanesForFociImage[iPlane] == "X")
			{
				// Creates 2D matrix
				G4float dy = (G4float)fResolutions[iRes];
				G4float dz = (G4float)fResolutions[iRes];

				G4int ny = std::floor((fymax-fymin) / dy + 1e-14);
				G4int nz = std::floor((fzmax-fzmin) / dz + 1e-14);

				G4float image2d[ny][nz];
				// Initialization
				for (G4int i = 0; i < ny; i++)
				{
					for (G4int j = 0; j < nz; j++)
						image2d[i][j] = 0;
				}

				// Convolves PSF for each DSB position
				for (G4int iDsb = 0; iDsb < dsbPositions.size(); iDsb++)
				{
					G4int idy = std::floor((dsbPositions[iDsb].y() - fymin) / dy + 1e-14);
					G4int idz = std::floor((dsbPositions[iDsb].z() - fzmin) / dz + 1e-14);
					G4int starty = (idy - twoSigmas) < 0 ? 0 : idy - twoSigmas;
					G4int startz = (idz - twoSigmas) < 0 ? 0 : idz - twoSigmas;
					G4int endy = (idy + twoSigmas) > ny ? ny : idy + twoSigmas;
					G4int endz = (idz + twoSigmas) > nz ? nz : idz + twoSigmas;
					for (G4int i = starty; i < endy; i++)
					{
						for (G4int j = startz; j < endz; j++)
							image2d[i][j] += psf[i - idy][j - idz];
					}
				}
				// Writes csv file
				G4String filename = "Foci2D_XPlane_" + std::to_string(int(fResolutions[iRes]*1e6)) + "nm.csv";
				std::fstream out;
				// Creates a new file
				out.open(filename, std::ios::out | std::ios::trunc);
				// Inserts headers
				out << "iy,iz,v\n";
				for (G4int i = 0; i < ny; i++)
				{
					for (G4int j = 0; j < nz; j++)
					{
						if (image2d[i][j] > 1e-9) out << i << "," << j << "," << image2d[i][j] << "\n";
						else out << i << "," << j << "," << 0.0 << "\n";
					}
				}
				out.close();
			}
		}
	}
}

G4float TsFociAnalysis::GetDistance(G4ThreeVector a, G4ThreeVector b)
{
	return std::sqrt(std::pow(a.x()-b.x(), 2) + std::pow(a.y()-b.y(), 2) + std::pow(a.z()-b.z(), 2));
}

G4float TsFociAnalysis::Gaussian3D(G4float x, G4float y, G4float z, G4float sigma)
{
	return std::exp(-(x*x+y*y+z*z)/(2*sigma*sigma));
}

G4float TsFociAnalysis::Gaussian2D(G4float x, G4float y, G4float sigma)
{
	return std::exp(-(x*x+y*y)/(2*sigma*sigma));
}
