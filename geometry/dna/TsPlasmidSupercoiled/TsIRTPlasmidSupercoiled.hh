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

#ifndef TsIRTPlasmidSupercoiled_hh
#define TsIRTPlasmidSupercoiled_hh

#include "TsVGeometryComponent.hh"
#include "G4RotationMatrix.hh"
#include "G4ThreeVector.hh"
#include "G4Orb.hh"
#include "G4Point3D.hh"

#include <map>
#include <vector>

class G4VSolid;
class G4Material;
class G4LogicalVolume;
class TsGeneratorManager;

class TsIRTPlasmidSupercoiled : public TsVGeometryComponent
{
public:
	TsIRTPlasmidSupercoiled(TsParameterManager* pM, TsExtensionManager* eM, TsMaterialManager* mM, TsGeometryManager* gM,
						 TsVGeometryComponent* parentComponent, G4VPhysicalVolume* parentVolume, G4String& name);
	~TsIRTPlasmidSupercoiled();
	
	G4VPhysicalVolume* Construct();

	std::vector<G4String> GetDNANames() {return fDNANames;}
	std::vector<G4double> GetDNATimes() {return fDNATimes;}
	std::vector<G4ThreeVector> GetDNAPositions() {return fDNAPositions;}
	std::vector<std::vector<G4int>> GetDNADetails() {return fDNADetails;}

	void AddPlasmidInformation(G4int, G4RotationMatrix*, G4ThreeVector*);
	
public:
	G4LogicalVolume* CreateLogicVolumeXYZ(G4int, G4RotationMatrix*, G4ThreeVector*);
	G4LogicalVolume* CreateLogicVolumeDNAFabric(G4int,G4RotationMatrix*, G4ThreeVector*);
	
private:
	
	struct TempMolecule
	{
		TempMolecule(std::string name, int copyNumber, G4ThreeVector position, double radius, double waterRadius, std::string material, int strand)
		{
			fName = name;
			fMaterial = material;
			fCopyNumber = copyNumber;
			fPosition = position;
			fRadius = radius;
			fRadiusWater = waterRadius;
			fStrand = strand;
		}
		
		std::string fName;
		std::string fMaterial;
		
		int fCopyNumber;
		int fStrand;
		
		G4ThreeVector fPosition;
		
		double fRadius;
		double fRadiusWater;
		
		bool operator<(const TempMolecule& str) const
		{
			return (fPosition.z() < str.fPosition.z() );
		}
	};

	std::string fGeoName;
	std::map<std::string, G4double> fRadiusMap;
	std::map<std::string, G4double> fWaterRadiusMap;
	
	std::vector<TempMolecule> fMolecules;
	G4Orb* fSolid;

	std::vector<G4String> fDNANames;
	std::vector<G4double> fDNATimes;
	std::vector<G4ThreeVector> fDNAPositions;
	std::vector<std::vector<G4int>> fDNADetails;
	std::vector<G4ThreeVector> fVertexes;

	std::vector<G4String> fSampleDNANames;
	std::vector<G4double> fSampleDNATimes;
	std::vector<G4ThreeVector> fSampleDNAPositions;
	std::vector<std::vector<G4int>> fSampleDNADetails;
	
	void ReadXYZFile(G4String);
	void ReadDNAFabricFile(G4String);
	G4VSolid* CreateCutSolid(G4Orb *solidOrbRef,
							 TempMolecule &molRef,
							 std::vector<TempMolecule> &molList,
							 G4bool in);
	G4double fSize;
	G4double fOffsetX;
	G4double fOffsetY;
	G4double fOffsetZ;
	G4double fXMin;
	G4double fXMax;
	G4double fYMin;
	G4double fYMax;
	G4double fZMin;
	G4double fZMax;
 
	std::vector<G4int> fTracks;
};

#endif
