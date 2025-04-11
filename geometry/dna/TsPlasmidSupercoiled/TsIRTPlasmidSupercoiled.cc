// Component for TsIRTPlasmidSupercoiled
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

#include "TsIRTPlasmidSupercoiled.hh"
#include "TsParameterManager.hh"

#include "G4VPhysicalVolume.hh"
#include "G4LogicalVolume.hh"
#include "G4PVPlacement.hh"
#include "G4RotationMatrix.hh"

#include "G4SubtractionSolid.hh"
#include "G4UnionSolid.hh"
#include "G4Box.hh"
#include "G4Orb.hh"
#include "G4Point3D.hh"
#include "G4Ellipsoid.hh"
#include "G4Scene.hh"

#include "G4SystemOfUnits.hh"
#include "G4PhysicalConstants.hh"
#include "G4UnionSolid.hh"
#include "G4VisAttributes.hh"
#include "G4Colour.hh"
#include "Randomize.hh"
#include "G4UIcommand.hh"

#include "G4NistManager.hh"
#include <sstream>

//TsIRTPlasmidSupercoiled* TsIRTPlasmidSupercoiled::theInstance = 0;

TsIRTPlasmidSupercoiled::TsIRTPlasmidSupercoiled(TsParameterManager* pM, TsExtensionManager* eM, TsMaterialManager* mM, TsGeometryManager* gM,
										   TsVGeometryComponent* parentComponent, G4VPhysicalVolume* parentVolume, G4String& name) :
TsVGeometryComponent(pM, eM, mM, gM, parentComponent, parentVolume, name) {;}


TsIRTPlasmidSupercoiled::~TsIRTPlasmidSupercoiled(){;}


G4VPhysicalVolume* TsIRTPlasmidSupercoiled::Construct()
{
	BeginConstruction();
	G4int nPlasdmids = fPm->GetIntegerParameter(GetFullParmName("NumberOfPlasmids"));
	G4double radius = fPm->GetDoubleParameter(GetFullParmName("R"),"Length");
	G4String fileName = fPm->GetStringParameter(GetFullParmName("InputFile"));
	
	fSolid = new G4Orb(fName, radius);
	fEnvelopeLog = CreateLogicalVolume(fSolid);
	fEnvelopePhys = CreatePhysicalVolume(fEnvelopeLog);
	
	fOffsetX = 0 * nm;
	fOffsetY = 0 * nm;
	fOffsetZ = 0 * nm;

	G4bool isDNAFabric = true;
	if   ( G4StrUtil::contains(fileName,".fab2g4dna")) 
		ReadDNAFabricFile(fileName);
	else { 
		ReadXYZFile(fileName);
		isDNAFabric = false; 
	}
	
	std::vector<G4RotationMatrix*> rotations;
	std::vector<G4ThreeVector*> translations;
	G4String theName = fPm->GetStringParameter(GetFullParmName("EnvelopeFileName"));
 
	std::ifstream envFile(theName);

    if (!envFile.is_open()) {
        G4String msg = theName+" could not be opened";
        G4Exception("TsIRTPlasmidSupercoiled::Construct()", "ENVELOPE_NotOpened", FatalException, msg);
    }

	G4double x, y, z, u, v, w;
	G4int n = 0, volID=0;
	while( n < nPlasdmids ) {
		envFile >> x >> y >> z >> u >> v >> w;
		if ( !envFile.good() )
			break;
		x *= nm;
		y *= nm;
		z *= nm;
		u *= deg;
		v *= deg;
		w *= deg;
		G4RotationMatrix* rot = new G4RotationMatrix();
		rot->rotateX(u);
		rot->rotateY(v);
		rot->rotateZ(w);
		G4ThreeVector* pos = new G4ThreeVector(x, y, z);
		
		rotations.push_back(rot);
		translations.push_back(pos);
		
		volID++;
		n++;
	}
	envFile.close();
	
	G4LogicalVolume* subComponent = 0;
	G4String plasmidName = "plasmid";
	for ( int p = 0; p < nPlasdmids; p++ ) {
		if ( p == 0 )
			if (isDNAFabric)
				subComponent = CreateLogicVolumeDNAFabric(p,rotations[p],translations[p]);
			else
				subComponent = CreateLogicVolumeXYZ(p,rotations[p],translations[p]);
		else
			CreatePhysicalVolume(plasmidName, p, true, subComponent, rotations[p], translations[p], fEnvelopePhys);
		AddPlasmidInformation(p, rotations[p], translations[p]);
	}

	InstantiateChildren(fEnvelopePhys);

	//theInstance = this;
	
	return fEnvelopePhys;
}


G4LogicalVolume* TsIRTPlasmidSupercoiled::CreateLogicVolumeXYZ(G4int copy, G4RotationMatrix* rot, G4ThreeVector* trans) {
	G4double des1 = 1.1344640137963142;  
	G4double des2 = des1 + (pi*.5);
	G4double ang = 0.6283185307179586;
	G4double bet1 = 0.6283185307179586 * 2;
	G4double posi = 1.0471975511965976;
	G4double sep = .1*angstrom;

	//Geometries Sizes
	G4double DeoxyriboseSize = 2.9389169420478556 * angstrom; 
	G4double PhosphateSize   = 2.7 * angstrom;
	G4double BaseSize        = 2.45 * angstrom;
 
	G4double xin = -170 * angstrom;
	G4double yin = -170 * angstrom;
	G4double zin = -170 * angstrom;
	G4double xfn =  170 * angstrom;
	G4double yfn =  170 * angstrom;
	G4double zfn =  170 * angstrom; 

	G4int nVertex = fVertexes.size();

	/////////////////////////////////////////////////////////////////////////////////////////////////
	//#####                                  Envelope                                          ####//
	/////////////////////////////////////////////////////////////////////////////////////////////////
	//Box

	std::string boxNameSolid = fGeoName + "_solid";
	G4Box* box_solid = new G4Box(boxNameSolid, 0.5*(fXMax-fXMin)+0.5*3.4*nm,
								               0.5*(fYMax-fYMin)+0.5*3.4*nm,
								               0.5*(fZMax-fZMin)+0.5*3.4*nm);

	G4String boxNameLogic       = fGeoName + "_logic";
	G4LogicalVolume* box_logic  = CreateLogicalVolume(boxNameLogic, box_solid);
	G4VPhysicalVolume* box_phys = CreatePhysicalVolume(boxNameLogic, copy, true, box_logic, rot, trans, fEnvelopePhys);

	/////////////////////////////////////////////////////////////////////////////////////////////////
	//#####                                    DNA                                             ####//
	/////////////////////////////////////////////////////////////////////////////////////////////////
	//Desoxyribose

	G4Ellipsoid* Deoxyribose = new G4Ellipsoid("deoxyribose1", 1.05*DeoxyriboseSize, 1.05*DeoxyriboseSize, 1.05*DeoxyriboseSize, -DeoxyriboseSize, .445*DeoxyriboseSize);
	G4LogicalVolume * LogicDeoxyribose = CreateLogicalVolume("deoxyribose1", Deoxyribose);

	/////////////////////////////////////////////////////////////////////////////////////////////////
	//Phosphoric Acid

	G4Ellipsoid* Phosphate = new G4Ellipsoid("phosphate1", 1.05*PhosphateSize, 1.05*PhosphateSize, 1.05*PhosphateSize, -PhosphateSize, .9 * angstrom);
	G4LogicalVolume * LogicPhosphate = CreateLogicalVolume("phosphate1", Phosphate);

	/////////////////////////////////////////////////////////////////////////////////////////////////
	//Base Pairs

	G4Ellipsoid* Base1 = new G4Ellipsoid("base_adenine", BaseSize, BaseSize, BaseSize, -BaseSize, 1.15 * angstrom);
	G4LogicalVolume* LogicBase1 = CreateLogicalVolume("base_adenine", Base1);

	G4Ellipsoid* Base2 = new G4Ellipsoid("base_thymine", BaseSize, BaseSize, BaseSize, -BaseSize, 1.15 * angstrom);
	G4LogicalVolume* LogicBase2 = CreateLogicalVolume("base_thymine", Base2);

	///////////////////////////////////////////////////////////////////////////////////////////////

	G4int index = 0; 
	G4double cAngle = 0;
	for (int vertex = 0; vertex < nVertex - 1; vertex++ ) {
		xin = fVertexes[vertex][0]-fOffsetX;
		yin = fVertexes[vertex][1]-fOffsetY;
		zin = fVertexes[vertex][2]-fOffsetZ;
		xfn = fVertexes[vertex+1][0]-fOffsetX;
		yfn = fVertexes[vertex+1][1]-fOffsetY;
		zfn = fVertexes[vertex+1][2]-fOffsetZ;

		//Calculo de los angulos para los segmentos
		G4double phi0   = std::atan2(zfn - zin,std::sqrt( ((xfn-xin) * (xfn-xin)) + ((yfn-yin) * (yfn-yin)) ));
		G4double theta0 = std::atan2(xfn-xin,yfn-yin);

		//Calculo del numero de segmentos 
		G4double lenght = std::sqrt(((xfn-xin)*(xfn-xin)) + ((yfn-yin)*(yfn-yin)) + ((zfn-zin)*(zfn-zin)));
		G4double dl = 1.0 / (lenght / (3.4*angstrom));

		G4int nChain = (fVertexes[vertex] - fVertexes[vertex+1]).mag() / (0.34 * nm);

		if (nChain == 0) {nChain++;}

		for (G4int nseg = 0; nseg < nChain ; nseg++) {
		  cAngle += ang;

		  //Posiciones originales de los segmentos
		  G4double theta = cAngle;
		  G4double x1 = 0;
		  G4double y1 = 0;
		  G4double z1 = ((2*BaseSize) + DeoxyriboseSize + sep);

		  G4double x2 = 0; 
		  G4double y2 = 0; 
		  G4double z2 = ((2*BaseSize) + DeoxyriboseSize + sep);

		  G4ThreeVector plus2 = G4ThreeVector(0,0,(.5 * DeoxyriboseSize) + PhosphateSize);
		  plus2.rotateX(-des1);
		  plus2.rotateZ(-posi);

		  G4ThreeVector plus2alt = G4ThreeVector(0,0,(.5 * DeoxyriboseSize) + PhosphateSize);
		  plus2alt.rotateX(-des1);
		  plus2alt.rotateZ(posi);

		  G4double x3 = 0   ;
		  G4double y3 = 0   ;
		  G4double z3 = BaseSize + sep;

		  G4ThreeVector position1i	= G4ThreeVector(x1,y1,z1);
		  G4ThreeVector position2i	= G4ThreeVector(x2,y2,z2) + plus2;
		  G4ThreeVector position2ialt = G4ThreeVector(x2,y2,-z2) - plus2alt;
		  G4ThreeVector position3i	= G4ThreeVector(x3,y3,z3);

		  //Rotacion de las posiciones originales
		  position1i.rotateY(theta);
		  position2i.rotateY(theta);
		  position2ialt.rotateY(theta);
		  position3i.rotateY(theta);

		  //Calculo de la posicion respect a la recta
		  G4double x = dl*nseg*(xfn-xin) + xin;
		  G4double y = dl*nseg*(yfn-yin) + yin;
		  G4double z = dl*nseg*(zfn-zin) + zin;

		  //Posiciones rotadas
			//Rotacion en el eje x
		  position1i.rotateX(phi0);
		  position2i.rotateX(phi0);
		  position2ialt.rotateX(phi0);
		  position3i.rotateX(phi0);

		  //Rotacion en el eje z
		  position1i.rotateZ(-theta0);
		  position2i.rotateZ(-theta0);
		  position2ialt.rotateZ(-theta0);
		  position3i.rotateZ(-theta0);

		  G4double yrot1 = theta, xrot1 = -des1;
		  G4RotationMatrix rotm1 = G4RotationMatrix();
		  rotm1.rotateX(xrot1);
		  rotm1.rotateZ(-posi);
		  rotm1.rotateY(yrot1);
		  rotm1.rotateX(phi0);
		  rotm1.rotateZ(-theta0);
		  G4ThreeVector position1 = position1i + G4ThreeVector(x,y,z);
	
		  G4double yrot1alt = theta + pi, xrot1alt = des1;
		  G4RotationMatrix rotm1alt = G4RotationMatrix();
		  rotm1alt.rotateX(xrot1alt);
		  rotm1alt.rotateZ(-posi);
		  rotm1alt.rotateY(yrot1alt);
		  rotm1alt.rotateX(phi0);
		  rotm1alt.rotateZ(-theta0);
		  G4ThreeVector position1alt = -position1i + G4ThreeVector(x,y,z);

		  G4double yrot2 = theta, xrot2 = -des2;
		  G4RotationMatrix rotm2 = G4RotationMatrix();
		  rotm2.rotateX(xrot2);
		  rotm2.rotateY(yrot2 - bet1 + 0.8726646259971648);
		  rotm2.rotateX(phi0);
		  rotm2.rotateZ(-theta0);
		  G4ThreeVector position2 = position2i + G4ThreeVector(x,y,z);

		  G4double yrot2alt = theta + pi, xrot2alt = des2;
		  G4RotationMatrix rotm2alt = G4RotationMatrix();
		  rotm2alt.rotateX(xrot2alt);
		  rotm2alt.rotateY(yrot2alt + bet1 - 0.8726646259971648);
		  rotm2alt.rotateX(phi0);
		  rotm2alt.rotateZ(-theta0);
		  G4ThreeVector position2alt = position2ialt + G4ThreeVector(x,y,z);

		  G4double yrot3 = theta;
		  G4RotationMatrix rotm3 = G4RotationMatrix();
		  rotm3.rotateX(-pi/2);
		  rotm3.rotateZ(-ang);
		  rotm3.rotateY(yrot3);
		  rotm3.rotateX(phi0);
		  rotm3.rotateZ(-theta0);
		  G4ThreeVector position3 = position3i + G4ThreeVector(x,y,z);

		  G4double yrot3alt = theta + pi;
		  G4RotationMatrix rotm3alt = G4RotationMatrix();
		  rotm3alt.rotateX(pi/2);
		  rotm3alt.rotateZ(-ang);
		  rotm3alt.rotateY(yrot3alt);
		  rotm3alt.rotateX(phi0);
		  rotm3alt.rotateZ(-theta0);
		  G4ThreeVector position3alt = -position3i + G4ThreeVector(x,y,z);

		 
		  CreatePhysicalVolume("deoxyribose1", index, true, LogicDeoxyribose,
		                        new G4RotationMatrix(rotm1), new G4ThreeVector(position1), box_phys);

		  CreatePhysicalVolume("phosphate1", index, true, LogicPhosphate,
		                        new G4RotationMatrix(rotm2), new G4ThreeVector(position2), box_phys);

		  CreatePhysicalVolume("base_adenine", index, true, LogicBase1, 
		  	                    new G4RotationMatrix(rotm3), new G4ThreeVector(position3), box_phys);

		  CreatePhysicalVolume("base_thymine", index, true, LogicBase2, 
		  	                    new G4RotationMatrix(rotm3alt), new G4ThreeVector(position3alt), box_phys);

		  CreatePhysicalVolume("phosphate2", index, true, LogicPhosphate, 
		  	                    new G4RotationMatrix(rotm2alt), new G4ThreeVector(position2alt), box_phys);

		  CreatePhysicalVolume("deoxyribose2", index, true, LogicDeoxyribose, 
		  	                    new G4RotationMatrix(rotm1alt), new G4ThreeVector(position1alt), box_phys);

		  G4ThreeVector Deoxy1 = (position1    + position2)    / 2;
		  G4ThreeVector Deoxy2 = (position1alt + position2alt) / 2;

		  fSampleDNANames.push_back("deoxyribose1");
		  fSampleDNATimes.push_back(1*ps);
		  fSampleDNAPositions.push_back(Deoxy1);
		  fSampleDNADetails.push_back({-1,index,1});

		  fSampleDNANames.push_back("deoxyribose2");
		  fSampleDNATimes.push_back(1*ps);
		  fSampleDNAPositions.push_back(Deoxy2);
		  fSampleDNADetails.push_back({-1,index,2});

		  fSampleDNANames.push_back("base1");
		  fSampleDNATimes.push_back(1*ps);
		  fSampleDNAPositions.push_back(position3);
		  fSampleDNADetails.push_back({-1,index,1});

		  fSampleDNANames.push_back("base2");
		  fSampleDNATimes.push_back(1*ps);
		  fSampleDNAPositions.push_back(position3alt);
		  fSampleDNADetails.push_back({-1,index,2});

		  index++;
		
		}
	}
	return box_logic;	
}


G4LogicalVolume* TsIRTPlasmidSupercoiled::CreateLogicVolumeDNAFabric(G4int copy, G4RotationMatrix* rot, G4ThreeVector* trans)
{
	G4String materialBox = fParentComponent->GetResolvedMaterialName();
	G4String materialDNA = GetResolvedMaterialName();
	
	std::string boxNameSolid = fGeoName+"_solid";
	G4Box* box_solid = new G4Box(boxNameSolid, 0.5*(fXMax-fXMin)+0.5*3.4*nm,
								 0.5*(fYMax-fYMin)+0.5*3.4*nm,
								 0.5*(fZMax-fZMin)+0.5*3.4*nm);
	G4String boxNameLogic = fGeoName+"_logic";
	G4LogicalVolume* box_logic = CreateLogicalVolume(boxNameLogic, box_solid);
	
	G4VPhysicalVolume* box_phys = CreatePhysicalVolume(boxNameLogic, copy, true, box_logic, rot, trans, fEnvelopePhys);
	
	G4ThreeVector offsetPosition(fOffsetX, fOffsetY, fOffsetZ);
	
	G4VisAttributes* purines = fPm->GetColor("red");
	purines->SetForceSolid();
	G4VisAttributes* pyrimidines = fPm->GetColor("yellow");
	pyrimidines->SetForceSolid();
	G4VisAttributes* deoxyriboses = fPm->GetColor("magenta");
	deoxyriboses->SetForceSolid();
	G4VisAttributes* phosphates = fPm->GetColor("grass");
	phosphates->SetForceSolid();
	G4VisAttributes* histones = fPm->GetColor("blue");
	histones->SetForceSolid();
	G4VisAttributes* wshell = new G4VisAttributes(G4Color(0,0,1,0.1));
	wshell->SetForceSolid();
	RegisterVisAtt(wshell);
	
	G4Orb* purineAndPyrimidine = 0;
	G4Orb* phosphate = 0;
	G4Orb* deoxyribose = 0;
	G4Orb* wpurineAndPyrimidine = 0;
	G4Orb* wphosphate = 0;
	G4Orb* wdeoxyribose = 0;
	G4Orb* histone = 0;
	G4bool solidPurinePyrimidine = true;
	G4bool solidPhosphate = true;
	G4bool solidDeoxyribose = true;
	G4bool solidHistone = true;
	
	G4cout << "#### Creating plasmid base pairs: begin" << G4endl;
	for(size_t i = 0; i < fMolecules.size(); i++)
	{
		G4String name        = fMolecules[i].fName;
		G4double radius      = fMolecules[i].fRadius;
		G4double waterRadius = fMolecules[i].fRadiusWater;
		G4ThreeVector moleculePosition = fMolecules[i].fPosition - offsetPosition;
		G4int copyNum = fMolecules[i].fCopyNumber;

		fSampleDNANames.push_back(name);
		fSampleDNATimes.push_back(1*ps);
		fSampleDNAPositions.push_back(moleculePosition);
		fSampleDNADetails.push_back({-1,copyNum,fMolecules[i].fStrand});
		
		if ( solidPhosphate || solidPurinePyrimidine || solidDeoxyribose || solidHistone ) {
			if ( G4StrUtil::contains(name,"phosphate") && solidPhosphate) {
				phosphate = new G4Orb("phosphate", radius);
				wphosphate = new G4Orb("wphosphate", waterRadius);
				solidPhosphate = false;
				G4cout << "####  Built phosphate" << G4endl;
			} else if ( G4StrUtil::contains(name,"base") && solidPurinePyrimidine ) {
				purineAndPyrimidine = new G4Orb("purineandpyrimidine", radius);
				wpurineAndPyrimidine = new G4Orb("wpurineandpyrimidine", waterRadius);
				solidPurinePyrimidine = false;
				G4cout << "####  Built purines and pyrimidines" << G4endl;
			} else if ( G4StrUtil::contains(name,"deoxyribose") && solidDeoxyribose ) {
				deoxyribose = new G4Orb("deoxyribose", radius);
				wdeoxyribose = new G4Orb("wdeoxyribose", waterRadius);
				solidDeoxyribose = false;
				G4cout << "####  Built deoxyriboses" << G4endl;
			} else if ( G4StrUtil::contains(name,"histone") && solidHistone ) {
				histone = new G4Orb("histone", radius);
				solidHistone = false;
			}
		}
		// Water hydration shell volume part
		G4VPhysicalVolume* moleculeWater_phys = 0;
		// If water radius != 0 then we have a water hydration shell
		G4double tol = 0.0001;
		G4VSolid* moleculeWaterCut_solid = 0;
		G4LogicalVolume* moleculeWater_logic = 0;
		if(waterRadius > (0 + tol)*nm)
		{
			G4String nameWaterSolid = name+"_waterShell";
			if ( G4StrUtil::contains(name,"phosphate") ) {
				moleculeWaterCut_solid = CreateCutSolid(wphosphate, fMolecules[i], fMolecules, false);
			} else if ( G4StrUtil::contains(name,"base") ) {
				moleculeWaterCut_solid = CreateCutSolid(wpurineAndPyrimidine, fMolecules[i], fMolecules, false);
			} else if ( G4StrUtil::contains(name,"deoxyribose") ) {
				moleculeWaterCut_solid = CreateCutSolid(wdeoxyribose, fMolecules[i], fMolecules, false);
			}
			
			//moleculeWater_logic = CreateLogicalVolume(nameWaterSolid, materialDNA, moleculeWaterCut_solid);
			moleculeWater_logic = CreateLogicalVolume(nameWaterSolid, moleculeWaterCut_solid);
			//moleculeWater_logic->SetVisAttributes(wshell);
			moleculeWater_phys = CreatePhysicalVolume(nameWaterSolid, copyNum, true, moleculeWater_logic, 0, new G4ThreeVector(moleculePosition), box_phys);
		}
		
		// Dna volume part
		G4VSolid* moleculeCut_solid = 0;
		G4String nameSolid = fMolecules[i].fName+"_solid";
		
		G4String nameLogic = name;
		G4LogicalVolume* molecule_logic = 0;
		
		if ( G4StrUtil::contains(name,"phosphate") ) {
			moleculeCut_solid = CreateCutSolid(phosphate, fMolecules[i], fMolecules, true);
		} else if ( G4StrUtil::contains(name,"base") ) {
			moleculeCut_solid = CreateCutSolid(purineAndPyrimidine, fMolecules[i], fMolecules, true);
		} else if ( G4StrUtil::contains(name,"deoxyribose") ) {
			moleculeCut_solid = CreateCutSolid(deoxyribose, fMolecules[i], fMolecules, true);
		} else if ( G4StrUtil::contains(name,"histone") ) {
			moleculeCut_solid = CreateCutSolid(histone, fMolecules[i], fMolecules, true);
		}
		
		molecule_logic = CreateLogicalVolume(nameLogic, moleculeCut_solid);
		
		G4ThreeVector position(0);
		G4String namePhys = name;
		if(waterRadius > (0 + tol)*nm)
			CreatePhysicalVolume(namePhys, copyNum, true, molecule_logic, 0, new G4ThreeVector(position), moleculeWater_phys);
		else
			CreatePhysicalVolume(namePhys, copyNum, true, molecule_logic, 0, new G4ThreeVector(moleculePosition), box_phys);
	}
	G4cout << "#### Creating plasmid base pairs: end" << G4endl;
	// Clear the containers

	fMolecules.clear();
	fRadiusMap.clear();
	fWaterRadiusMap.clear();
	return box_logic;
}

void TsIRTPlasmidSupercoiled::AddPlasmidInformation(G4int copy, G4RotationMatrix* rot, G4ThreeVector* trans) {
	G4Point3D aPoint         = G4Point3D(trans->x(), trans->y(), trans->z());
	G4RotationMatrix aInvRot = G4RotationMatrix(rot->inverse());

	for (size_t i = 0; i < fSampleDNANames.size(); i++) {
		if ( G4StrUtil::contains(fSampleDNANames[i],"deoxyribose") ) {
			fDNATimes.push_back(fSampleDNATimes[i]);

			std::vector<G4int> Details = fSampleDNADetails[i];
			Details[0] = copy;
			fDNADetails.push_back(Details);

			fDNANames.push_back("deox");

			G4double x = fSampleDNAPositions[i].x();
			G4double y = fSampleDNAPositions[i].y();
			G4double z = fSampleDNAPositions[i].z();

			G4Point3D newPoint = G4Translate3D(aInvRot * G4Point3D(x,y,z)) * (aPoint);

			fDNAPositions.push_back(G4ThreeVector(newPoint.x(), newPoint.y(), newPoint.z()));
		}
	}
}

void TsIRTPlasmidSupercoiled::ReadXYZFile(G4String fileName)
{
	G4double x, y, z;
	fXMin = 1*mm, fYMin = 1*mm, fZMin = 1*mm;
	fXMax=0.0, fYMax=0.0, fZMax=0.0;
	std::ifstream plasmidFile(fileName); //"envelopes.xyz");
	while(true) {
		plasmidFile >> x >> y >> z;
		if ( !plasmidFile.good() ) break;
		x *= nm;
		y *= nm;
		z *= nm;
		fVertexes.push_back(G4ThreeVector(x, y, z));
		if ( fXMin > x ) fXMin = x;
		if ( fXMax < x ) fXMax = x;
		if ( fYMin > y ) fYMin = y;
		if ( fYMax < y ) fYMax = y;
		if ( fZMin > z ) fZMin = z;
		if ( fZMax < z ) fZMax = z;
	}
	plasmidFile.close();
	fOffsetX = (fXMin + fXMax)*0.5;
	fOffsetY = (fYMin + fYMax)*0.5;
	fOffsetZ = (fZMin + fZMax)*0.5;
	fGeoName = "VoxelStraight";
}

void TsIRTPlasmidSupercoiled::ReadDNAFabricFile(G4String fileName)
{
	// Clear the containers
	fMolecules.clear();
	fRadiusMap.clear();
	fWaterRadiusMap.clear();

	fSampleDNAPositions.clear();
	fSampleDNATimes.clear();
	fSampleDNANames.clear();
	fSampleDNADetails.clear();

	fDNAPositions.clear();
	fDNATimes.clear();
	fDNANames.clear();
	fDNADetails.clear();
	
	//G4bool addMoleFlag = true;
	
	// Setup the input stream
	std::ifstream file(fileName.c_str());
	
	// Check if the file was correctly opened
	if(!file.is_open())
	{
		G4String msg = fileName+" could not be opened";
		G4Exception("PhysGeoImport::ReadDNAFabricFile()", "Geo_InputFileNotOpened", FatalException, msg);
	}
	
	fXMin = 1*mm, fYMin = 1*mm, fZMin = 1*mm;
	fXMax=0.0, fYMax=0.0, fZMax=0.0;
	
	// Define the line string variable
	std::string line;
	// Read the file line per line
	while(std::getline(file, line) )
	{
		// Check the line to determine if it is empty
		if(line.empty()) continue; // skip the line if it is empty
		
		// Data string stream
		std::istringstream issLine(line);
		
		// String to determine the first letter/word
		std::string firstItem;
		
		// Put the first letter/word within the string
		issLine >> firstItem;
		
		// Check first letter to determine if the line is data or comment
		if(firstItem=="#") continue; // skip the line if it is comment
		
		// Use the file
		else if(firstItem=="_Name")
		{
			std::string name;
			issLine >> name;
			
			fGeoName = name;
		}
		else if(firstItem=="_Size")
		{
			G4double size;
			issLine >> size;
			size *= nm;
			
			fSize = size;
		}
		else if(firstItem == "_Version") {;}
		else if(firstItem=="_Number") {;}
		else if(firstItem=="_Radius")
		{
			std::string name;
			issLine >> name;
			
			G4double radius;
			issLine >> radius;
			radius *= nm;
			
			G4double waterRadius;
			issLine >> waterRadius;
			waterRadius *= nm;
			
			fRadiusMap[name] = radius;
			fWaterRadiusMap[name] = waterRadius;
		}
		else if(firstItem=="_pl")
		{
			std::string name;
			issLine >> name;
			
			std::string material;
			issLine >> material;
			
			G4int strand;
			issLine >> strand;
			
			G4int copyNumber;
			issLine >> copyNumber;
			
			G4double x;
			issLine >> x;
			x *= nm;
			
			G4double y;
			issLine >> y;
			y *= nm;
			
			G4double z;
			issLine >> z;
			z *= nm;
			
			TempMolecule molecule(name, copyNumber, G4ThreeVector(x, y, z), fRadiusMap[name], fWaterRadiusMap[name], material, strand);
			fMolecules.push_back(molecule);
			if ( fXMin > x )
				fXMin = x;
			if ( fXMax < x )
				fXMax = x;
			if ( fYMin > y )
				fYMin = y;
			if ( fYMax < y )
				fYMax = y;
			if ( fZMin > z )
				fZMin = z;
			if ( fZMax < z )
				fZMax = z;
		}
		else
		{
			// Geant4 exception
			G4String msg = firstItem+" is not defined in the parser. Check the input file: "+fileName+".";
			G4Exception("PhysGeoImport::ReadDNAFabricFile()", "Geo_WrongParse", FatalException, msg);
		}
	}
	fOffsetX = (fXMin + fXMax)*0.5;
	fOffsetY = (fYMin + fYMax)*0.5;
	fOffsetZ = (fZMin + fZMax)*0.5;
	file.close();
}


G4VSolid* TsIRTPlasmidSupercoiled::CreateCutSolid(G4Orb *solidOrbRef,
											   TempMolecule &molRef,
											   std::vector<TempMolecule> &molList,
											   G4bool in)
{
	// The idea behing this method is to cut overlap volumes by selecting one of them (the reference) and checking all the other volumes (the targets).
	// If a reference and a target volumes are close enough to overlap they will be cut.
	// The reference is already selected when we enter this method.
	
	// Use the tiny space to differentiate the frontiers (may not be necessary)
	G4double tinySpace = 0.001*nm;
	
	// Cutted solid to be returned
	G4SubtractionSolid* solidCut(NULL);
	
	// Some flags
	G4bool isCutted = false;
	G4bool isOurVol = false;
	
	// Radius of the molecule to cut
	G4double radiusRef;
	if(molRef.fRadiusWater==0)
		radiusRef = molRef.fRadius;
	else
		radiusRef = molRef.fRadiusWater;
	
	// Reference volume position
	G4ThreeVector posRef = molRef.fPosition;
	
	// Look the other volumes of the voxel
	// Loop on all the target volumes (other volumes with potential overlaps)
	for(int i=0, ie=molList.size(); i<ie; ++i)
	{
		G4ThreeVector posTar = molList[i].fPosition;
		
		G4double rTar = posRef.z();
		G4double zTar = posTar.z();
		
		if(zTar>rTar+20*nm)
		{
			break;
		}
		else if(zTar<rTar-20*nm)
		{
			continue;
		}
		
		
		// Retrieve current target sphere informations
		G4double radiusTar;
		if(molList[i].fRadiusWater==0) radiusTar = molList[i].fRadius;
		else radiusTar = molList[i].fRadiusWater;
		
		// Compute the distance reference-target
		G4double distance = std::abs( (posRef - posTar).getR() );
		
		// Use the distance to check if the current target is also the reference.
		// This can only happen once per loop.
		if(distance==0 && !isOurVol)
		{
			// Target volume is also reference volume.
			
			// Set the flag
			isOurVol = true;
			
			// Next iteration
			continue;
		}
		// If the condition is correct more than one time then there is a mistake somewhere.
		else if(distance == 0 && isOurVol)
		{
			G4cerr<<"********************* Fatal Error **************************"<<G4endl;
			G4cerr<<"DetectorConstruction::CreateCutSolid: Two volumes are placed at the same position."<<G4endl;
			exit(EXIT_FAILURE);
		}
		
		// If the volumes are differents then we want to know if they are
		// close enough to overlap and, thus, to intiate a cut.
		else if(distance <= radiusRef+radiusTar)
		{
			// Volumes are close enough, there will be a cut
			
			// Box used to cut
			G4Box* solidBox = new G4Box("solid_box_for_cut", 2*radiusTar, 2*radiusTar, 2*radiusTar);
			
			// This part is tricky.
			// The goal is to calculate the position of the intersection center
			
			// diff vector to from ref to tar
			G4ThreeVector diff = posTar - posRef;
			
			// Find the intersection point and add to it half the length of the box used to cut
			G4double d = (pow(radiusRef,2)-pow(radiusTar,2)+pow(distance,2) ) / (2*distance) + solidBox->GetZHalfLength() - tinySpace;
			
			// If we are in another volume we choose to double the tinySpace to differentiate without ambiguities the inner and outer volume frontiers.
			// (may not be necessary)
			if(in) d -= 2*tinySpace;
			
			// Position of the box in order to achieve the cut.
			// "* ( diff/diff.getR() )" is necessary to get a vector in the right direction as output
			G4ThreeVector pos = d *( diff/diff.getR() );
			
			// Get the rotation angles because the box used to cut needs to be rotated
			// to give the right "cut angle".
			G4double phi = std::acos(pos.getZ()/pos.getR());
			G4double theta = std::acos( pos.getX() / ( pos.getR()*std::cos(M_PI/2.-phi) ) );
			
			if(pos.getY()<0) theta = -theta;
			
			G4ThreeVector rotAxisForPhi(1*nm,0.,0.);
			rotAxisForPhi.rotateZ(theta+M_PI/2);
			
			// Create the rotation matrix
			G4RotationMatrix *rotMat = new G4RotationMatrix;
			rotMat->rotate(-phi, rotAxisForPhi);
			
			// Rotate it again
			G4ThreeVector rotZAxis(0.,0.,1*nm);
			rotMat->rotate(theta, rotZAxis);
			
			// If the volume is cutted for the first time
			if(!isCutted) solidCut = new G4SubtractionSolid("solidCut", solidOrbRef, solidBox, rotMat, pos);
			
			// For the other times
			else solidCut = new G4SubtractionSolid("solidCut", solidCut, solidBox, rotMat, pos);
			
			// Set the cut flag
			isCutted = true;
		}
	}
	
	//delete rotMat;
	if(isCutted) return solidCut;
	
	// Otherwise, we return the original volume
	else return solidOrbRef;
}



