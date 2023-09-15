#include "RADiCAL.hh"

#include "G4Material.hh"
#include "G4NistManager.hh"

#include "G4Box.hh"
#include "G4LogicalVolume.hh"
#include "G4PVPlacement.hh"
#include "G4PVReplica.hh"
#include "G4GlobalMagFieldMessenger.hh"
#include "G4AutoDelete.hh"

#include "G4SDManager.hh"
#include "G4SDChargedFilter.hh"
#include "G4MultiFunctionalDetector.hh"
#include "G4VPrimitiveScorer.hh"
#include "G4PSEnergyDeposit.hh"
#include "G4PSTrackLength.hh"

#include "G4VisAttributes.hh"
#include "G4Colour.hh"

#include "G4PhysicalConstants.hh"
#include "G4SystemOfUnits.hh"

#include "G4AnalysisManager.hh"
#include "G4RunManager.hh"
#include "G4Event.hh"
#include "G4SDManager.hh"
#include "G4HCofThisEvent.hh"
#include "G4UnitsTable.hh"

#include "Randomize.hh"
#include <iomanip>

#include "G4Run.hh"

#include "G4LogicalVolumeStore.hh"
#include "G4ParticleGun.hh"
#include "G4ParticleTable.hh"
#include "G4ParticleDefinition.hh"
#include "Randomize.hh"


namespace RAD {
  //=============================================================================
  //============================> DetectorConstruction <=========================
  G4ThreadLocal
  G4GlobalMagFieldMessenger* DetectorConstruction::fMagFieldMessenger = nullptr;
  //=============================================================================
  G4VPhysicalVolume* DetectorConstruction::Construct() {
    DefineMaterials();
    return DefineVolumes();
  }
  //=============================================================================
  void DetectorConstruction::DefineMaterials() {
    auto nistManager = G4NistManager::Instance();
    nistManager->FindOrBuildMaterial("G4_W");

    // Optical Properties
    const G4int num = 20;
    G4double ene[num]   =  {1.79*eV, 1.85*eV, 1.91*eV, 1.97*eV,
			    2.04*eV, 2.11*eV, 2.19*eV, 2.27*eV,
			    2.36*eV, 2.45*eV, 2.56*eV, 2.67*eV,
			    2.80*eV, 2.94*eV, 3.09*eV, 3.25*eV,
			    3.44*eV, 3.65*eV, 3.89*eV, 4.16*eV};
    G4double rAir[num]  =  {1.00, 1.00, 1.00, 1.00,
			    1.00, 1.00, 1.00, 1.00,
			    1.00, 1.00, 1.00, 1.00,
			    1.00, 1.00, 1.00, 1.00,
			    1.00, 1.00, 1.00, 1.00};
    //G4MaterialPropertiesTable* mpt2 = new G4MaterialPropertiesTable();
    //mpt2->AddProperty("RINDEX" , ene, rAir , num);

    G4double fast[num]  =  {0.01, 0.10, 0.20, 0.50,
			    0.90, 1.70, 2.90, 5.00,
			    8.30, 12.5, 17.0, 22.9,
			    26.4, 25.6, 16.8, 4.20,
			    0.30, 0.20, 0.10, 0.01};
    G4double rLyso[num] =  {1.81, 1.81, 1.81, 1.81,
			    1.81, 1.81, 1.81, 1.81,
			    1.81, 1.81, 1.81, 1.81,
			    1.81, 1.81, 1.81, 1.81,
			    1.81, 1.81, 1.81, 1.81};
    G4double abs[num]   =  {3.5*m, 3.5*m, 3.5*m, 3.5*m,
			    3.5*m, 3.5*m, 3.5*m, 3.5*m,
			    3.5*m, 3.5*m, 3.5*m, 3.5*m,
			    3.5*m, 3.5*m, 3.5*m, 3.5*m,
			    3.5*m, 3.5*m, 3.5*m, 3.5*m};
    //G4MaterialPropertiesTable* mpt = new G4MaterialPropertiesTable();
    //mpt->AddProperty("FASTCOMPONENT", ene, fast, num);
    //mpt->AddProperty("RINDEX", ene, rLyso , num);
    //mpt->AddProperty("ABSLENGTH", ene, abs, num);
    //mpt->AddConstProperty("SCINTILLATIONYIELD",32/keV);
    //mpt->AddConstProperty("RESOLUTIONSCALE", 1);
    //mpt->AddConstProperty("FASTTIMECONSTANT",41*ns);


    // Liquid argon material
    G4double a;  // mass of a mole;
    G4double z;  // z=mean number of protons;
    G4double density;
    new G4Material("liquidArgon", z=18., a= 39.95*g/mole, density= 1.390*g/cm3);
    // The argon by NIST Manager is a gas with a different density

    // Vacuum material
    G4Material* vac = new G4Material("Galactic", z=1., a=1.01*g/mole,density= universe_mean_density,
				     kStateGas, 2.73*kelvin, 3.e-18*pascal);
    //vac->SetMaterialPropertiesTable(mpt2);

    // LYSO material
    G4double prelude_density = 7.4*g/cm3;
    G4Material* prelude = new G4Material("prelude", prelude_density, 4);
    prelude->AddElement(nistManager->FindOrBuildElement("Lu"), 71*perCent);
    prelude->AddElement(nistManager->FindOrBuildElement("Si"),  7*perCent);
    prelude->AddElement(nistManager->FindOrBuildElement("O"),  18*perCent);
    prelude->AddElement(nistManager->FindOrBuildElement("Y"),   4*perCent);
    G4Material* lyso = new G4Material("G4_lyso",    prelude_density, 2);
    lyso->AddMaterial(prelude,                       99.81*perCent);
    lyso->AddElement(nistManager->FindOrBuildElement("Ce"),  0.19*perCent);
    //lyso->SetMaterialPropertiesTable(mpt);

    // Print materials
    G4cout << *(G4Material::GetMaterialTable()) << G4endl;
  }
  //=============================================================================
  G4VPhysicalVolume* DetectorConstruction::DefineVolumes() {
    // Geometry parameters
    G4int    nofLayers = 28; // number of lyso layers
    G4double absoThickness =  2.5*mm;
    G4double lysoThickness =  1.5*mm;
    G4double calorSizeXY   = 14.0*mm;

    auto  layerThickness = absoThickness + lysoThickness;
    auto  calorThickness = nofLayers * layerThickness + absoThickness; //ends in W
    auto  worldSizeXY = 1.2 * calorSizeXY;
    auto  worldSizeZ  = 1.2 * calorThickness;

    auto defaultMaterial  = G4Material::GetMaterial("Galactic");
    auto absorberMaterial = G4Material::GetMaterial("G4_W");
    //auto lysoMaterial     = G4Material::GetMaterial("liquidArgon");
    auto lysoMaterial     = G4Material::GetMaterial("G4_lyso");
    //auto lysoMaterial     = G4Material::GetMaterial("G4_W");
    
    if ( ! defaultMaterial || ! absorberMaterial || ! lysoMaterial ) {
      G4ExceptionDescription msg;
      msg << "Cannot retrieve materials already defined.";
      G4Exception("DetectorConstruction::DefineVolumes()",
		  "MyCode0001", FatalException, msg);
    }

    auto worldS = new G4Box("World", worldSizeXY/2, worldSizeXY/2, worldSizeZ/2);
    auto worldLV = new G4LogicalVolume(worldS, defaultMaterial, "World");
    auto worldPV = new G4PVPlacement(nullptr,                  // no rotation
				     G4ThreeVector(),          // at (0,0,0)
				     worldLV,                  // its logical volume
				     "World",                  // its name
				     nullptr,                  // its mother  volume
				     false,                    // no boolean operation
				     0,                        // copy number
				     fCheckOverlaps);          // checking overlaps
    
    auto calorimeterS = new G4Box("Calorimeter", calorSizeXY/2, calorSizeXY/2, calorThickness/2);
    auto calorLV = new G4LogicalVolume(calorimeterS, defaultMaterial, "Calorimeter");
    new G4PVPlacement(nullptr,                  // no rotation
		      G4ThreeVector(),          // at (0,0,0)
		      calorLV,                  // its logical volume
		      "Calorimeter",            // its name
		      worldLV,                  // its mother  volume
		      false,                    // no boolean operation
		      0,                        // copy number
		      fCheckOverlaps);          // checking overlaps

    auto absorberS = new G4Box("Abso", calorSizeXY/2, calorSizeXY/2, absoThickness/2);
    auto absorberLV = new G4LogicalVolume(absorberS, absorberMaterial, "AbsoLV");

    auto lysoS = new G4Box("LYSO", calorSizeXY/2, calorSizeXY/2, lysoThickness/2);
    auto lysoLV = new G4LogicalVolume(lysoS, lysoMaterial, "LYSOLV");

    if(0) {
      auto layerS = new G4Box("Layer", calorSizeXY/2, calorSizeXY/2, layerThickness/2);
      auto layerLV = new G4LogicalVolume(layerS, defaultMaterial, "Layer");
    
      new G4PVPlacement(nullptr,                                   // no rotation
			G4ThreeVector(0., 0., -lysoThickness / 2), // its position
			absorberLV,                                // its logical volume
			"Abso",                                    // its name
			layerLV,                                   // its mother  volume
			false,                                     // no boolean operation
			0,                                         // copy number
			fCheckOverlaps);                           // checking overlaps
      
      new G4PVPlacement(nullptr,                                   // no rotation
			G4ThreeVector(0., 0., absoThickness / 2),  // its position
			lysoLV,                                    // its logical volume
			"LYSO",                                    // its name
			layerLV,                                   // its mother  volume
			false,                                     // no boolean operation
			0,                                         // copy number
			fCheckOverlaps);                           // checking overlaps
      
      new G4PVReplica("Layer",          // its name
		      layerLV,          // its logical volume
		      calorLV,          // its mother
		      kZAxis,           // axis of replication
		      nofLayers,        // number of replica
		      layerThickness);  // witdth of replica
    } else {
      G4double startZ = -calorThickness/2;
      for(int i=0; i!=nofLayers+1; ++i) {
	G4double placeAbs = startZ + i * layerThickness;
	G4double placeLys = startZ + absoThickness + i * layerThickness;
	new G4PVPlacement(nullptr,                                   // no rotation
			  G4ThreeVector(0., 0., placeAbs),           // its position
			  absorberLV,                                // its logical volume
			  "Abso",                                    // its name
			  calorLV,                                   // its mother  volume
			  false,                                     // no boolean operation
			  i,                                         // copy number
			  fCheckOverlaps);                           // checking overlaps
	if( i == nofLayers ) continue;
	new G4PVPlacement(nullptr,                                   // no rotation
			  G4ThreeVector(0., 0., placeLys),           // its position
			  lysoLV,                                    // its logical volume
			  "LYSO",                                    // its name
			  calorLV,                                   // its mother  volume
			  false,                                     // no boolean operation
			  i,                                         // copy number
			  fCheckOverlaps);                           // checking overlaps
      }
    }

    G4cout
      << G4endl
      << "------------------------------------------------------------" << G4endl
      << "---> The calorimeter is " << nofLayers << " layers of: [ "
      << absoThickness/mm << "mm of " << absorberMaterial->GetName()
      << " + "
      << lysoThickness/mm << "mm of " << lysoMaterial->GetName() << " ] " << G4endl
      << "------------------------------------------------------------" << G4endl;
    
    worldLV->SetVisAttributes (G4VisAttributes::GetInvisible());
    
    auto simpleBoxVisAtt= new G4VisAttributes(G4Colour(1.0,1.0,1.0));
    simpleBoxVisAtt->SetVisibility(true);
    calorLV->SetVisAttributes(simpleBoxVisAtt);
    
    // Always return the physical World
    return worldPV;
  }
  //=============================================================================
  void DetectorConstruction::ConstructSDandField() {
    G4SDManager::GetSDMpointer()->SetVerboseLevel(1);

    // Scorers
    G4VPrimitiveScorer* primitive;
    auto charged = new G4SDChargedFilter("chargedFilter");

    // declare Absorber as a MultiFunctionalDetector scorer
    auto absDetector = new G4MultiFunctionalDetector("AbsorberDet");
    primitive = new G4PSEnergyDeposit("Edep"); // create a Prim Scorer for EneDep
    absDetector->RegisterPrimitive(primitive);
    primitive = new G4PSTrackLength("TrackLength"); // create a Prim Scorer for TrackLength
    primitive ->SetFilter(charged);
    absDetector->RegisterPrimitive(primitive);
    G4SDManager::GetSDMpointer()->AddNewDetector(absDetector);
    SetSensitiveDetector("AbsoLV",absDetector);

    // declare Lyso as a MultiFunctionalDetector scorer
    auto lysoDetector = new G4MultiFunctionalDetector("LYSODet");
    primitive = new G4PSEnergyDeposit("Edep");
    lysoDetector->RegisterPrimitive(primitive);
    primitive = new G4PSTrackLength("TrackLength");
    primitive ->SetFilter(charged);
    lysoDetector->RegisterPrimitive(primitive);
    G4SDManager::GetSDMpointer()->AddNewDetector(lysoDetector);
    SetSensitiveDetector("LYSOLV",lysoDetector);
    
    // Magnetic field
    G4ThreeVector fieldValue;
    fMagFieldMessenger = new G4GlobalMagFieldMessenger(fieldValue);
    fMagFieldMessenger->SetVerboseLevel(1);
    G4AutoDelete::Register(fMagFieldMessenger);
  }
  //=============================================================================
  //================================> EventAction <==============================
  G4THitsMap<G4double>* EventAction::GetHitsCollection(G4int hcID,
						       const G4Event* event) const {
    auto hitsCollection = static_cast<G4THitsMap<G4double>*>(event->GetHCofThisEvent()->GetHC(hcID));
    if ( ! hitsCollection ) {
      G4ExceptionDescription msg;
      msg << "Cannot access hitsCollection ID " << hcID;
      G4Exception("EventAction::GetHitsCollection()",
		  "MyCode0003", FatalException, msg);
    }
    return hitsCollection;
  }
  //=============================================================================
  G4double EventAction::GetSum(G4THitsMap<G4double>* hitsMap) const {
    G4double sumValue = 0.;
    for ( auto it : *hitsMap->GetMap() ) {
      // hitsMap->GetMap() returns the map of std::map<G4int, G4double*>
      sumValue += *(it.second);
    }
    return sumValue;
  }
  //=============================================================================
  void EventAction::PrintEventStatistics(G4double absoEdep, G4double absoTrackLength,
					 G4double lysoEdep, G4double lysoTrackLength) const {
    G4cout
      << "   Absorber: total energy: "
      << std::setw(7) << G4BestUnit(absoEdep, "Energy")
      << "       total track length: "
      << std::setw(7) << G4BestUnit(absoTrackLength, "Length")
      << G4endl
      << "        Lyso: total energy: "
      << std::setw(7) << G4BestUnit(lysoEdep, "Energy")
      << "       total track length: "
      << std::setw(7) << G4BestUnit(lysoTrackLength, "Length")
      << G4endl;
  }
  //=============================================================================
  void EventAction::BeginOfEventAction(const G4Event* /*event*/) {}
  //=============================================================================
  void EventAction::EndOfEventAction(const G4Event* event) {
    // Get hist collections IDs
    if ( fAbsoEdepHCID == -1 ) {
      fAbsoEdepHCID = G4SDManager::GetSDMpointer()->GetCollectionID("AbsorberDet/Edep");
      fLysoEdepHCID = G4SDManager::GetSDMpointer()->GetCollectionID("LYSODet/Edep");
      fAbsoTrackLengthHCID = G4SDManager::GetSDMpointer()->GetCollectionID("AbsorberDet/TrackLength");
      fLysoTrackLengthHCID = G4SDManager::GetSDMpointer()->GetCollectionID("LYSODet/TrackLength");
    }
    // Get sum values from hits collections
    auto absoEdep = GetSum(GetHitsCollection(fAbsoEdepHCID, event));
    auto lysoEdep = GetSum(GetHitsCollection(fLysoEdepHCID, event));
    auto absoTrackLength = GetSum(GetHitsCollection(fAbsoTrackLengthHCID, event));
    auto lysoTrackLength = GetSum(GetHitsCollection(fLysoTrackLengthHCID, event));
    
    // get analysis manager
    auto analysisManager = G4AnalysisManager::Instance();
    
    // fill histograms
    analysisManager->FillH1(0, absoEdep);
    analysisManager->FillH1(1, lysoEdep);
    analysisManager->FillH1(2, absoTrackLength);
    analysisManager->FillH1(3, lysoTrackLength);

    // fill profiles
    auto lysoHits = *GetHitsCollection(fLysoEdepHCID, event)->GetMap();
    for ( auto it : lysoHits ) {
      analysisManager->FillP1( 0, it.first, *(it.second) );
    }
    
    // fill ntuple
    analysisManager->FillNtupleDColumn(0, absoEdep);
    analysisManager->FillNtupleDColumn(1, lysoEdep);
    analysisManager->FillNtupleDColumn(2, absoTrackLength);
    analysisManager->FillNtupleDColumn(3, lysoTrackLength);
    analysisManager->AddNtupleRow();
    
    //print per event (modulo n)
    auto eventID = event->GetEventID();
    auto printModulo = G4RunManager::GetRunManager()->GetPrintProgress();
    if ( ( printModulo > 0 ) && ( eventID % printModulo == 0 ) ) {
      G4cout << "---> End of event: " << eventID << G4endl;
      PrintEventStatistics(absoEdep, absoTrackLength, lysoEdep, lysoTrackLength);
    }
  }
  //=============================================================================
  //================================> RunAction <================================
  RunAction::RunAction() {
    // set printing event number per each event
    G4RunManager::GetRunManager()->SetPrintProgress(1);

    // Create analysis manager
    auto analysisManager = G4AnalysisManager::Instance();
    // Create directories
    analysisManager->SetHistoDirectoryName("histograms");
    analysisManager->SetNtupleDirectoryName("ntuple");
    analysisManager->SetVerboseLevel(1);
    analysisManager->SetNtupleMerging(true);
    // Note: merging ntuples is available only with Root output
    
    // Book histograms, ntuple
    // Creating histograms
    analysisManager->CreateH1("Eabs","Edep in absorber",   100, 0., 120000*MeV);
    analysisManager->CreateH1("Elyso","Edep in lyso",      100, 0., 120000*MeV);
    analysisManager->CreateH1("Labs","trackL in absorber", 100, 0., 50*cm);
    analysisManager->CreateH1("Llyso","trackL in lyso",    100, 0., 50*cm);
    // Creating profiles
    analysisManager->CreateP1("LysoAbs","Edep in Lyso Layer", 30, -0.5, 29.5);
    
    // Creating ntuple
    analysisManager->CreateNtuple("RADiCAL", "Edep and TrackL");
    analysisManager->CreateNtupleDColumn("Eabs");
    analysisManager->CreateNtupleDColumn("Elyso");
    analysisManager->CreateNtupleDColumn("Labs");
    analysisManager->CreateNtupleDColumn("Llyso");
    analysisManager->FinishNtuple();
  }
  //=============================================================================
  void RunAction::BeginOfRunAction(const G4Run* /*run*/) {
    //inform the runManager to save random number seed
    //G4RunManager::GetRunManager()->SetRandomNumberStore(true);
    
    // Get analysis manager
    auto analysisManager = G4AnalysisManager::Instance();
    
    // Open an output file
    G4String fileName = "RADout.root";
    // Other supported output types:
    // G4String fileName = "RADout.csv";
    // G4String fileName = "RADout.hdf5";
    // G4String fileName = "RADout.xml";
    analysisManager->OpenFile(fileName);
    G4cout << "Using " << analysisManager->GetType() << G4endl;
  }
  //=============================================================================
  void RunAction::EndOfRunAction(const G4Run* /*run*/) {
    // print histogram statistics
    auto analysisManager = G4AnalysisManager::Instance();
    if ( analysisManager->GetH1(1) ) {
      G4cout << G4endl << " ----> print histograms statistic ";
      if(isMaster) {
	G4cout << "for the entire run " << G4endl << G4endl;
      }
      else {
	G4cout << "for the local thread " << G4endl << G4endl;
      }
      
      G4cout << " EAbs : mean = "
	     << G4BestUnit(analysisManager->GetH1(0)->mean(), "Energy")
	     << " rms = "
	     << G4BestUnit(analysisManager->GetH1(0)->rms(),  "Energy") << G4endl;
      
      G4cout << " ELyso : mean = "
	     << G4BestUnit(analysisManager->GetH1(1)->mean(), "Energy")
	     << " rms = "
	     << G4BestUnit(analysisManager->GetH1(1)->rms(),  "Energy") << G4endl;
      
      G4cout << " LAbs : mean = "
	     << G4BestUnit(analysisManager->GetH1(2)->mean(), "Length")
	     << " rms = "
	     << G4BestUnit(analysisManager->GetH1(2)->rms(),  "Length") << G4endl;
      
      G4cout << " LLyso : mean = "
	     << G4BestUnit(analysisManager->GetH1(3)->mean(), "Length")
	     << " rms = "
	     << G4BestUnit(analysisManager->GetH1(3)->rms(),  "Length") << G4endl;
    }
    
    // save histograms & ntuple
    analysisManager->Write();
    analysisManager->CloseFile();
  }
  //=============================================================================
  //=========================> PrimaryActionGenerator <==========================
  PrimaryGeneratorAction::PrimaryGeneratorAction() {
    G4int nofParticles = 1;
    fParticleGun = new G4ParticleGun(nofParticles);

    // default particle kinematic
    auto particleDefinition = G4ParticleTable::GetParticleTable()->FindParticle("e-");
    fParticleGun->SetParticleDefinition(particleDefinition);
    fParticleGun->SetParticleMomentumDirection(G4ThreeVector(0.,0.,1.));
    fParticleGun->SetParticleEnergy(120.*GeV);
  }
  //=============================================================================
  PrimaryGeneratorAction::~PrimaryGeneratorAction() {
    delete fParticleGun;
  }
  //=============================================================================
  void PrimaryGeneratorAction::GeneratePrimaries(G4Event* anEvent) {
    // This function is called at the begining of event
    // In order to avoid dependence of PrimaryGeneratorAction
    // on DetectorConstruction class we get world volume
    // from G4LogicalVolumeStore
    G4double worldZHalfLength = 0.;
    auto worldLV = G4LogicalVolumeStore::GetInstance()->GetVolume("World");
    
    // Check that the world volume has box shape
    G4Box* worldBox = nullptr;
    if ( worldLV ) {
      worldBox = dynamic_cast<G4Box*>(worldLV->GetSolid());
    }
    if ( worldBox ) {
      worldZHalfLength = worldBox->GetZHalfLength();
    } else {
      G4ExceptionDescription msg;
      msg << "World volume of box shape not found." << G4endl;
      msg << "Perhaps you have changed geometry." << G4endl;
      msg << "The gun will be place in the center.";
      G4Exception("PrimaryGeneratorAction::GeneratePrimaries()",
		  "MyCode0002", JustWarning, msg);
    }
    
    // Set gun position
    fParticleGun->SetParticlePosition(G4ThreeVector(0., 0., -worldZHalfLength));
    fParticleGun->GeneratePrimaryVertex(anEvent);
  }
  //=============================================================================
  //==========================> ActionInitialization <===========================
  void ActionInitialization::BuildForMaster() const {
    SetUserAction(new RunAction);
  }
  //=============================================================================
  void ActionInitialization::Build() const {
    SetUserAction(new PrimaryGeneratorAction);
    SetUserAction(new RunAction);
    SetUserAction(new EventAction);
  }

}
