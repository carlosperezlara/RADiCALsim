#ifndef RADiCAL_h
#define RADiCAL_h 1

#include "G4VUserDetectorConstruction.hh"
#include "G4UserEventAction.hh"
#include "G4THitsMap.hh"
#include "G4UserRunAction.hh"
#include "G4VUserPrimaryGeneratorAction.hh"
#include "G4VUserActionInitialization.hh"
#include "globals.hh"

class G4VPhysicalVolume;
class G4GlobalMagFieldMessenger;
class G4Run;
class G4ParticleGun;
class G4Event;

namespace RAD {
  //=====================> G4UserDetectorConstruction <=================
  class DetectorConstruction : public G4VUserDetectorConstruction {
  public:
    DetectorConstruction() = default;
    ~DetectorConstruction() override = default;

  public:
    G4VPhysicalVolume* Construct() override;
    void ConstructSDandField() override;

  private:
    void DefineMaterials();
    G4VPhysicalVolume* DefineVolumes();
    
    static G4ThreadLocal G4GlobalMagFieldMessenger*  fMagFieldMessenger;
    G4bool fCheckOverlaps = true;
  };

  //======================> G4UserEventAction <=========================
  class EventAction : public G4UserEventAction {
  public:
    EventAction() = default;
    ~EventAction() override = default;
    
    void  BeginOfEventAction(const G4Event* event) override;
    void    EndOfEventAction(const G4Event* event) override;
    
  private:
    G4THitsMap<G4double>* GetHitsCollection(G4int hcID,
					    const G4Event* event) const;
    G4double GetSum(G4THitsMap<G4double>* hitsMap) const;
    void PrintEventStatistics(G4double absoEdep, G4double absoTrackLength,
			      G4double lysoEdep, G4double lysoTrackLength) const;
    
    G4int fAbsoEdepHCID = -1;
    G4int fLysoEdepHCID = -1;
    G4int fAbsoTrackLengthHCID = -1;
    G4int fLysoTrackLengthHCID = -1;
  };

  //=========================> G4UserRunAction <========================
  class RunAction : public G4UserRunAction {
  public:
    RunAction();
    ~RunAction() override = default;
    
    void BeginOfRunAction(const G4Run*) override;
    void   EndOfRunAction(const G4Run*) override;
  };

  //===================> G4UserPrimaryGeneratorAction <=================
  class PrimaryGeneratorAction : public G4VUserPrimaryGeneratorAction {
  public:
    PrimaryGeneratorAction();
    ~PrimaryGeneratorAction() override;

    void GeneratePrimaries(G4Event* event) override;
    void SetRandomFlag(G4bool value);

  private:
    G4ParticleGun* fParticleGun = nullptr; // G4 particle gun
  };

  //=====================> G4UserActionInitialization <=================
  class ActionInitialization : public G4VUserActionInitialization {
  public:
    ActionInitialization() = default;
    ~ActionInitialization() override = default;

    void BuildForMaster() const override;
    void Build() const override;
  };

}

#endif
