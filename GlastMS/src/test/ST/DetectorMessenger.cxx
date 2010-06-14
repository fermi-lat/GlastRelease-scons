#include "DetectorMessenger.hh"

#include "DetectorConstruction.hh"
#include "G4UIdirectory.hh"
#include "G4UIcmdWithAString.hh"
#include "G4UIcmdWithAnInteger.hh"
#include "G4UIcmdWithADoubleAndUnit.hh"
#include "G4UIcmdWithoutParameter.hh"

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

DetectorMessenger::DetectorMessenger(DetectorConstruction *
Det)
:Detector(Det)
{ 
  detDir = new G4UIdirectory("/material/");
  detDir->SetGuidance("Material properties control.");
      
  MaterCmd = new G4UIcmdWithAString("/material/setMat",this);
  MaterCmd->SetGuidance("Select Material of the box.");
  MaterCmd->SetParameterName("choice",false);
  MaterCmd->AvailableForStates(G4State_Idle);
  
  ThickCmd = new G4UIcmdWithADoubleAndUnit("/material/setThick",this);
  ThickCmd->SetGuidance("Set Thickness of the box");
  ThickCmd->SetParameterName("Size",false);
  ThickCmd->SetRange("Size>=0.");
  ThickCmd->SetUnitCategory("Length");
  ThickCmd->AvailableForStates(G4State_Idle);
  
  SizeYCmd = new G4UIcmdWithADoubleAndUnit("/material/setSizeY",this);
  SizeYCmd->SetGuidance("Set Y size of the box");
  SizeYCmd->SetParameterName("Size",false);
  SizeYCmd->SetRange("Size>0.");
  SizeYCmd->SetUnitCategory("Length");    
  SizeYCmd->AvailableForStates(G4State_Idle);

  SizeZCmd = new G4UIcmdWithADoubleAndUnit("/material/setSizeZ",this);
  SizeZCmd->SetGuidance("Set Z size of the box");
  SizeZCmd->SetParameterName("Size",false);
  SizeZCmd->SetRange("Size>0.");
  SizeZCmd->SetUnitCategory("Length");    
  SizeZCmd->AvailableForStates(G4State_Idle);
  

  UpdateCmd = new G4UIcmdWithoutParameter("/material/update",this);
  UpdateCmd->SetGuidance("Update the box geometry.");
  UpdateCmd->SetGuidance("This command MUST be applied before \"beamOn\" ");
  UpdateCmd->SetGuidance("if you changed geometrical value(s).");
  UpdateCmd->AvailableForStates(G4State_Idle);
      
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

DetectorMessenger::~DetectorMessenger()
{
  delete MaterCmd;   
  delete ThickCmd; 
  delete SizeYCmd; 
  delete SizeZCmd;   
  delete UpdateCmd;
  delete detDir;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

void DetectorMessenger::SetNewValue(G4UIcommand* command,G4String
newValue)
{ 
  if( command == MaterCmd )
   { Detector->SetMaterial(newValue);}
  
  if( command == ThickCmd )
   { Detector->SetThickness(ThickCmd->GetNewDoubleValue(newValue));}

  if( command == SizeYCmd )
   { Detector->SetSizeY(SizeYCmd->GetNewDoubleValue(newValue));}

  if( command == SizeZCmd )
   { Detector->SetSizeZ(SizeZCmd->GetNewDoubleValue(newValue));}
  
  if( command == UpdateCmd )
   { Detector->UpdateGeometry(); }

}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....
