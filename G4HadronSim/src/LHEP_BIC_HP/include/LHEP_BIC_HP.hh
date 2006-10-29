//
// ********************************************************************
// * DISCLAIMER                                                       *
// *                                                                  *
// * The following disclaimer summarizes all the specific disclaimers *
// * of contributors to this software. The specific disclaimers,which *
// * govern, are listed with their locations in:                      *
// *   http://cern.ch/geant4/license                                  *
// *                                                                  *
// * Neither the authors of this software system, nor their employing *
// * institutes,nor the agencies providing financial support for this *
// * work  make  any representation or  warranty, express or implied, *
// * regarding  this  software system or assume any liability for its *
// * use.                                                             *
// *                                                                  *
// * This  code  implementation is the  intellectual property  of the *
// * GEANT4 collaboration.                                            *
// * By copying,  distributing  or modifying the Program (or any work *
// * based  on  the Program)  you indicate  your  acceptance of  this *
// * statement, and all its terms.                                    *
// ********************************************************************
//
// $Id$
// GEANT4 tag $Name$
//
//---------------------------------------------------------------------------
//
// ClassName:  LHEP_BIC_HP
//
// Author: 2002 J.P. Wellisch
//
// Modified:
// 30.11.2005 G.Folger: migration to non static particles
//
//----------------------------------------------------------------------------
//
#ifndef TLHEP_BIC_HP_h
#define TLHEP_BIC_HP_h 1

#include "G4VModularPhysicsList.hh"
#include "globals.hh"
#include "CompileTimeConstraints.hh"

template<class T>
class TLHEP_BIC_HP: public T
{
public:
  TLHEP_BIC_HP(G4int ver = 1);
  virtual ~TLHEP_BIC_HP();
  
public:
  // SetCuts() 
  virtual void SetCuts();
private:
  enum {ok = CompileTimeConstraints::IsA<T, G4VModularPhysicsList>::ok };

};
#include "LHEP_BIC_HP.icc"
typedef TLHEP_BIC_HP<G4VModularPhysicsList> LHEP_BIC_HP;

// 2002 by J.P. Wellisch

#endif



