//++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
// PrimaryGenerator.cc
//++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
#include "PrimaryGenerator.hh"
#include "PltParticleGun.hh"
//#include "G4ParticleGun.hh"

//------------------------------------------------------------------------------
  PrimaryGenerator::PrimaryGenerator()
  : fpParticleGun(0) 
//------------------------------------------------------------------------------
{
  fpParticleGun = new PltParticleGun();
}

//------------------------------------------------------------------------------
  PrimaryGenerator::~PrimaryGenerator()
//------------------------------------------------------------------------------
{
  delete fpParticleGun;
}

//------------------------------------------------------------------------------
  void PrimaryGenerator::GeneratePrimaries(G4Event* anEvent)
//------------------------------------------------------------------------------
{
  fpParticleGun->GeneratePrimaryVertex(anEvent);
}

