#include "PltParticleGun.hh"
//#include "PltParticleGunMessenger.hh"

//#include "PltDetectorParameter.hh"

#include "G4Event.hh"
#include "G4ParticleGun.hh"
#include "G4ParticleTable.hh"
#include "G4ParticleDefinition.hh"
#include "G4RotationMatrix.hh"
#include "Randomize.hh"
#include "G4PhysicalConstants.hh"
#include "G4SystemOfUnits.hh"

using namespace CLHEP;

//PltParticleGun::PltParticleGun(const PltDetectorParameter* pDP)
PltParticleGun::PltParticleGun()
  :G4ParticleGun(1),
   particleTable(G4ParticleTable::GetParticleTable()),
   //pltDP(pDP),
   positionFlag(1), // default 0
   particleFlag(2), // default 0
   //cosThetaMax(-10./17.),
   cosThetaMax(-10./17.)
   //pMessenger(0)
{
  //pMessenger = new PltParticleGunMessenger(this);
}


PltParticleGun::~PltParticleGun()
{
  //delete pMessenger;
}



void PltParticleGun::GeneratePrimaryVertex(G4Event* anEvent)
{
 //G4ThreeVector vPos;
 // if (positionFlag==Target)        vPos = PutTraget();
 // else if (positionFlag==Source)   vPos = PutSource();
 // else                             vPos = particle_position;
 
	G4ThreeVector vPos = G4ThreeVector(0.0,0.0,0.0); 
 
  particleFlag = Gam3;

  // create a new vertex
  G4PrimaryVertex* vertex = 
    new G4PrimaryVertex(vPos,particle_time);

  if (particleFlag==Gam2) { 
    G4PrimaryParticle* particle[2]={0,0}; 
    GenerateTwoGamma(particle);
    vertex->SetPrimary( particle[0] );
    vertex->SetPrimary( particle[1] );

  } else if (particleFlag==Gam3) {

    G4PrimaryParticle* particle[3]={0,0,0}; 
    GenerateThreeGamma(particle);
    vertex->SetPrimary( particle[0] );
    vertex->SetPrimary( particle[1] );
    vertex->SetPrimary( particle[2] );

  } else if (particleFlag==Gam1275) {
    G4PrimaryParticle* particle[1]={0}; 
    GenerateGamma1275(particle);
    vertex->SetPrimary( particle[0] );

 } else if (particleFlag==Positron) {
    G4PrimaryParticle* particle[1]={0}; 
    GeneratePositron(particle);
    vertex->SetPrimary( particle[0] );

 } else {
    // create new primaries and set them to the vertex 
    //   (same as G4ParticleGun)
    G4double mass =  particle_definition->GetPDGMass();
    for( G4int i=0; i<NumberOfParticlesToBeGenerated; i++ ){
      G4PrimaryParticle* particle =
	new G4PrimaryParticle(particle_definition);
      particle->SetKineticEnergy( particle_energy );
      particle->SetMass( mass );
      particle->SetMomentumDirection( particle_momentum_direction );
      particle->SetCharge( particle_charge );
      particle->SetPolarization(particle_polarization.x(),
				particle_polarization.y(),
				particle_polarization.z());
      vertex->SetPrimary( particle );
    }
  } 

  anEvent->AddPrimaryVertex( vertex );

}

void PltParticleGun::GenerateTwoGamma(G4PrimaryParticle* gamma[2])
{
  G4ParticleDefinition* pD = particleTable->FindParticle("gamma");
  gamma[0] = new G4PrimaryParticle(pD);
  gamma[1] = new G4PrimaryParticle(pD);
  gamma[0]->SetMass(0.);
  gamma[1]->SetMass(0.);
  gamma[0]->SetCharge(0.);
  gamma[1]->SetCharge(0.);
  gamma[0]->SetKineticEnergy(electron_mass_c2);
  gamma[1]->SetKineticEnergy(electron_mass_c2);

  G4double px, py, pz;
  G4double cs, sn, phi;
  cs    =  RandFlat::shoot(-1.,1.);   
  sn    =  std::sqrt((1.0-cs)*(1.0+cs));   
  phi   =  RandFlat::shoot(0., CLHEP::twopi);   
  px    =  sn*std::cos(phi);
  py    =  sn*std::sin(phi);
  pz    =  cs;
  gamma[0]->SetMomentumDirection(G4ThreeVector(px, py, pz));
  gamma[0]->SetPolarization(G4ThreeVector(px, py, pz));
  gamma[1]->SetMomentumDirection(G4ThreeVector(-1.*px, -1.*py, -1.*pz));
  gamma[1]->SetPolarization(G4ThreeVector(px, py, pz));
}

void PltParticleGun::GenerateThreeGamma(G4PrimaryParticle* gamma[3])   
{
  G4ParticleDefinition* pD = particleTable->FindParticle("gamma");
  gamma[0] = new G4PrimaryParticle(pD);
  gamma[1] = new G4PrimaryParticle(pD);
  gamma[2] = new G4PrimaryParticle(pD);
  gamma[0]->SetMass(0.);
  gamma[1]->SetCharge(0.);
  gamma[2]->SetCharge(0.);

  // Determine 3 gamma energy  (based on Tanioka's program)
  const G4double eps=0.001;//w1 and w2 under line
  const G4double Sigmax=1.0;//sigma max
  G4double w0 = RandFlat::shoot(eps,1.0);//energy gamma[0]
  G4double w1 = RandFlat::shoot(eps,1.0);//energy gamma[1]
  G4double w2 = 2.0-w0-w1;                      //energy gamma[2]
  G4double z  = RandFlat::shoot(0.,Sigmax);
  while ( w2 >1. || z > sigma(w0,w1) ){
    w0 = RandFlat::shoot(eps,1.0);//energy gamma[0]
    w1 = RandFlat::shoot(eps,1.0);//energy gamma[1]
    w2 = 2.0-w0-w1;                      //energy gamma[2]
    z  = RandFlat::shoot(0.,Sigmax);
  }

  G4double theta =  std::acos(RandFlat::shoot(-1.,1.));
  G4double phi   =  RandFlat::shoot(0., CLHEP::twopi);   

  G4ThreeVector pMom0(1.0, 0.0, 0.0);
  pMom0.rotateY(theta);
  pMom0.rotateZ(phi);
  gamma[0]->SetKineticEnergy(electron_mass_c2*w0);
  gamma[0]->SetMomentumDirection(pMom0);

  G4double p1 = phi1(w0,w1); // angle between gamma[0] and gamma[1]
  G4ThreeVector pMom1(std::cos(p1), std::sin(p1), 0.0);
  pMom1.rotateY(theta);
  pMom1.rotateZ(phi);
  gamma[1]->SetKineticEnergy(electron_mass_c2*w1);
  gamma[1]->SetMomentumDirection(pMom1);

  G4double p2 = phi2(w0,w1); // angle between gamma[0] and gamma[2]
  G4ThreeVector pMom2(std::cos(p2), std::sin(p2), 0.0);
  pMom2.rotateY(theta);
  pMom2.rotateZ(phi);
  gamma[2]->SetKineticEnergy(electron_mass_c2*w2);
  gamma[2]->SetMomentumDirection(pMom2);

  //G4cout << w0 << "( " << w0 <<","<< w0*0.<< ")" << G4endl;
  //G4cout << w1 << "( " << w1*std::cos(p1) <<","<< w1*std::sin(p1) << ")" << G4endl;
  //G4cout << w2 << "( " << w2*std::cos(p2) <<","<< w2*std::sin(p2) << ")" << G4endl;
  //G4cout << "(" << theta <<","<< phi <<")"<< G4endl;
  
}


void PltParticleGun::GenerateGamma1275(G4PrimaryParticle* gamma[1])
{
  G4ParticleDefinition* pD = particleTable->FindParticle("gamma");
  gamma[0] = new G4PrimaryParticle(pD);
  gamma[0]->SetMass(0.);
  gamma[0]->SetCharge(0.);
  gamma[0]->SetKineticEnergy(1274.6*keV);
  G4double px, py, pz;
  G4double cs, sn, phi;
  cs    =  RandFlat::shoot(-1.,0.);  // only downward    
  sn    =  std::sqrt((1.0-cs)*(1.0+cs));   
  phi   =  RandFlat::shoot(0., CLHEP::twopi);   
  px    =  sn*std::cos(phi);
  py    =  sn*std::sin(phi);
  pz    =  cs;
  gamma[0]->SetMomentumDirection(G4ThreeVector(px, py, pz));
			 
}

void PltParticleGun::GeneratePositron(G4PrimaryParticle* positron[1])  
{
  G4cout <<"GeneratePositron : cosThetaMax=" << cosThetaMax << G4endl;

  G4ParticleDefinition* pD = particleTable->FindParticle("e+");
  positron[0] = new G4PrimaryParticle(pD);
  positron[0]->SetMass(electron_mass_c2);
  positron[0]->SetCharge(1.);
  G4double px, py, pz;
  G4double cs, sn, phi;
  cs    =  RandFlat::shoot(-1.,cosThetaMax);  // only downward    
  sn    =  std::sqrt((1.0-cs)*(1.0+cs));   
  phi   =  RandFlat::shoot(0., CLHEP::twopi);   
  px    =  sn*std::cos(phi);
  py    =  sn*std::sin(phi);
  pz    =  cs;
  positron[0]->SetMomentumDirection(G4ThreeVector(px, py, pz)); // down

  const G4double emax = 0.543 * MeV; 
  const G4double pmax = std::sqrt((2*electron_mass_c2+emax)*emax);
  G4double pep=RandFlat::shoot(0.,pmax );
  G4double prob = RandFlat::shoot(0.,1.5*beta(pmax*0.5));
  while (  prob > beta(pep) ){
    pep=RandFlat::shoot(0.,pmax);
    prob = RandFlat::shoot(0.,1.5*beta(pmax*0.5));
  }
  positron[0]->SetKineticEnergy(std::sqrt(electron_mass_c2*electron_mass_c2+pep*pep)-electron_mass_c2); // kinetic energy
}

const G4ThreeVector& PltParticleGun::PutSource()
{
  static G4ThreeVector vPos(0.,0.,0.);

//  G4double x = 9999.;
//  G4double y = 9999.;
//  while (x*x+y*y > pltDP->Na22R()*pltDP->Na22R()) {
//    x = RandFlat::shoot(-1.*pltDP->Na22R(),pltDP->Na22R());
//    y = RandFlat::shoot(-1.*pltDP->Na22R(),pltDP->Na22R());
//  }
//  G4double z = RandFlat::shoot(pltDP->Na22Z()-pltDP->Na22Thick()/2.,
//			       pltDP->Na22Z()+pltDP->Na22Thick()/2.-pltDP->Na22CoverThick());

  G4double x = 0.0*mm;
  G4double y = 0.0*mm;
  G4double z = 0.0*mm;

  vPos.setX(x);
  vPos.setY(y);
  vPos.setZ(z);

  return vPos;  
}   

const G4ThreeVector& PltParticleGun::PutTraget()
{
  static G4ThreeVector vPos(0.,0.,0.);

//  G4double x = 9999.;
//  G4double y = 9999.;
//  while (x*x+y*y > pltDP->TargetR()*pltDP->TargetR()) {
//    x = RandFlat::shoot(-1.*pltDP->TargetR(),pltDP->TargetR());
//    y = RandFlat::shoot(-1.*pltDP->TargetR(),pltDP->TargetR());
//  }
//  G4double z = RandFlat::shoot(pltDP->TargetZ()-pltDP->TargetThick()/2.,
//			       pltDP->TargetZ()+pltDP->TargetThick()/2.);

  G4double x = 0.0*mm;
  G4double y = 0.0*mm;
  G4double z = 0.0*mm;

  vPos.setX(x);
  vPos.setY(y);
  vPos.setZ(z);

  return vPos;  
}



