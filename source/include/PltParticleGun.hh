#ifndef PltParticleGun_h
#define PltParticleGun_h 1

#include "G4VUserPrimaryGeneratorAction.hh"
#include "globals.hh"

#include "G4ParticleGun.hh"
#include "G4PhysicalConstants.hh"
#include "G4SystemOfUnits.hh"

class G4ParticleTable;
class G4Event;
//class PltDetectorParameter;
//class PltParticleGunMessenger;

class PltParticleGun : public G4ParticleGun
{
 friend class  PltParticleGunMessenger;
 public:
  //PltParticleGun(const PltDetectorParameter*);    
  PltParticleGun();    
  ~PltParticleGun();
  
 public:
  virtual void GeneratePrimaryVertex(G4Event*);

 protected:  
  void GenerateTwoGamma(G4PrimaryParticle* gamma[2]);   
  void GenerateThreeGamma(G4PrimaryParticle* gamma[3]);   
  void GenerateGamma1275(G4PrimaryParticle* gamma[1]);   
  void GeneratePositron(G4PrimaryParticle* positron[1]);   

  const G4ThreeVector& PutSource();   
  const G4ThreeVector& PutTraget();   

 private:  
  G4double sigma(G4double x,G4double y) const {
    return (x+y-1.)*(x+y-1.)/(x*x*y*y);
  }

  G4double phi1(G4double x,G4double y) const {
    return std::acos((2.-2.*x -2.*y + x*y)/(x*y));
  }

  G4double phi2(G4double x,G4double y) const {
    return -1.*std::acos(-1.*(x*x + x*y -2.*x -2.*y + 2.)/(x*(2.-x-y)));
  }
 
  G4double beta(G4double x){
    const G4double me2=electron_mass_c2 * electron_mass_c2;
    const G4double emax=543 * keV;
    const G4double p2=(2.*electron_mass_c2+emax)*emax;
    G4double val = (std::sqrt(me2+p2)-std::sqrt(me2+x*x))*x;
    return val*val;
  }

 private:
  G4ParticleTable*                particleTable;
  //const PltDetectorParameter*     pltDP;
  
  G4int positionFlag;
  enum { UserPos=0, Target, Source};
  G4int particleFlag;
  enum { User=0, Gam2, Gam3, Gam1275, Positron};
  
  G4double                     cosThetaMax; 
  //PltParticleGunMessenger*     pMessenger;
  
};


#endif


