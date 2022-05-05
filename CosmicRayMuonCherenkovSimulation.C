//Cherenkov Radiation Emission from Atmospheric Muons.
//John Campbell-08408891, summer 2011.


# ifndef __CINT__ 



#include <fstream>
#include <string>
#include<vector>  
#include<iostream>
#include <stdio.h>
#include <stddef.h>
#include <stdlib.h>
#include <ctype.h>
#include <TTree.h>
#include <TRandom3.h>
#include <TStyle.h>
#include <TROOT.h>
#include<TStopwatch.h>
#include<TH1.h>
#include<TH2.h>
#include<TH2Poly.h>
#include<TFile.h>
#include<TMath.h>  
#include<TNtuple.h>  
#include<TCanvas.h>
#include<TRandom.h>
#include<TFormula.h>
#include<TGraph.h>
#include<TPad.h>


#endif   


using namespace std;


int PoissonRandomNumber(const double lambda)
{
if (lambda==0) return 0;

int k =0; 
 const int max_k = 1000; // upper limit
double p = gRandom->Uniform (0,1); //uniform random number
double P = exp(-lambda); // probability
double sum=P ; // cumulant
if (sum>=p ) return 0; // done
for (k =1; k<max_k ; ++k ) { // Loop over k's
P*=lambda/(double)k; // Calculate next probability
sum+=P ; // Increase cumulant
if (sum>=p) break; // Leave loop
}
return k; // return random number
}
double DEG( const double xrad){
return xrad*(180./TMath::Pi());
}
double RAD( const double xdeg ){
return xdeg*(TMath::Pi()/180.);
}

static const double CritEnergy = 8.4e-2; // 84 MeV

int AP2 = 1;

int MUON = 1 ; //Simulate a single Muon



//Simulate Bremstrahlung Photons using a probability amplitude vertex function for 1 and 2 loop photon processes.


double EnergySplitting(double E){

  
if (AP2==1) return E ;
else if (AP2==2){
     
double SOLV=0;
for (int i= 0; ; ){
double u = gRandom->Uniform(0,1); // returns a random number between 0 and 1
double u2 = pow(u,2);
double SOLV=(0.5 + 0.132283 *pow((648.*pow((u2 -(1./3.)*u+0.0347222),0.5)+648.*u
-108.),(1./3.))-1.88988/pow((648.*pow((u2-(1./3.)*u+0.0347222),0.5)+648.*u
-108.),(1./3.)));
if (SOLV<=1 && SOLV>=0) return SOLV*E;
}
}
}


//calculate max radiation length


int MaxRadiation(double InitEnergy){
return TMath::Nint(TMath::Log(InitEnergy/CritEnergy)/TMath::Log(2));
}


//calculate height at start of muon path

double x0(double hstart)
{
const double rho0=1.24e-3; //Volume Density of Atmosphere in gcm-3
const double xicm = 1./8.33e5 ; // cm-1
const double xi = 1./8.33 ; // km
const double mu = 37.3 ; // Mean Free path in gcm-2
double a = mu*xicm/rho0;
double b = TMath::Exp(-xi * hstart);
double c = pow(xi,-1) + hstart;

return TMath::Log(a+b)*pow(xi,-1) + hstart; //return in km
}

int NumberOfElectrons(int step){
int e = 0 ;
if (step == 0) return 0;
if (step == 1) return 2;
if (step == 2) return 2;
else if (step > 2) {
     
for (unsigned int i = 3; i<=step; i++){
    
e = pow(2,i) - e;
}
}
return e;
}

int NumberOfPhotons (int step){
int photon = 0;
if (step == 0) return 1;
if (step == 1) return 0;
if (step >= 2){
photon = pow(2,step) - NumberOfElectrons(step);
return photon;
}
}


double RefractiveIndex(double height, double x0){


// I would like to be able to create an array of refractive indices using the Sellmiere Equation, however a code I wrote to do this turned out to conflict with a declaration function in ROOT. It may have had too many arrays.
//For now I just input a constant refractive index.


//I input an arbitrary value for the refractive index at seas level and modify it as the muon emits Cherenkov Photons as it passes down through the atmosphere, eventually reaching this value



  double n =  (1.002) + ((TMath::Exp(-height/x0))/100);



cout << " Refractive Index: " << n <<endl;

return TMath::Abs(n);
}




double EnergyToBeta (double E, double& nSoll)
{


  //EnteR Muon mass values in GeV/c^2:

if (MUON==1) const double m = 0.10565837; 
if (MUON==0) const double m = 0.10565837; 


//relativistic pulse

double p=pow ((pow(E,2)-pow(m,2)),0.5);
double v=pow ((pow((m/p),2)+1), -0.5);

nSoll = 1./v;
return v;
}



//return the random number of photons generated at a given height and energy from the Frank-Tamm function:

double FrankTammFormula(const double height, const double energy, const double x0){
       

double nSoll = 0;
double beta = EnergyToBeta(energy ,nSoll);


//could insert a simple version of the sellmiere formula code here...


 double nIsAtHeight = RefractiveIndex(height, x0);
double dNdz = 0;


cout << " beta : " << beta << " nSoll: " << nSoll << endl;


if (beta *nIsAtHeight>1){

//Frank-Tamm Formula for energy lost by Muon in the form of Cherenkov Photons
//multiplied by the wavelength sensativity of the PMT detector corresponding to the maximum Quantum Efficieny, which is 20% Q.E at 320 nm for a PMT used in a VERITAS detector.


dNdz = ((2 *TMath::Pi())/137.)*(1-pow(beta*nIsAtHeight,-2))*pow(320e-9,-1);




}
else dNdz=0;

cout << "dNdz(number of photons/m) : " << dNdz << endl;

return dNdz;
}

double CherenkovAngle (const double height, const double energy, const double x0){
double nSoll = 0;
double beta = EnergyToBeta(energy,nSoll);
 double n = RefractiveIndex(height,x0);
if (beta*n>1){
              
double costheta = pow(beta*n,-1);
double thetaRad = TMath::ACos(costheta);
return thetaRad;
}
else return 0 ;
}


void CosmicRayMuonCherenkovSimulation(){
double Ezero = 1e3; 
double hstart = 80.;
double height = hstart;

int choice = 1;

double aperture = 0;


double detectorheight = 0;


 double incidentangle =0;

 double focal =0 ;

 double pixelradius =0;

 double PMTpixels =0;

 double QuantumEff =0;

 double Distance =0;

//Terminal Interface...


cout<< " _____________________________________________ "<< endl;

cout<< " __Cherenkov Emission from Atmospheric Muons__ "<< endl;
cout<< " _____________________________________________ "<< endl;
cout<< " ________John_Campbell___Summer_2011__________ "<< endl;
cout<< " _____________________________________________ "<< endl;


 cout<<" The first set of simulation routines should take"<<endl;
 cout<<" about 10-15 seconds to complete, after this wait"<<endl;
 cout<<" untill the detector simulation finsihes(may take about 2 to 3 mins)"<<endl;

choice=2;

MUON = 0;



cout<< "Enter Muon Kinetic Energy in GeV(e.g 8 GeV): " << endl;
cin>>Ezero;
cout<< "Enter height where shower starts in km(e.g 8 KM): " << endl;
cin>>height;


cout<<"Enter Incident Angle"<< endl;
 cin>>incidentangle;



cout<< "Enter height of the detector above seas level in km(1 kilometer for VERITAS): " << endl;
cin>>detectorheight;

 cout<<"Enter the Distance of the detector from the shower core(in meters) to calculate the impact parameter"<< endl;
 cin>>Distance;


cout<<"Enter Radius of Detector aperture in meters(12m for VERITAS)"<< endl;
 cin>>aperture;

 cout<<"Enter the focal point of the detector in meters(12m for VERITAS)"<< endl;
 cin>>focal;


 cout<<"Enter the radius of an individual PMT used in the detector array in mm(29mm for VERITAS PMTs)"<< endl; 
 cin>>pixelradius;

 cout<<"Enter the Quantum Efficiency of an individual PMT in % (~20% for VERTIAS PMTs)"<< endl;
 cin>>QuantumEff;

 cout<<"Enter the number of PMT pixels in the camera(499 for VERITAS)"<< endl; 
 cin>>PMTpixels;

cout<< " _________________________________ "<< endl;




// timing...

TStopwatch processtime;
processtime.Start();


//declarations


double heightEval = 0;
double heightStart = 0;
double heightstep = 0;
int SumElectrons = 0;
int SumPhotons = 0;
int Nlength = MaxRadiation (Ezero);
double hstep = 0.01;
int ChPhotons = 0;
double EnergyFraction = 0;
double Split =0;
double dNdz = 0;


int weight = 0;
int NumberOfMCPhotons = 0;
int NumberOfCherenkovPhotons = 0;
double EnergyForCherenkovPhotons = 0;
double thetaCh = 0;
double RadiusCh = 0;
double xpos,ypos,phi;
double CurrentHeight = 0;


// Generate ROOT tree Ntuples of main data



TNtuple *NTcherenkov = new TNtuple( "NTcherenkov " , " tree containing cherenkov production data " , " height:dNdz:Photons:thetaCh:RadiusCh " );
TNtuple *NTChPhoton = new TNtuple( "NTChPhoton " , " tree containing single photon information " , " height:RadiusCh:thetaCh:xpos:ypos " ) ;

// Store data as histograms and Graphs.

TH1D* h1ChPhotonProduction = new TH1D( " h1ChPhotonProduction " , " Created Cherenkov Photons " ,30,0,30);

TH1D* h1ChRadius = new TH1D( " h1ChRadius " , " Radius of Cherenkov photons " ,30,0,160);


TH2D*h2RadiusLength = new TH2D("Cherenkov Radius over length", "Radius of Cherenkov Photons per unit length",30,0,160,300,-150,150);


TH1D* h1ChPhi = new TH1D( " h1ChPhi " , " Free angle of Cherenkov Photons " ,200,0,360);


TH2D* h2ChPhiposition = new TH2D("h1ChPhiposition", "Free angle of Chrenkov Photons over length x",200,0,360,300,-150,150);



TH1D* h1ChTheta = new TH1D( " Cherenkov angle Theta " , "Angle of Cherenkov photons " ,300,0,1);


TH2D*h2AngDistr = new TH2D("Cherenkov angle Theta vs length","Angular Distribution of Cherenkov photons per unit length ",30,0,1,300,-150,150); 



TH3D* h3AngDistrArea = new TH3D( "Cherenkov angle Theta vs Area " , " Angular Distribution of Cherenkov photons per Area " ,30,0,1,300,-150,150,30,0,height);




TH2D* h2DistrPhotons = new TH2D( " Photon Distribution per unit area " , " Photons on the ground ",300,-150,150,300,-150,150);

TH3D* h3CherenkovCone = new TH3D( "Cherenkov Cone " , " Cherenkov Photons along Path " ,300,-150,150,300,-150,150,30,0,height);


 TH3D* h3CherenkovRadius= new TH3D( "Cherenkov ConE Radius " , " Cherenkov Photons along cone " ,300,-150,150,300,-150,150,30,0,height);


TH3D* h3RadiusArea = new TH3D("Cherenkov Cone over Area", "Cherenkov 3D Cone",30,0,160,300,-150,150,300,-150,150);



TH3D* MuonPath = new TH3D("Muon path xpos", "Muon path 2D",30,0,height,300,-150,150,300,-150,150);



TH2D* h2RadiusArea = new TH2D("Cherenkov Cone over Area", "Cherenkov 3D Cone",30,0,160,300,-150,150);



TH2D* h2MuonPath = new TH2D("Muon path xpos", "Cherenkov Ring in aperture",30,0,height,300,-150,150);

 TH3D* h3MuonPath = new TH3D("Muon path xpos", "Muon path 2D",30,0,height,300,-150,150,300,-150,150);



double BIN[4];
BIN[0] = 0.5;
BIN[1] = 1.5;
BIN[2] = 2.5;
BIN[3] = 3.5;


cout<< "We perform " << Nlength << " steps. " << endl;


//main loop

if (choice == 1){
for (unsigned int i = 0; i<=Nlength; i++){
heightStart = height;
height = height - x0(height);
heightEval = heightStart-height;

cout<< "we start at: " << heightStart << " and go to " << height << endl;
cout<< " Step: " << i << " Height: " << height << " refract. index: " <<
RefractiveIndex(height)<<endl;



EnergyFraction = Ezero/pow(2,i);
h1EnergyFraction->Fill(EnergyFraction);
SumPhotons = SumPhotons + NumberOfPhotons (i);
SumElectrons = SumElectrons + NumberOfElectrons(i);
NTems->Fill(height,NumberOfElectrons(i),NumberOfPhotons(i),SumElectrons,
SumPhotons,EnergyFraction,x0(height));


heightstep = 0; 
if (i ==Nlength-1) hstep = 0.2;
if (i >0){
while (heightEval>=heightstep){
heightstep = heightstep + hstep;


EnergyForCherenkovPhotons = EnergySplitting(EnergyFraction);
h1AP->Fill(AP);
h1EnergyFraction->Fill(EnergyForCherenkovPhotons/EnergyFraction);
dNdz = FrankTammFormula(heightStart-heightstep,EnergyForCherenkovPhotons);
thetaCh = CherenkovAngle(heightStart-heightstep,EnergyForCherenkovPhotons);
RadiusCh = TMath::Tan(thetaCh )*(heightStart-heightstep)*1e3; // t o g e t i t i n m.
NTcherenkov->Fill(heightStart-heightstep,dNdz,NumberOfCherenkovPhotons,DEG(thetaCh),RadiusCh);
weight = NumberOfElectrons(i);



//Looping over all steps...


for (int k = 0; k<weight; k++){


//Monte-Carlo of dNdz to generate an int number of photons from the Frank-Tamm formula
// this is a crucial step and was one which was overlooked in an 
//attempt at a matlab simulation of this process, which caused it to fail.


NumberOfMCPhotons = PoissonRandomNumber(dNdz);
for (int j = 0; j<NumberOfMCPhotons; j++){


// Monte-Carlo of Phi angles
phi = gRandom->Uniform(0,2*TMath::Pi());
xpos = RadiusCh*TMath::Cos(phi)+ RadiusCh*TMath::Sin(incidentangle);
ypos = RadiusCh*TMath::Sin(phi)+ RadiusCh*TMath::Sin(incidentangle);

//ntuples


NTChPhoton->Fill(height,RadiusCh,DEG(thetaCh),xpos,ypos);


//histograms

h1ChPhi->Fill(DEG(phi));
h2DistrPhotons->Fill(xpos,ypos);
}
NumberOfCherenkovPhotons=NumberOfCherenkovPhotons+NumberOfMCPhotons;
h1ChTheta->Fill(DEG(thetaCh),NumberOfMCPhotons);
h1ChRadius->Fill(RadiusCh,NumberOfMCPhotons);
h1ChPhotonProduction->Fill(double(NumberOfMCPhotons));
}
}
}
}
cout<< "Total photons: " << SumPhotons << " Conversion electrons: " << SumElectrons << " Cherenkov Photons: " << NumberOfCherenkovPhotons<< endl;

}
else if (choice == 2)
{ 

//simulate a single muon

hstep = height/1000.;
CurrentHeight = height;
for(int i = 0; i<1000; i++){
CurrentHeight = CurrentHeight - hstep;



 dNdz = FrankTammFormula(CurrentHeight,Ezero, height);



 thetaCh = CherenkovAngle(CurrentHeight,Ezero, height) + incidentangle;



RadiusCh = TMath::Tan(thetaCh-incidentangle)*CurrentHeight*1e3; 

NTcherenkov->Fill(CurrentHeight,dNdz,NumberOfCherenkovPhotons,DEG(thetaCh),RadiusCh);


cout << " Cherenkov Angle : " << 1+thetaCh << " Cone Radius: " << RadiusCh << endl;




NumberOfMCPhotons = PoissonRandomNumber (dNdz); 
for (int j = 0; j<NumberOfMCPhotons; j++){


  //Monte-Carlo of Phi angles

phi = gRandom->Uniform(0,2*TMath::Pi());
xpos = RadiusCh*TMath::Cos(phi)+ RadiusCh*TMath::Sin(incidentangle);
ypos = RadiusCh*TMath::Sin(phi)+ RadiusCh*TMath::Sin(incidentangle);


//ntuples


NTChPhoton->Fill(CurrentHeight,RadiusCh,DEG(thetaCh),xpos,ypos);


//histograms

h1ChPhi->Fill(DEG(phi));

h2DistrPhotons->Fill(xpos,ypos);

h3CherenkovCone->Fill(height,xpos,ypos);

h2RadiusLength->Fill(RadiusCh,xpos);



h3RadiusArea->Fill(RadiusCh,xpos,ypos);


MuonPath->Fill(CurrentHeight,xpos,ypos);


h3CherenkovRadius->Fill(RadiusCh,xpos,ypos);

h2ChPhiposition->Fill(DEG(phi),xpos);


h3AngDistrArea->Fill(DEG(thetaCh),xpos,ypos);
 
h2AngDistr->Fill(DEG(thetaCh),xpos);



//Graphs 


}
NumberOfCherenkovPhotons=NumberOfCherenkovPhotons+NumberOfMCPhotons;
h1ChTheta->Fill(DEG(thetaCh),NumberOfMCPhotons);
h1ChRadius->Fill(RadiusCh,NumberOfMCPhotons);

h1ChPhotonProduction->Fill(double(NumberOfMCPhotons));
}
}




//Simulate High Energy Coulomb Scattering using a form of the Klein-Nishina formula...
//I do this to compare the High Energy Coulomb scattering simulated in the same way as the bremstrahlung photons are simulated using a vertex function.


//this is just to make sense of some of the numbers involved in the detected number of photons.


int i;
  int Emin,Estep,Egamma;
  double Emax,Emuon,PI,DSigma_DE,R,Bohr_Radius_sqrd,Eo_muon,s,psi;

  PI=3.1415927;

  Bohr_Radius_sqrd=2.817e-13*2.817e-13/1e-24;/*barns*/
   




  int Emin = 83000; //Critical Energy

  int Egamma = 80000000000;  //Typical Energy

  int Estep = height*1000;  //Height Steps



TH2D*KleinNishina = new TH2D("Differential cross section","Klein-Nishina DiffXsect DSigma/DE",300,-50,50,300,-50,50);




  Eo_muon= 105658370;   /* rest mass of muon in eV */

  psi=Egamma/Eo_muon;

  Emax=Egamma*2*psi/(1+2*psi);  /* compton edge */


  fprintf(stderr,"Calculating Differential cross section (in barns/eV) for \n\t %d < E < %g in steps of %d and %d eV incident Photons\n",Emin,Emax,Estep,Egamma);


 Emuon=Emin-Estep;

 for(i=Emin;i<=Emax; i+=Estep)
   {

     Emuon+=Estep;

    

     s=Emuon/Egamma;

     DSigma_DE=s*(s-2/psi)/(1-s);

     DSigma_DE=DSigma_DE+s*s/(1-s)/(1-s)/psi/psi;

     DSigma_DE=2.0+DSigma_DE;

     DSigma_DE=PI*Bohr_Radius_sqrd*DSigma_DE/Eo_muon/psi/psi;
     

     KleinNishina->Fill(-Egamma/Emuon,DSigma_DE);



   }





 cout<<"Simulating PMT detector response..."<<endl;
 cout<<"Please Wait........................"<<endl;





//Part 2: Photon Collection, Reflection and PMT Quantum Efficiencey
//Begin Photon collection simulation:





//Placing circle to collect photons: 

	    
	    
	    unsigned int hits = 0;
    
    for (unsigned int z = 0; z < NumberOfCherenkovPhotons; ++z);
    
		if (xpos * xpos + ypos * ypos < aperture)
	    ++hits;
    
    	   
	 
		const Float_t Events = hits *dNdz/pixelradius*0.1;


  const int     BinNumbers  =       PMTpixels; // corresponds to the number of PMTs in a Veritas Detector.

  const Float_t Production  =   NumberOfCherenkovPhotons; // Number of generated photon per event
  
  const Float_t DecayLenght =   MaxRadiation; // Maximum Radiation depth of Cherenkov photons along the mfp of the muon.



  Double_t DetectorHalfWidth =    1000.0 ; // mm
  Double_t DetectorYposition =    1000.0 ; // mm
  
  // I use TRandom3
  TRandom3 *myRandom3 = new TRandom3();
  myRandom3->SetSeed(0);
 
  gROOT->SetStyle("Plain"); 
  TStyle *myStyle = gROOT->GetStyle("Plain");


  myStyle->SetPalette(1);


  //Enter Custom Palette Here....

  myStyle->cd();

  TCanvas *PMTView = new TCanvas("c1","c1");

  c1->Divide(1,2);
  // Create a tree
  
  TTree *tree = new TTree("T"," ROOT tree");


  TH2D*histo = new TH2D("Reflected Cherenkov Photons","Cherenkov Rays",BinNumbers,-1000,1000,BinNumbers,-1000,1000);
  

  // TH2*h4 = new TH2("histo","histo",BinNumbers,-1000,1000,BinNumbers,-1000,1000);


  //TH2I *hc = new TH2I("histo","histo",BinNumbers,-1000,1000,BinNumbers,-1000,1000);
  


  // Create Branches In order to plot the stopping distribution in the 3d space
  Float_t px,py,pz;  // position in which the photon stops
  Float_t dx,dy,dz;  // direction vector of the photon
  Float_t lx,ly,lz;  // 



  tree->Branch("px",&px,"px\F");
  tree->Branch("py",&py,"py\F");
  tree->Branch("pz",&pz,"pz\F");



  
  Double_t theta1,phi1,sigma1,rx,ry,rz,qx,qy,qz,zpos,r,zfr;
  Double_t u2,v2,theta2,phi2;
  
  for (int i=0;i<Events;i++){




    //generate incident angles... 

    theta1 = (90-thetaCh)+incidentangle;
    phi1   = (90-phi)+incidentangle;
    sigma1 = (incidentangle);
    

    //use ray tracing to find the z-position on a parabolic mirror:

    zpos=(xpos*xpos +ypos*ypos)/aperture ; 




    // stop point of the photon


    rx     = xpos*RadiusCh * TMath::Cos(thetaCh) * TMath::Sin(phi) ;
    ry     = ypos*RadiusCh * TMath::Sin(thetaCh) * TMath::Sin(phi) ; 
    rz     = zpos*RadiusCh * TMath::Cos(phi) ;
   

   

    //define full rotational symmetry and obtain position on the z-axiz as a function of R:


    r= sqrt(rx*rx + ry*ry + rz*rz);



    zfr=(r*r)/4*focal ; 




    //Pass Photons through a Rotation Matrix...



    qx =  ((rz*TMath::Cos(phi1)*TMath::Cos(thetaCh))-(TMath::Sin(phi1)*TMath::Sin(thetaCh))) ;


    qy =  ((rz*TMath::Cos(phi1)*TMath::Sin(thetaCh))-(TMath::Sin(phi1)*TMath::Cos(thetaCh))) ;


    qz =  (TMath::Cos(phi)*TMath::Cos(phi1)) ;



    px      = (rx * TMath::Cos(theta1)) + TMath::Sin(theta1) * qx ;


    py      = (ry * TMath::Cos(theta1)) + TMath::Sin(theta1) * qy ;


    pz      = (rz * TMath::Cos(theta1)) - TMath::Sin(theta1) * qz ;






    tree->Fill();




    
    // For each point generate a projection of theta and phi.
    
    for (int j=0;j<Production;j++){
      u2 = myRandom3->Uniform(0,1);
      v2 = myRandom3->Uniform(0,1);


      theta2 =     (2*TMath::Pi()* u2*aperture);  //Generate reflected angles over the aperture.
      phi2   = TMath::ACos(2*u2-1);
      
      // direction on the unitary sphere
      dx = TMath::Cos(theta2) *  TMath::Sin(phi2) ;
      dy = TMath::Sin(theta2) *  TMath::Sin(phi2) ; 
      dz = TMath::Cos(phi2) ;
      
     

      //Check if the photon arrive on the surface of the detector and construct the 
      //Cherenkov Ring from the muon track using the Hit-Postition Matching technique.
     


      ly=DetectorYposition-py;



      if (dy>0){       // if dy<Dtheta the photon goes in the other direction... 
	lx=ly/(dy/dx); // simple proportion - x axis
	lz=ly/(dy/dz); // simple proportion - z axis

	//cout << (px+lx) << " " << (pz+lz) << endl;

	if( ((px+lx) <  DetectorHalfWidth) &&
	    ((px+lx) > -DetectorHalfWidth) &&
	    ((pz+lz) <  DetectorHalfWidth) &&
	    ((pz+lz) > -DetectorHalfWidth)    )



	  
	  if (myRandom3->Uniform(0,1)<(QuantumEff/100)) histo->Fill(px+lx,pz+lz); // Enter Quantum Efficiency of PMTs
	  


h2RadiusArea->Fill(thetaCh,pz);



h2MuonPath->Fill(px+lx+(zfr*incidentangle),pz+lz+(zfr*incidentangle));




h3MuonPath->Fill(RadiusCh,px+lx,pz+lz);


      }
    }
  }
  

  //Show the Light yield...

  cout << "The light yield is " <<  100*histo->GetEntries()/(Production) << "%" << endl;
  
 

  PMTView->cd(1);
  tree->Draw("px:py:pz");
  
  PMTView->cd(2);
  histo->Draw("COLZ");
  histo->GetXaxis()->SetTitle("X [mm]");
  histo->GetYaxis()->SetTitle("Z [mm]");
  histo->SetTitle("Photomultiplier Surface");
  
  PMTView->SaveAs("PhotonData.pdf");  






  //10-Dynode PMT Simulation Routine.

  //this Routine simulates a camera of 499 PMT pixels from the response of a single PMT and performing a monte carlo disribution over the mean field of view to simulate the entire detector response.



gRandom->SetSeed(0);



TH1D *h1 = new TH1D("h1", "Number of electrons after 10 dynodes, dNdz Number of Cherenkov Photons per PMT ;n", 100, 0, 3000);
TH1D *h2 = new TH1D("h2", "Mean out distribution, 499 times Number of Cherenkov Photons", 20, 0.0e+03, 40.0e+03);
TH1D *h3 = new TH1D("h3", "Fraction of zeros, 499 times Number of Cherenkov Photons", 20, -1.0e+4, 1.0e+4);
TH1D *h4 = new TH1D("h4", "Variances, 499 times Number of Cherenkov Photons",20, -10.0e+06, 100.0e+06);


TCanvas *c1 = new TCanvas();

c1->Divide(2,2);

 int const npseudos = dNdz/pixelradius*0.1;   //Number of Cherenkov Photons per pixel 

 int const ntest = PMTpixels;  // number of pixels in the PMT camera

 float const mean = 3.5; // mean field of view for the PMT camera in degrees



//Generate a PMT with 10 dynodes (as used in a VERITAS camera) and generate a Poisson distribution of electrons at each dynode, adding on the previous distribution.


int stage1;
int stage2;
int stage3;
int stage4;
int stage5;
int stage6;
int stage7;
int stage8;
int stage9;
int stage10;

int final;


for(int tests=0;tests<ntest;tests++) {

if(tests % 100 == 0) cout << "Ntest: " << tests << endl;
int mgain[npseudos];

int nsum = 0;

for(int loop=0;loop<npseudos;loop++) {

final = 0;

stage1 = gRandom->Poisson(mean);

for(int i=0;i<stage1;i++) {

stage2 = gRandom->Poisson(mean);

for(int ii=0;ii<stage2;ii++) {

stage3 = gRandom->Poisson(mean);

for(int iii=0;iii<stage3;iii++) {

stage4 = gRandom->Poisson(mean);

for(int iv=0;iv<stage4;iv++) {

stage5 = gRandom->Poisson(mean);

for(int v=0;v<stage5;v++) {

stage6 = gRandom->Poisson(mean);

for(int v=0;v<stage6;v++) {

stage7 = gRandom->Poisson(mean);

for(int v=0;v<stage7;v++) {

stage8 = gRandom->Poisson(mean);

for(int v=0;v<stage8;v++) {

stage9 = gRandom->Poisson(mean);

for(int v=0;v<stage9;v++) {

stage10 = gRandom->Poisson(mean);



final += stage10;

}
}
}
}
}
}
}
}
}




//sum up the number of electrons after 10 dynodes.


if(tests == 0) h1->Fill(final);

mgain[loop] = final;

nsum += final;


}

//Generate Poisson distributions for the mean, varience, sum and fraction of zero values.

float xmean = 0.0;

float xvar = 0.0;

float xsum = 0.0;

int nzero = 0;

xmean = (float)nsum / (float)npseudos;

for(int ivar=0;ivar<npseudos;ivar++) {

xsum += ((float)mgain[ivar] - xmean)*((float)mgain

[ivar] - xmean);

if(mgain[ivar] <= 0) nzero++;

}
xvar = xsum/((float)npseudos - 1.0);

float xzero = (float)nzero / (float)npseudos;




h2->Fill(xmean);
h3->Fill(xzero);
h4->Fill(xvar);


if(tests == 0) {


cout << "Mean: " << xmean << endl;
cout << "Variance: " << xvar << endl;
cout << "Frac zeros: " << xzero << endl;


c1->cd(1);
h1->Draw();

 }
 

c1->Update();

c1->cd(2);
h2->Draw();

c1->cd(3);
h3->Draw();

c1->cd(4);
h4->Draw();

c1->cd();


}


 //Part 3 Cherenkov Ray tracing routine...
//this takes the output distribution of photons detected and places them in their corresponding bins in a detector image


  cout<<"Performing Ray-Tracing Simulation of Cherenkov Photons"<<endl;



int choice2 = 1;

double incidentangle =0 ;


//The Simulation ray tracing technique currently allows the Coulomb scattering part of the code to be viewed as a function    of impact parameter and the bremsstrahlung part of the vertex function to be viewed as a function of incident angle.         combining the 2 together generates strange numbers and is limited in the values that can be input to generate images.





//Uncomment the section below to switch on the vertex function for Bremsstrahlung and High Energy Coulomb Scattering.


AP2=2;  
choice2 = 1;




TStopwatch processtime;
processtime.Start();

//declaring variables used in generated shower...
double heightEval = 0;
double heightStart = 0;
double heightstep = 0;
int SumElectrons = 0;
int SumPhotons = 0;
int Nlength = MaxRadiation (Ezero);
double hstep = 0.01;
int ChPhotons = 0;
double EnergyFraction = 0;
double Split =0;
double dNdz = 0;


int weight = 0;
int NumberOfMCPhotons = 0;
int NumberOfCherenkovPhotons = 0;
double EnergyForCherenkovPhotons = 0;
double thetaCh = 0;
double RadiusCh = 0;
double xpos2,ypos2,phi;
double CurrentHeight = 0;



TNtuple *NTems = new TNtuple("NTems","tree containing the EMS in formation " , "height:Electrons:Photons:SumElectrons:SumPhotons:EnergyFraction:x0 " );
TNtuple *NTcherenkov = new TNtuple( "NTcherenkov " , " tree containing cherenkov production data " , " height:dNdz:Photons:thetaCh:RadiusCh " );
TNtuple *NTChPhoton = new TNtuple( "NTChPhoton " , " tree containing single photon information " , " height:RadiusCh:thetaCh:xpos:ypos " ) ;




TH2D* DetectorPhotons = new TH2D( " h2DistrPhotons " , " Photons in the Detector ",499,-height*1000,height*1000,499,-height*1000,height*1000);

double BIN[4];
BIN[0] = 0.5;
BIN[1] = 1.5;
BIN[2] = 2.5;
BIN[3] = 3.5;


TH1D* h1AP2 = new TH1D( "h1AP" , "AP switch (cross check) " ,3,BIN);


cout<< "performing " << Nlength << " steps. " << endl;


if (choice2 == 1){
for (unsigned int d = 0; d<=Nlength; d++){



heightStart = height;


height = focal - x0(height);


heightEval = heightStart-height;



EnergyFraction = Ezero/pow(2,d);

SumPhotons = SumPhotons + NumberOfPhotons (d);



SumElectrons = SumElectrons + NumberOfElectrons(d);


NTems->Fill(height,NumberOfElectrons(d),NumberOfPhotons(d),SumElectrons,


SumPhotons,EnergyFraction,x0(height));



heightstep = 0;
if (d ==Nlength-1) hstep = 0.2;
if (d >0){
while (heightEval>=heightstep){
heightstep = heightstep + hstep;



EnergyForCherenkovPhotons = EnergySplitting(EnergyFraction);
h1AP2->Fill(AP2);

 dNdz = FrankTammFormula(heightStart-heightstep,EnergyForCherenkovPhotons, height);


 thetaCh = CherenkovAngle(heightStart-heightstep,EnergyForCherenkovPhotons, height);


 RadiusCh = TMath::Tan(thetaCh)*(heightStart-heightstep)*1e3; // t o g e t i t i n m.


 NTcherenkov->Fill(heightStart-heightstep,dNdz,NumberOfCherenkovPhotons,DEG(thetaCh),RadiusCh);



weight = NumberOfElectrons(i);


for (int k = 0; k<weight; k++){

NumberOfMCPhotons = PoissonRandomNumber(dNdz); 
for (int j = 0; j<NumberOfMCPhotons; j++){

phi = gRandom->Uniform(0,2*TMath::Pi());



 xpos2 = RadiusCh*TMath::Cos(phi) + (sigma1/10)*(RadiusCh*TMath::Cos(phi)+ RadiusCh*TMath::Cos(thetaCh))

   + Distance*(RadiusCh*TMath::Cos(theta2) + RadiusCh*TMath::Cos(phi2));
 

 
 
 ypos2 = RadiusCh*TMath::Sin(phi) ;



NTChPhoton->Fill(height,RadiusCh,DEG(thetaCh),xpos2,ypos2);

DetectorPhotons->Fill(xpos2,ypos2);
}
NumberOfCherenkovPhotons=NumberOfCherenkovPhotons+NumberOfMCPhotons;


}
}
}
}
cout<< "Total photons: " << SumPhotons << " total electrons: " << SumElectrons << " Cherenkov Photons: " << NumberOfCherenkovPhotons<< endl;

}



//Generate ROOT files of the ntuples,trees, ect for portable analysis and viewing of the graphs


TFile *f12 = new TFile("Detector_data.root", "RECREATE");



NTems->Write();
NTcherenkov->Write();
NTChPhoton->Write();

DetectorPhotons->Write();

h1AP2->Write();
f12->Write();
processtime.Stop();
cout<<"Ray-Tracing Monte Carlo Complete. "<<endl;
processtime.Print();




TFile *f1 = new TFile("Cherenkov_data.root", "RECREATE");


NTcherenkov->Write();
NTChPhoton->Write();
h1ChPhotonProduction->Write();
h1ChRadius->Write();
h1ChPhi->Write();
h1ChTheta->Write();
h2DistrPhotons->Write();
h3CherenkovCone->Write();
h3CherenkovRadius->Write();
h3AngDistrArea->Write();
h2AngDistr->Write();
h2RadiusLength->Write();
h2ChPhiposition->Write();


h3RadiusArea->Write();



MuonPath->Write();


f1->Write();



cout<<"Cherenkov Simulation Complete. "<<endl;


//Input Graphics routine here...




//Generating Graphics

  

  TCanvas *c2 = new TCanvas("c2","Cherenkov Photon Plots",200,10,700,900);

   //c1.Divide(2,3); //this is for dividing the canvas into 2 sets of 3 plots.



   c2->SetFillColor(42);
   gStyle->SetFrameFillColor(42);
   

   //Set Custom made Colour Palette(adds more contrast)

gStyle->SetPalette

void set_plot_style(){


    const Int_t NRGBs = 5;
    const Int_t NCont = 255;

    Double_t stops[NRGBs] = { 0.00, 0.34, 0.61, 0.84, 1.00 };
    Double_t red[NRGBs]   = { 0.00, 0.00, 0.87, 1.00, 0.51 };
    Double_t green[NRGBs] = { 0.00, 0.81, 1.00, 0.20, 0.00 };
    Double_t blue[NRGBs]  = { 0.51, 1.00, 0.12, 0.00, 0.00 };
    TColor::CreateGradientColorTable(NRGBs, stops, red, green, blue, NCont);
    gStyle->SetNumberContours(NCont);
   }




   
   title = new TPaveText(.2,0.96,.8,.995);
   title->SetFillColor(33);
   title->AddText("Cherenkov Photon Distributions on the Ground");
   title->Draw();

   

   pad1 = new TPad("pad1","Color mesh",0.03,0.50,0.98,0.95,21);

   pad2 = new TPad("pad2","Color mesh",0.03,0.02,0.98,0.48,21);
   

   pad1->Draw();

   pad2->Draw();
   //
   



   h2DistrPhotons->SetContour(48);
   h2DistrPhotons->SetFillColor(45);


// Draw the data in pad1 with color mesh option
   pad1->cd();
  
  
   h2DistrPhotons->Draw("surf1");
   

   // Draw the data in pad2 with colz option
   pad2->cd();
  
  
   h2DistrPhotons->Draw("colz");


   


  //Finishing graphics for 1st canvas...





//Generating second canvas for the 3D plots: 




 TCanvas *c3 = new TCanvas("c3","3D Surfaces Plots",200,10,700,900);

   //c1.Divide(2,3); //this is for dividing the canvas into 2 sets of 3 plots.



   c3->SetFillColor(42);
   gStyle->SetFrameFillColor(42);
   

   //Set Custom made Colour Palette

gStyle->SetPalette

void set_plot_style(){


    const Int_t NRGBs = 5;
    const Int_t NCont = 255;

    Double_t stops[NRGBs] = { 0.00, 0.34, 0.61, 0.84, 1.00 };
    Double_t red[NRGBs]   = { 0.00, 0.00, 0.87, 1.00, 0.51 };
    Double_t green[NRGBs] = { 0.00, 0.81, 1.00, 0.20, 0.00 };
    Double_t blue[NRGBs]  = { 0.51, 1.00, 0.12, 0.00, 0.00 };
    TColor::CreateGradientColorTable(NRGBs, stops, red, green, blue, NCont);
    gStyle->SetNumberContours(NCont);
   }




   
   title = new TPaveText(.2,0.96,.8,.995);
   title->SetFillColor(33);
   title->AddText("3D Plots of Cherenkov Cones");
   title->Draw();

   

   pad1 = new TPad("pad1","Color mesh",0.03,0.50,0.98,0.95,21);

   pad2 = new TPad("pad2","Color mesh",0.03,0.02,0.98,0.48,21);
   

   pad1->Draw();

   pad2->Draw();
   //
   


   MuonPath->SetContour(48);
   MuonPath->SetFillColor(45);


   h3RadiusArea->SetContour(48);
   h3RadiusArea->SetFillColor(45);


// Draw the data in pad1 with Lego2 option
   pad1->cd();
  
  
   h3RadiusArea->SetMarkerColor(kBlue);

   h3RadiusArea->Draw("Lego1");
   

// Draw the data in pad2 with Lego2 option
   pad2->cd();
  
   
   MuonPath->SetMarkerColor(kBlue);

   MuonPath->Draw("Lego1");





   //generate a Fouth Canvas(for testing purposes): 








 TCanvas *c4 = new TCanvas("c4","Test Plots",200,10,700,900);

   //c1.Divide(2,3); //this is for dividing the canvas into 2 sets of 3 plots.



   c4->SetFillColor(42);
   gStyle->SetFrameFillColor(42);
   

   //Set Custom made Colour Palette

gStyle->SetPalette

void set_plot_style(){


    const Int_t NRGBs = 5;
    const Int_t NCont = 255;

    Double_t stops[NRGBs] = { 0.00, 0.34, 0.61, 0.84, 1.00 };
    Double_t red[NRGBs]   = { 0.00, 0.00, 0.87, 1.00, 0.51 };
    Double_t green[NRGBs] = { 0.00, 0.81, 1.00, 0.20, 0.00 };
    Double_t blue[NRGBs]  = { 0.51, 1.00, 0.12, 0.00, 0.00 };
    TColor::CreateGradientColorTable(NRGBs, stops, red, green, blue, NCont);
    gStyle->SetNumberContours(NCont);
   }




   
   title = new TPaveText(.2,0.96,.8,.995);
   title->SetFillColor(33);
   title->AddText("Test Plots");
   title->Draw();

   

   pad1 = new TPad("pad1","Color mesh",0.03,0.50,0.98,0.95,21);

   pad2 = new TPad("pad2","Color mesh",0.03,0.02,0.98,0.48,21);
   

   pad1->Draw();

   pad2->Draw();
   //
   



   h2MuonPath->SetContour(48);
   h2MuonPath->SetFillColor(45);


   h3MuonPath->SetContour(48);
   h3MuonPath->SetFillColor(45);


// Draw the data in pad1 with Surf2 option
   pad1->cd();
  
  
  KleinNishina->Draw("surf2");

   
   

// Draw the data in pad2 with Lego2 option
   pad2->cd();
  
   
   

   h2MuonPath->Draw("pollego2z");




//This ends the 3D plot section...





//Input Graphics routine here...
//Generating Graphics

  


TCanvas *c5 = new TCanvas("c5","Cherenkov Photon Plots",200,10,700,900);

   //c1.Divide(2,3); //this is for dividing the canvas into 2 sets of 3 plots.



   c5->SetFillColor(42);
   gStyle->SetFrameFillColor(42);
   

   //Set Custom made Colour Palette

gStyle->SetPalette

void set_plot_style(){


    const Int_t NRGBs = 5;
    const Int_t NCont = 255;

    Double_t stops[NRGBs] = { 0.00, 0.34, 0.61, 0.84, 1.00 };
    Double_t red[NRGBs]   = { 0.00, 0.00, 0.87, 1.00, 0.51 };
    Double_t green[NRGBs] = { 0.00, 0.81, 1.00, 0.20, 0.00 };
    Double_t blue[NRGBs]  = { 0.51, 1.00, 0.12, 0.00, 0.00 };
    TColor::CreateGradientColorTable(NRGBs, stops, red, green, blue, NCont);
    gStyle->SetNumberContours(NCont);
   }




   title = new TPaveText(.2,0.96,.8,.995);
   title->SetFillColor(33);
   title->AddText("Photon Distributions in the Detector Array");
   title->Draw();

   

   pad1 = new TPad("pad1","Color mesh",0.03,0.50,0.98,0.95,21);

   pad2 = new TPad("pad2","Color mesh",0.03,0.02,0.98,0.48,21);
   

   pad1->Draw();

   pad2->Draw();
   //
   

   DetectorPhotons->SetContour(48);
   DetectorPhotons->SetFillColor(45);


// Draw the data in pad1 with color mesh option
   pad1->cd();
  
  
   DetectorPhotons->Draw("surf1");
   

   // Draw the data in pad2 with colz option
   pad2->cd();
  
  
   DetectorPhotons->Draw("colz");



//This ends the 3D plot section...





//End Graphics Routine here...
//Program end...




}



