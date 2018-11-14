#define zdecay_cxx
#include "zdecay.h"
#include <TH1.h>
#include <TStyle.h>
#include <TCanvas.h>
#include <TLorentzVector.h>
#include <TFormula.h>
#include <TF1.h>
#include <TList.h>
#include "RooRealVar.h"
#include "RooDataSet.h"
#include "RooDataHist.h"
#include "RooGaussian.h"
#include "RooExponential.h"
#include "RooFitResult.h"
#include "TCanvas.h"
#include "RooPlot.h"
#include "TAxis.h"
#include <RooCBShape.h>
using namespace RooFit ;



class particles{
public:

   //*Constructor
   particles(double elePhi, double eleEta, double elePt, double eleCharge){
      phi = elePhi;
      eta = eleEta;
      Pt = elePt;
      charge = eleCharge;
   }

   //* destructor
   ~particles(){
    }

    //* Accessors
   double getphi(){
      return phi;
   }
   double geteta(){
      return eta;
   }
   double getPt(){
      return Pt;
   }
   double getE(){
      return Pt/(sin(2.0*atan(exp(-eta))));
   }
   double getPx(){
      return Pt*cos(phi);
   }
   double getPy(){
      return Pt*sin(phi);
   }
   double getPz(){
      return Pt*cos(phi)/sin(phi);
   }
   double getcharge(){
      return charge;
   }


   //* Mutators
   void setphi(double Phi){
      phi = Phi;
   }
   void seteta(double Eta){
      eta = Eta;
   }
   void setPt(double Pt){
      Pt = Pt;
   }
   void setcharge(double Charge){
      charge = Charge;
   }

private:
   double E,phi,theta,eta,Pt,Px,Py,Pz,charge;
};


//* This is the comparing function. It takes two values and returns False/True based on the expression used.
bool compare_trans(particles &particle1, particles &particle2){
    return particle1.getPt() > particle2.getPt();
   }


//* This is the main loop.
void zdecay::Loop()
{

   
   // cout << "Begin" << endl;
   if (fChain == 0) return;  
   
   TH1F *h1 = new TH1F("h1","Histogram;Mass;Count",150,0,600);
   // h1.FillRandom("gaus",250000);
   
   Long64_t nentries = fChain->GetEntriesFast();
   cout<<"Total entries "<<nentries<<endl<<endl;
   Long64_t nbytes = 0, nb = 0;

   for (Long64_t jentry=0; jentry<nentries;jentry++) {
      Long64_t ientry = LoadTree(jentry);
      //if (ientry < 0) break;
      nb = fChain->GetEntry(jentry);
      nbytes += nb;
      
      //* Creating a vector of objects of class particles.
      vector<particles> electrons;
  
      
            // Acquisition loop. If the event has at least one electron (or particle) then it makes an electron object
            // which is inserted into the vector - electrons

      if(nEle>1) {
         for (int i=0 ;i<nEle; i++){
            particles electron(elePhi -> at(i),eleEta->at(i),elePt->at(i),eleCharge->at(i));
            electrons.push_back(electron);
         }
         std::sort(electrons.begin(),electrons.end(), compare_trans);
         
         // Contructing the four vector for the first two particles with highest Pt - Summing them in V3 for display.
         TLorentzVector v1;
         v1.SetXYZT(electrons[0].getPx(),electrons[0].getPy(),electrons[0].getPz(),electrons[0].getE());
         TLorentzVector v2;
         v2.SetXYZT(electrons[1].getPx(),electrons[1].getPy(),electrons[1].getPz(),electrons[1].getE());
         TLorentzVector v3 = v1+v2;
         h1 -> Fill(v3.M());
      
               }
      
                                                    }

      // for ( int i = 14; i<=28; i++){
      //    double a = h1->GetBinCenter(i);
      //    h1 -> SetBinContent(i,100.);
      // }
      

   RooRealVar x("x","mass",0.,1000.);
   // x.setRange("Range1",35.75,50.13) ;
   x.setRange("Range2",110.101,200.) ;
   RooRealVar lambda("lambda", "slope", -0.001, -1., 1.);
   RooExponential expo("expo", "expo", x, lambda);
   RooExponential expo1("expo1", "expo", x, lambda);
   RooDataHist dh("dh","e-e+ peak histo",x,Import(*h1)) ;
   // RooFitResult* r = 
   RooExtendPdf background("background","Side Bands",expo,x,"FULL") ;
   RooFitResult* r = expo.fitTo(dh,Range("Range2"), Save(kTRUE)) ;
   r -> Print();
   
   RooPlot* frame = x.frame() ;
   RooPlot* frame1 = x.frame() ;

      
   dh.plotOn(frame,DrawOption("B"),DataError(RooAbsData::None),FillColor(kGray));
   expo.plotOn(frame);
   // background.plotOn(frame,LineStyle(kDashed),LineColor(kRed));
   

   RooDataSet *h2 = expo.generate(x,2500);
   
   
   h2-> plotOn(frame1, LineColor(kBlue));
   
   TH1* h2_hist = h2->createHistogram("h2_hist", x);
   RooDataHist dh2("dh2","e-e+ peak histo",x,Import(*h2_hist)) ;
   RooFitResult* r2 = expo1.fitTo(dh2, Save(kTRUE)) ;
   r2 -> Print();
   
   RooPlot *frame3 = x.frame();
   dh2.plotOn(frame3);
   expo1.plotOn(frame3);
   
   TCanvas* c = new TCanvas("canvas","canvas",800,800) ;
   c -> Divide(2,2);
   c -> cd(1);
   frame -> Draw();
   c -> cd(2);
   frame1 -> Draw();
   c -> cd(3);
   frame3 -> Draw();
}









      // Old fit functions and saving routine

      //! The crystalbal+expo fit routine
      // TF1 *fit = new TF1("fit", "crystalball(0) + expo(5)", 0., 1000.);
      // fit ->SetParameters(160.0,90.0,9.0,1.0,220.0, 0.01);
      // for (Int_t i = 0; i < fit->GetNpar(); i++) std::cout << i << " : " << fit->GetParName(i) << std::endl;
      // h1->Fit("fit", "", "", 18., 200.);
      

      // TF1 *fline = new TF1("fline","expo(0)",18.,550.);
      // fline->SetParameters(220.,-0.01);
   
      // h1->Fit("fline","","",18.,65.);
      // h1->Fit("fline","","",105.,600.);
      

      //! This part is to create a file and save the plots.

      // TFile f("Histo.root","recreate");
      // h1 -> Write();
      // fit -> Write();
      // f.Close();
      // c1 -> cd();
      // h1 ->Draw();
      

      
      
     
   
