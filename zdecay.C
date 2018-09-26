#define zdecay_cxx
#include "zdecay.h"
#include <TH1.h>
#include <TStyle.h>
#include <TCanvas.h>
#include <TLorentzVector.h>
#include <TFormula.h>




// using namespace std;

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

   //* Creating a vector of objects of class particles.
   // TCanvas c1("c1","Histogram here",600,800);
   cout << "Begin" << endl;
   // int verbose_level = 1;
   // int need = 0;

   if (fChain == 0) return;  
   // TCanvas c1("c1","c1",800,1000);
   TCanvas *c1 = new TCanvas("c1","Plots",200,10,900,900);

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
      // cout << "processing event "<<jentry <<endl;
      vector<particles> electrons;
      // cout<<"event no "<<jentry<<" eta of first electron "<<nEle<<endl;
      
      
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
      // if (v3.M()>60.0 & v3.M() <140.0)
      h1 -> Fill(v3.M());
      // if (v3.M()>60.0 & v3.M() <140){
      //    need += 1;
      // }
      // //if (verbose_level) cout << v3.M() <<endl;
      // if(verbose_level)cout << v3.M() <<endl;
      //delete electrons;
      //cout << "Pt,Phi,Eta for the top 2 highest Pt:"<<endl<<endl;
      // if (v3.M()>0) cout <<v3.M()<<jentry<<endl;
      // h1.Fill(v3.M(),1);
      // cout << v3.M()<<endl;
      // h1.Fill(v3.M());
      // if (verbose_level > 1){
      //    for(int i=0;i<2;i++){   
      //       cout << electrons[i].getPt()<<" "<< electrons[i].getphi()<<" "<< electrons[i].geteta()<<" "<<electrons[i].getcharge()<<endl;
      //    }  
      //    cout << endl<< endl<<"Four vectors(E,Px,Py,Pz)"<<endl<<endl;
      //    cout << "   e-"<<"          "<<"e+"<<"        "<<"sum"<<endl;

      //    for (int i=0;i<4;i++){
      //    cout << v1[i] <<"    " << v2[i] <<"    "<<v3[i]<< endl;
      //    }
      // }
      }

      
      // /cout <<endl<<endl;
      // cout << v3.M()<< endl;
      // cout <<endl;
      
      } //end jentry loop
      // h1.Draw();

      TF1 *fit = new TF1("fit", "[0]*exp(-0.5*pow((x-[1])/[2],2)) + [3]*exp(-[4]*x)", 0., 1000.);
      fit -> SetNpx(10000);
      fit ->SetParNames("amp_gauss", "Mean_Gaus", "sigma_Gaus", "amp_exp", "exp_arg");
      fit ->SetParameters(110.0, 90.0, 10.0, 220.0, 0.01);
      // fit.Draw();

      h1->Fit("fit", "", "", 18., 150.);

      TFile f("Histo.root","recreate");
      h1 -> Write();
      fit -> Write();
      f.Close();
      c1 -> cd();
      h1 ->Draw();
      TPad *newpad=new TPad("newpad","a transparent pad",0,0,1,1);
      newpad -> SetFillStyle(4000);
      newpad ->Draw();
      newpad ->cd();
      TPaveLabel *title = new TPaveLabel(0.1,0.5,0.3,0.3,"Mass = 91.9155");
      title->SetFillColor(0);
      title->SetTextFont(52);
      title->Draw();

      // h1.Print("all"); > histo.txt
      // cout << need<<endl;
      cout << "end";
      // TBrowser b;
      // h1.Draw();
   }



      
    