#include <iostream>
#include <fstream>
#include <string>
#include <sstream>
#include <cmath>
#include <ctime>
#include <vector>
#include <sstream>
#include <cstdlib>
#include <TFile.h>
#include <TDirectory.h>
#include <TCanvas.h>
#include <TTree.h>
#include <TH1F.h>
#include <TH2D.h>
#include <TF1.h>
#include <TMath.h>
#include <TCanvas.h>
#include <TLatex.h>
#include <TLegend.h>
#include <TStyle.h>
#include <TVector3.h>
#include <TString.h>
#include <TObjString.h>
#include <TObjArray.h>
#include <TKey.h>
#include <TIterator.h>
#include <TChain.h>
#include <TApplication.h>
#include <vector>
#include <TH3D.h>
#include <TRandom3.h>
//TApplication app("GUI",0,NULL);
 
using namespace std;

static double jd (time_t t) {
     struct tm *s = gmtime (&t);   // if t==0, year is 1970
     double year = s->tm_year + 1900;
     double month = s->tm_mon + 1;
     double day = s->tm_mday;
     double hour = s->tm_hour + s->tm_min / 60.;
     //cout<<year<<" "<<month<<" "<<day<<" "<<hour<<endl;
     return 367 * year - floor (7. / 4 * (floor ((month + 9) / 12) + year)) + floor (275 * month / 9) + day - 730531.5  + hour / 24;  //converts the time in days with repect to year 2000 (Julian epoch).
     }

static const double pi=3.1415926535897932384626433832795029L;
static double _2r (double grad) { return grad / 180 * pi; }
static double _2g (double rad) { return rad * 180 / pi; }

void get_sunposition(){
   
    int n_years = 10; 
    Long64_t t = (2024-1970)*365.25*24*60*60; // this should be given w.r.t 1970 to use jd()
    Long64_t t_end = (n_years+2024-1970)*365.25*24*60*60; // this should be given w.r.t 1970 to use jd()
   

    float altitude;
    float azimuth;
    int count=0;
    
    //TH2D* plot = new TH2D("plot","plot",260,-180,180,360,-180,180);
    TH2D* plot = new TH2D("plot","plot",200,-1,1,200,-1,1);
    TH1D* cos_theta = new TH1D("cos_theta","cos_theta",200,-1,1);
    TH1D* phi_h = new TH1D("phi_h","phi_h",200,-1,1);
    TH1D* theta_h = new TH1D("theta_h","theta_h",180,0,180);


    TString outfilename = "sunposition_since1Jan2024_for10years_JUNOlongandlat_forPDFcreation.root";
    //TString outfilename = "sunposition_since1Jan2024_for6years_eqandJUNOlong_forPDFcreation.root";

    //TString outfilename = "sunposition_since1Jan2019_for1years_eqandJUNOlong.root";    
    //TString outfilename = "sunposition_since1Jan2019_for1years_bxlatandlong.root";
    //TString outfilename = "sunposition_since1Jan2019_for1years_bxlatandJUNOlong.root";
    //TString outfilename = "sunposition_since1Jan2019_for1years_ArcticlatandJUNOlong.root";
    //TString outfilename = "sunposition_since1Jan2019_for1years_PrimeMeridianlongandJUNOlat.root";
    //TString outfilename = "sunposition_since1Jan2019_for1years_RandomlatandJUNOlong.root";
    //TString outfilename = "sunposition_since1Jan2019_for1years_NorthpolelatandJUNOlong.root";
    //TString outfilename = "sunposition_since1Jan2019_for1years.root";
    TFile *outfile=new TFile(outfilename,"RECREATE");
    TCanvas *c1 = new TCanvas("c1","c1");
      
    ofstream ofile("sunposition_since1Jan2024_for10years_PDFcreation.txt"); //, ofstream::app);
    //ofstream ofile("sunposition_since1Jan2024_for6years_PDFcreation.txt"); //, ofstream::app);
 
    while (t < t_end){
       
       //const double latitude = 42.421;    // bx
       //const double longitude = 13.51466; // bx
       
       //const double latitude = 42.421;    // bx
       //const double longitude = 112.51805555555555; // JUNO
       
       const double longitude = 112.51805555555555; // JUNO is at E 112° 31' 05"
       const double latitude = 22.118055555555557;    // JUNO is at N 22° 07' 05"
       
       //const double longitude = 112.51805555555555; // JUNO is at E 112° 31' 05"
       //const double latitude = 0;    // equator
       
       //const double longitude = 112.51805555555555; // JUNO is at E 112° 31' 05"
       //const double latitude = 90;    // North pole
       
       //const double longitude = 112.51805555555555; // JUNO is at E 112° 31' 05"
       //const double latitude = 23.5;    // Tropic of Cancer
       
       //const double longitude = 112.51805555555555; // JUNO is at E 112° 31' 05"
       //const double latitude = 12;    // Random
       
        //const double longitude = 112.51805555555555; // JUNO is at E 112° 31' 05"
        //const double latitude = 63.5;    // Arctic Circle
      
       //const double longitude = 0; // Prime Meridian
       //const double latitude = 22.118055555555557;    // JUNO is at N 22° 07' 05"
       
       double jj = jd (t);
       double T = jj / 36525; /* Julian Century */
        
       double M = fmod (357.5291 + 35999.0503 * T - 0.0001559 * T * T - 0.00000045 * T * T * T, 360);
       double Lo = fmod (280.46645 + 36000.76983 * T + 0.0003032 * T * T,  360);
       double DL = (1.9146 - 0.004817 * T - 0.000014 * T * T) * sin (_2r (M)) + (0.019993 - 0.000101 * T)* sin (_2r (2 * M)) + 0.00029 * sin (_2r (3 * M));
       double L = Lo + DL;
       
       double eps = 23.43999 - 0.013 * T;
       double delta=  _2g (asin (sin (_2r (L)) * sin (_2r (eps))));
       double RA = fmod (_2g (atan2 (cos (_2r (eps)) * sin (_2r (L)), cos (_2r (L)))) + 360, 360);
        
       double GMST = fmod (280.46061837 + 360.98564736629 * jj + 0.000387933 * T * T - T * T * T /38710000 + 360,  360);
       double LMST = GMST + longitude;
       double H = LMST - RA;
       azimuth = fmod (_2g (atan2 (-sin (_2r (H)), cos (_2r (latitude)) * tan (_2r (delta))- sin (_2r (latitude)) * cos (_2r (H)))) + 360,  360);
       altitude = _2g (asin (sin (_2r (latitude)) * sin (_2r (delta)) + cos (_2r (latitude)) * cos (_2r (delta))*cos(_2r (H))));
       
       //double phi = TMath::Pi()*(360-56.7-azimuth)/180.0;
       double phi = TMath::Pi()*(180-56.7-azimuth)/180.0;
       double theta = acos(cos(TMath::Pi()*(90.0 - altitude)/180.0));
       
       ofile<<azimuth<<" "<<altitude<<" "<<theta<<" "<<phi<<endl;
      
       t = t + 934.06;
       //t = t + 1000;
       
       double Xdir = -sin(theta)*cos(phi);
       double Ydir = -sin(theta)*sin(phi);
       double Zdir = -cos(theta);
      
       //if (cos(theta)>0.99 || cos(theta)<-0.99){
       //cout<<count<<" theta: "<<_2g(theta)<<" phi: "<<_2g(phi)<<endl;
       //cout<<count<<" altitude: "<<altitude<<" azimuth: "<<azimuth<<endl;
       //cout<<count<<" Xdir: "<<Xdir<<" Ydir: "<<Ydir<<" Zdir: "<<Zdir<<endl;
       //}
        
       //plot->Fill(_2g(TMath::ATan2(Ydir, Xdir)),_2g(theta));
       plot->Fill(TMath::ATan2(Ydir, Xdir)/TMath::Pi(),Zdir);
       cos_theta->Fill(Zdir);
       phi_h->Fill(TMath::ATan2(Ydir, Xdir)/TMath::Pi());
       theta_h->Fill(180*theta/TMath::Pi());
       count++;
     }
    
    
    
    c1->cd();
    plot->Draw("colz");
    c1->Update();
    outfile->cd();
    c1->Write();
    cos_theta->Write();
    phi_h->Write();
    theta_h->Write();
    outfile->Close();
    
}

int main(){
 
    //g++ -std=c++11 `root-config --cflags --glibs` get_sunposition.C -o get 
     
    get_sunposition();
 
    //cout<<jd(0)<<endl;
    //cout<<jd((2024-1970)*365.25*24*60*60)<<endl;
      
    return 0;
}