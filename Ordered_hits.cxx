#include <fstream>
#include <string>
#include <iostream>
#include <cstdlib>
#include <sstream>
#include <iomanip>

#include <TApplication.h>
#include <TVector3.h>
#include <TFile.h>
#include <TChain.h>
#include <TH1.h>
#include "TMath.h"
#include <TF1.h>
#include "TH1D.h"
#include "TH1.h"
#include "TH2F.h"
#include "TH3F.h"
#include <TROOT.h>
#include "TStyle.h"
#include <TCanvas.h>
#include <TRandom3.h>
#include <TSpectrum.h>
#include <TMultiGraph.h>
#include <TGraph.h>
#include <TGraphErrors.h>
#include <TLegend.h>
#include <TVirtualFitter.h>

#include <TRandom3.h>
#include <TFractionFitter.h>
#include <vector>
#include <TFitResult.h>
#include <TMatrixDSym.h>
#include "TMinuit.h"
#include "Functions.h"

#include <TString.h>

#ifndef ROOT_TH1D
#endif

#define PI 3.141592653589793
#include "TMath.h"
#include "Math/Util.h"

using namespace std;

double cos_alpha (double XPMT, double YPMT, double ZPMT, double XVertex, double YVertex, double ZVertex, double x2, double y2, double z2) {

    double PhX = XPMT/1000; //because in the PMT file the distances are stored in millimeters, in my tree in meters
    double PhY = YPMT/1000;
    double PhZ = ZPMT/1000;

    double x1 = PhX - XVertex;
    double y1 = PhY - YVertex;
    double z1 = PhZ - ZVertex;

    double modules = mod(x1,y1,z1)*mod(x2,y2,z2);
    double scalar = x1*x2 + y1*y2 + z1*z2 ;

    return  scalar/modules;
}

int main(int argc, char** argv) {
        
    if(argc!=5) {
            cout << "\n     USAGE:  Input_rootfile Output_rootfile Nth_hits Is_Bkg \n" << endl;
            return 1;
    }


    //Reading PMT positions

    std::vector<vector<double>> PMT_Position;

    ifstream ReadPMTPosition;
	ReadPMTPosition.open("PMTPos_CD_LPMT_onlyHama.csv");
	double blank;
	int Index;
	double x_PMT,y_PMT,z_PMT, r_PMT, theta_PMT, phi_PMT;
    int PMTNumber = 4996;

	for(int PMT=0;PMT<PMTNumber;PMT++){		
		ReadPMTPosition >> Index;
		ReadPMTPosition >> x_PMT;
		ReadPMTPosition >> y_PMT;
		ReadPMTPosition >> z_PMT;
		ReadPMTPosition >> blank >> blank;
		PMT_Position.push_back({x_PMT,y_PMT,z_PMT});			
	}	

    string Input_rootfile = argv[1];
    string Output_Rootfile = argv[2];
    int Nth_hits = atoi(argv[3]);
    int Is_Bkg = atoi(argv[4]);

    double XVertex,YVertex,ZVertex;
    double nu_theta, nu_phi;
    double nu_x,nu_y,nu_z;

    std::vector<double>* Start_Time = new std::vector<double>;
    std::vector<bool>* Type = new std::vector<bool>;
    std::vector<bool>* Hit = new std::vector<bool>;
    std::vector<double>* PhX = new std::vector<double>;
    std::vector<double>* PhY = new std::vector<double>;
    std::vector<double>* PhZ = new std::vector<double>;
    std::vector<int>* ClosestPMT = new std::vector<int>;

    cout << "Reading " << Input_rootfile  << endl;
    cout << "Nth_hits = " << Nth_hits << endl << endl;

    TFile *f = new TFile (Input_rootfile.c_str());
    TTree* tree= (TTree*)f->Get("t");

    tree -> SetBranchAddress("Start_Time",&Start_Time);
    tree -> SetBranchAddress("Type",&Type);
    tree -> SetBranchAddress("Hit",&Hit);
    tree -> SetBranchAddress("Ph_x_AtPMT",&PhX);
    tree -> SetBranchAddress("Ph_y_AtPMT",&PhY);
    tree -> SetBranchAddress("Ph_z_AtPMT",&PhZ);
    tree -> SetBranchAddress("Int_Vertex_x",&XVertex);
    tree -> SetBranchAddress("Int_Vertex_y",&YVertex);
    tree -> SetBranchAddress("Int_Vertex_z",&ZVertex);
    tree -> SetBranchAddress("Closest_PMT",&ClosestPMT);
    tree -> SetBranchAddress("Neutrino_theta",&nu_theta);
    tree -> SetBranchAddress("Neutrino_phi",&nu_phi);

    int TotalEvents = tree -> GetEntries();

    TH1F *ChScRatio = new TH1F("ChScRatio","ChScRatio",Nth_hits,0,Nth_hits);
    TH1F *Cherenkov_cos_alpha = new TH1F("Cherenkov_cos_alpha","Cherenkov_cos_alpha",60,-1,1);
    TH1F *Scint_cos_alpha = new TH1F("Scint_cos_alpha","Scint_cos_alpha",60,-1,1);
    TH1F *Scint_cos_alpha_all = new TH1F("Scint_cos_alpha_all","Scint_cos_alpha_all",60,-1,1);
    TH1F *Cherenkov_cos_alpha_firstn = new TH1F(TString::Format("Cherenkov_cos_alpha_first%i",Nth_hits),TString::Format("Cherenkov_cos_alpha_first%i",Nth_hits),60,-1,1);
    TH1F *Scint_cos_alpha_firstn = new TH1F(TString::Format("Scint_cos_alpha_first%i",Nth_hits),TString::Format("Scint_cos_alpha_first%i",Nth_hits),60,-1,1);
    TH1F *Cherenkov_cos_alpha_8th_hit = new TH1F("Cherenkov_cos_alpha_8th_hit","Cherenkov_cos_alpha_8th_hit",60,-1,1);

    TH1F *Nth_hit_cos_alpha[Nth_hits];

    for (int i=0;i<Nth_hits;i++) {
        if (Is_Bkg == 0) {
        	Nth_hit_cos_alpha[i] = new TH1F(TString::Format("Nsolar_%i",i+1),TString::Format("Nsolar_%i",i+1),60,-1,1);
	    }
	    else {
		    Nth_hit_cos_alpha[i] = new TH1F(TString::Format("Nbkg_%i",i+1),TString::Format("Nbkg_%i",i+1),60,-1,1);
	    }
    }

    int* FirstTen = new int[Nth_hits]; //contains the number of times a Cherenkov photon arrived in the n-th position
	double* FirstTenValues = new double[Nth_hits]; //contains the time at which the first n-th photon arrived (could be Cherenkov or scintillation)
	int* FirstTenPlaces = new int[Nth_hits]; //contains the place in the vector of the n-th photon to arrive (Cher or scint)

	for (int i=0; i<Nth_hits; i++) { //cleans the counter
		FirstTen[i] = 0;
	}


    for (int i=0;i<TotalEvents;i++) {

        tree -> GetEntry(i);

        int vec_lenght = Start_Time -> size();

        vector <double> val_cos_alpha;
        
        for (int k=0; k<Nth_hits; k++) {

		    FirstTenValues[k] = 10000;
		    FirstTenPlaces[k] = -1; //useless, just to be sure that it is updating correctly into the cycle
        }
       
        for (int j=0;j<vec_lenght;j++) {
           
            SphericalToCartesian(nu_x,nu_y,nu_z,1.,nu_theta,nu_phi);

            if (Type -> at(j) == 0 && Hit -> at(j) == 1) {
                Scint_cos_alpha -> Fill(cos_alpha(PMT_Position[ClosestPMT -> at(j)][0],PMT_Position[ClosestPMT -> at(j)][1],PMT_Position[ClosestPMT -> at(j)][2],XVertex,YVertex,ZVertex,nu_x,nu_y,nu_z));
            } else if (Type -> at(j)==1 && Hit -> at(j) == 1) {
                Cherenkov_cos_alpha -> Fill(cos_alpha(PMT_Position[ClosestPMT -> at(j)][0],PMT_Position[ClosestPMT -> at(j)][1],PMT_Position[ClosestPMT -> at(j)][2],XVertex,YVertex,ZVertex,nu_x,nu_y,nu_z));
            } else if (Type -> at(j) == 0 ) {
                Scint_cos_alpha_all -> Fill(cos_alpha(PhX -> at(j)*1000.,PhY -> at(j)*1000. ,PhZ -> at(j)*1000.,XVertex,YVertex,ZVertex,nu_x,nu_y,nu_z));
            }

            if (ClosestPMT -> at(j) != -1 ) val_cos_alpha.push_back(cos_alpha(PMT_Position[ClosestPMT -> at(j)][0],PMT_Position[ClosestPMT -> at(j)][1],PMT_Position[ClosestPMT -> at(j)][2],XVertex,YVertex,ZVertex,nu_x,nu_y,nu_z));
            else val_cos_alpha.push_back(0);

            for (int k=0; k<Nth_hits; k++) {	//looks if a photon is in the first N to arrive and updates the arrays accordingly

					if (Start_Time -> at(j) < FirstTenValues[k] && Hit -> at(j) ==1) {
					
						if (k != Nth_hits -1 ) {
							
							int* provvint = new int[Nth_hits-1];
							double* provvdouble = new double[Nth_hits-1];
							int cost =0;

							for (int h=k;h<Nth_hits-1;h++) {			
								provvint[cost] = FirstTenPlaces[h];
								provvdouble[cost] = FirstTenValues[h];
								cost ++;
								
							}
							cost=0;
							for (int h=k;h<Nth_hits-1;h++) {
								FirstTenPlaces[h+1] = provvint[cost];
								 FirstTenValues[h+1] = provvdouble[cost];
								cost ++;

							}
						}	
						

						FirstTenPlaces[k] = j;
						FirstTenValues[k] = Start_Time -> at(j) ;
						break;
					}
				
				}

        } 
        
        for (int k=0;k<Nth_hits;k++) { //updates the counter

			if ( Type -> at(FirstTenPlaces[k]) == 1) {
				FirstTen[k]++;
                Cherenkov_cos_alpha_firstn -> Fill(val_cos_alpha[FirstTenPlaces[k]]);
			} else {
                Scint_cos_alpha_firstn -> Fill(val_cos_alpha[FirstTenPlaces[k]]);
            }
            
            for (int h=0;h<Nth_hits;h++) {
                if (k == h) {
                    Nth_hit_cos_alpha[h] -> Fill(val_cos_alpha[FirstTenPlaces[k]]);

                    if (k==7 && Type -> at(FirstTenPlaces[k]) == 1) {
                        Cherenkov_cos_alpha_8th_hit ->  Fill(val_cos_alpha[FirstTenPlaces[k]]);
                    }

                } 
            }

		}

        val_cos_alpha.clear();
        
    }

    for (int k=0;k<Nth_hits;k++) {

		ChScRatio -> Fill(k,double(FirstTen[k])/(TotalEvents));

    }

    cout << "Creating " << Output_Rootfile << endl; 

    TFile* out = new TFile(Output_Rootfile.c_str(),"recreate");

    ChScRatio -> Write();
    Cherenkov_cos_alpha -> Write();
    Scint_cos_alpha -> Write();
    Scint_cos_alpha_all -> Write();
    ChScRatio -> Write();
    Scint_cos_alpha_firstn -> Write();
    Cherenkov_cos_alpha_firstn -> Write();
    Cherenkov_cos_alpha_8th_hit -> Write();

    
    for (int i=0; i<Nth_hits;i++) {
        Nth_hit_cos_alpha[i] -> Write();
    }

    

    out -> Close();
    
    return 0;
    
} 