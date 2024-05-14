#include <fstream>
#include <string>
#include <iostream>
#include <cstdlib>
#include <sstream>
#include <iomanip>

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
#include <TApplication.h>

#include <TRandom3.h>
#include <TFractionFitter.h>
#include <vector>
#include <TFitResult.h>
#include <TMatrixDSym.h>
#include "TMinuit.h"

#include <TString.h>

#include "Functions.h"

#ifndef ROOT_TH1D
#endif

#define PI 3.141592653589793
#include "TMath.h"
#include "Math/Util.h"

int PMTNumber;
double PMTRadius = 0.25;
double JUNORadius = 19.434;
double FV;
int LY;
double ChScRatio;
int NEvents;
int TotalPhotons = 0;
TH1D** Time_PDFs = new TH1D*[3]; //container for the time PDFs, pos 0 for Cherenkov, pos 1 for scintillation, pos 3 for DN
TH1D** B8_PDFs = new TH1D*[1]; 

//Useful values

double n = 1.5 ; //refraction index
double c = 299792458 ; // m/s
double m_e = 0.51099895; //MeV    electron mass
double Be7_energy = 0.862; //MeV    energy of a 7Be neutrino
double pep_energy = 1.44; //MeV  energy of a pep neutrino
double G_F = 1.1663787*pow(10,-11); //MeV^-2  Fermi constant
double sin2_thetaW = 0.23121; //Weinberg angle

double min_eEnergy; // minimum deposited energy from a neutrino scattering -> set globally by the energy window 

double RefractionIndex = 1.5;

int PEatMeV; //photoelectron emitted @ 1 MeV
int DN_per_event;

double El_Direction_x_t,El_Direction_y_t,El_Direction_z_t,phi_t,theta_t,Start_Time_t,Electron_Energy_t,Neutrino_Energy_t,type_t;
double Ph_x_AtPMT_t,Ph_y_AtPMT_t,Ph_z_AtPMT_t, Ph_r_AtPMT_t, Ph_theta_AtPMT_t, Ph_phi_AtPMT_t, TravelledDistance_t;
double Int_Vertex_x_t,Int_Vertex_y_t,Int_Vertex_z_t;
int Type_t, Closest_PMT_t;
bool IsFirst_t;
int NEvent_t;
double Cher_angle_t,Elec_angle_t;
double Neutrino_theta_t;
double Neutrino_phi_t;
double Solar_theta_t, Solar_phi_t;

std::vector <double> *Ph_x_AtPMT_v;
std::vector <double> *Ph_y_AtPMT_v;
std::vector <double> *Ph_z_AtPMT_v;
std::vector <double> *TravelledDistance_v;
std::vector <double> *Ph_theta_AtPMT_v;
std::vector <double> *Ph_phi_AtPMT_v;
std::vector <int> *Closest_PMT_v;
std::vector <double> *Start_Time_v;
std::vector <int> *Type_v;
std::vector <bool> *Hit_v;

std::vector<std::vector<double>> PMT_Position_Spherical;	
std::vector<std::vector<std::tuple<int,double,double,double>>> PMT_Position_Fast;	

bool Hit_t;
bool fastmode; //saves only a limited number of hits
double TimeCut; //cut photons emitted at times lower than timecut
bool FixedSun, IsBackgrounds;
int EventLimit; //the number of scintillation hits produced
std::string typenu;
double WindowEndpoint;

using namespace std;

double CalculateNuEnergy () {

	double neutrino_energy, max_eEnergy;

	if (typenu == "Be7") {

		neutrino_energy = Be7_energy;

	} else if (typenu == "pep") {

		neutrino_energy = pep_energy;

	} else if (typenu == "B8" ) {

		while (true) {

			neutrino_energy = B8_PDFs[0] -> GetRandom();
			max_eEnergy = 2*pow(neutrino_energy,2)/(m_e+2*neutrino_energy); 

			if (max_eEnergy > min_eEnergy) { //checks if the neutrino has enough energy to emit an electron inside the ROI
				break;
			} 

		}
		
	} else {
		cout << "ERROR : invalid nu type in Configuration File" << endl;
		exit(1);
	}

	return neutrino_energy;

}

//calculate 1st order cross section for a neutrino-electron elastic scattering
double cross_section (double T, double nu_energy) {
	double gl = 0.5 - sin2_thetaW -1 ;
	double gr = - sin2_thetaW;

	return ((2*(pow(G_F,2))*m_e)/M_PI)*(pow(gl,2) + pow(gr,2)*pow(1-T/nu_energy,2) - gl*gr*m_e*T/pow(nu_energy,2));
}

double cross_section_nuX (double T, double nu_energy) {
	double g_L = 0.5 - sin2_thetaW;
	double g_R = -sin2_thetaW;

	return 2*m_e*G_F*G_F/M_PI * (pow(g_L,2) + pow(g_R,2)*pow(1-T/nu_energy,2) - g_L*g_R*m_e/nu_energy*T/nu_energy);
}

double survival_probability (double nu_energy) {

	double E_nu = nu_energy*pow(10,6);

	double V_x, Delta_V_x2__V_x2;
	//datas from pep, see Xuefeng technote
	if (typenu == "pep") {
		V_x = 5.13e-12;
		Delta_V_x2__V_x2 = 0.076;
	} else if (typenu == "Be7") {
		V_x = 6.16e-12;
		Delta_V_x2__V_x2 = 0.029;
	} else if (typenu == "B8") {
		V_x = 6.81e-12;
		Delta_V_x2__V_x2 = 0.010;
	} else {
		throw "ERROR: invalid typenu";
	}

	double theta12 = asin(sqrt(0.307));
	double theta13 = asin(sqrt(0.0220));
	double Delta_m2_12 = 7.53e-5; 

	double epsilon12 = 2*V_x*E_nu/Delta_m2_12;
	double alpha = cos(2*theta12) - ( pow(cos(theta13),2) * epsilon12 );
	double cos_2theta_m12 = alpha / sqrt(alpha*alpha + pow(sin(2*theta12),2) );
	double d_x = 3/2 * (pow(epsilon12,2)*pow(sin(2*theta12),2))/(pow(pow(cos(2*theta12) - epsilon12 ,2) + pow(sin(2*theta12),2),2)) * Delta_V_x2__V_x2;
	double Pad_2 = 0.5 + 0.5*(1-d_x)*cos(2*theta12)*cos_2theta_m12;

	return pow(cos(theta13),4)*Pad_2 + pow(sin(theta13),4);

}
double total_cross_section (double T, double nu_energy) {

	return survival_probability(nu_energy)*cross_section(T,nu_energy) + (1-survival_probability(nu_energy))*cross_section_nuX(T,nu_energy);

}

//sample cross_section with an accept-reject method
double CalculateEventEnergy (double nu_energy) {
	double EventEnergy, test;
	bool flag = false;

	double max_eEnergy = 2*pow(nu_energy,2)/(m_e+2*nu_energy); 

	if (max_eEnergy > nu_energy ) {
		cerr << "SOMETHING WENT WRONG: max electron energy > neutrino energy " << endl;
		exit(1);
	}

	if (IsBackgrounds == true) {
		return gRandom -> Uniform(min_eEnergy,WindowEndpoint);
	}

	while (flag == false) {
		flag = false;
		EventEnergy = gRandom -> Uniform(min_eEnergy,max_eEnergy);
		test = gRandom -> Uniform(0.,total_cross_section(0.,nu_energy)); //the maximum cross section is at T=0
		if (test < total_cross_section(EventEnergy,nu_energy)) {
			flag = true;
		}
	}

	return EventEnergy;
}

int ClosestPMTIndex (double x_Event,double y_Event,double z_Event) {

	double Distance_Temp;
	double r_Event, theta_Event, phi_Event;
	double x_PMT,y_PMT,z_PMT, r_PMT, theta_PMT, phi_PMT;

	CartesianToSpherical(r_Event, theta_Event, phi_Event,x_Event, y_Event, z_Event);

	//std::cout << "Event_in = " << theta_Event << "  " << phi_Event << endl; 
		
	for(int PMT_ring=0;PMT_ring<PMT_Position_Fast.size();PMT_ring++){

		//r_PMT = get<0>(PMT_Position_Fast[PMT_ring][0]);
		//theta_PMT = get<1>(PMT_Position_Fast[PMT_ring][0]);
		phi_PMT = get<3>(PMT_Position_Fast[PMT_ring][0]);
		//Distance_Temp = DistanceOnASphere(r_PMT, theta_PMT, phi_PMT,theta_Event, phi_Event);

		if ( (abs (phi_PMT - phi_Event) + 1e-5)*JUNORadius <= PMTRadius ) {

			//cout << "IN ring #" << PMT_ring << endl << endl;

			for (int PMT_NuminRing=0; PMT_NuminRing < PMT_Position_Fast[PMT_ring].size(); PMT_NuminRing++) {

				r_PMT = get<1>(PMT_Position_Fast[PMT_ring][PMT_NuminRing]);
				theta_PMT = get<2>(PMT_Position_Fast[PMT_ring][PMT_NuminRing]);
				phi_PMT = get<3>(PMT_Position_Fast[PMT_ring][PMT_NuminRing]);
				Distance_Temp = DistanceOnASphere(theta_PMT, phi_PMT,theta_Event, phi_Event);

				//cout << "       From #" << PMT_NuminRing << " =  " << Distance_Temp*r_PMT << endl;

					if (Distance_Temp*r_PMT/1000 <= PMTRadius ) {
						return get<0>(PMT_Position_Fast[PMT_ring][PMT_NuminRing]);
					}
			}
		}
	}

	return -1;

}

//outputs the position (x,y,z) on the sphere for a photon generated in (j,k,l) with direction (a,b,c)

void Position_on_PMTsphere (double & x,double& y , double& z, double j, double k, double l, double a, double b , double c ) {

	double sroot,num,denom,t;

	denom = 2*(a*a + b*b + c*c);
	sroot = pow(2*j*a + 2*k*b + 2*l*c,2) - 4*(a*a + b*b + c*c)*(j*j + k*k + l*l - JUNORadius*JUNORadius);
	num= -(2*j*a + 2*k*b + 2*l*c) + sqrt(sroot);

	t = num/denom;
	x = j+a*t;
	y = k+b*t;
	z = l+c*t;

	return;
}

//generates a randomly distributed photon with an angle to the beginning direction
std::tuple<double,double> Generate_Cone (double theta_0, double phi_0, double angle) {

	std::tuple<double,double> out;
    double phi_out, theta_out;
    double theta, phi, cos_th1, cos_ph1, sin_th1, sin_ph_sin_th, sin_ph_cos_th;

	//generate random values on a cone centered in 0,0
    theta = gRandom->TRandom::Uniform(2*PI); 
    phi = angle; 

	//some exceptions
    if (phi_0 == 0) {

        phi_out = phi;
        theta_out = theta;
        
    } else if (phi_0 == M_PI) {

        phi_out = M_PI - phi;
        theta_out = theta;

	// traslate the cone on the original system
    } else { 

        cos_ph1 = -sin(phi_0)*cos(theta)*sin(phi) + cos(phi_0)*cos(phi);
        sin_ph_cos_th = cos(phi_0)*cos(theta_0)*cos(theta)*sin(phi) - sin(theta_0)*sin(theta)*sin(phi) + cos(theta_0)*sin(phi_0)*cos(phi);
        sin_ph_sin_th = sin(theta_0)*cos(phi_0)*cos(theta)*sin(phi) + cos(theta_0)*sin(theta)*sin(phi) + sin(theta_0)*sin(phi_0)*cos(phi);        

        phi_out = acos(cos_ph1);
        cos_th1 = sin_ph_cos_th/sin(phi_out);
        sin_th1 = sin_ph_sin_th/sin(phi_out);

        if (sin_th1 >= 0) {
            theta_out = acos(cos_th1) ;
        } else if (sin_th1 < 0) {
            theta_out = Pbc_theta(-acos(cos_th1));
        }

    }

	out = make_tuple(theta_out,phi_out);

	return out;
}

int CheckHit (int SeenPhotons) {

	if (Closest_PMT_t != -1) {
		Hit_t = 1;
		SeenPhotons++;

	} else {Hit_t = 0;}


	return SeenPhotons;
}

std::vector<tuple<double,int>> GenerateTimes (int NPhotonsScint, int NPhotonsCher, int NPhotonsDN) {

	vector<tuple<double,int>> Times;
	double Test_Time;
	for (int i = 0; i < NPhotonsScint; i++) {
		do {
			Test_Time = Time_PDFs[0] -> GetRandom();
		} while (Test_Time < TimeCut);
		Times.push_back(make_tuple(Test_Time,0));
		
	}
	for (int i = 0; i < NPhotonsCher; i++) {
		do {
			Test_Time = Time_PDFs[1] -> GetRandom();
		} while (Test_Time < TimeCut);
		Times.push_back(make_tuple(Test_Time,1));
		
	}
	for (int i = 0; i < NPhotonsDN; i++) {
		do {
			Test_Time = Time_PDFs[2] -> GetRandom();
		} while (Test_Time < TimeCut);
		Times.push_back(make_tuple(Test_Time,2));
		
	}
	std::sort(Times.begin(),Times.end());

	/*for (int itr = 0; itr < 40; itr++) {
		cout << get<1>(Times[itr]) << "  " << get<0>(Times[itr]) << endl;
	}*/


	return Times;
}

void GenerateScintHit () {

	double Old_Closest_PMT, New_Closest_PMT;

	double ph_Direction_x,ph_Direction_y,ph_Direction_z,trash;

	double theta_vers = gRandom->TRandom::Uniform(2*PI);
	double phi_vers = TMath::ACos(-1.+2.*gRandom->TRandom::Uniform(0,1));
	SphericalToCartesian(ph_Direction_x,ph_Direction_y,ph_Direction_z,1,theta_vers,phi_vers);

	Position_on_PMTsphere(Ph_x_AtPMT_t,Ph_y_AtPMT_t,Ph_z_AtPMT_t,Int_Vertex_x_t,Int_Vertex_y_t,Int_Vertex_z_t,ph_Direction_x,ph_Direction_y,ph_Direction_z);

	CartesianToSpherical(trash,theta_t,phi_t,Ph_x_AtPMT_t,Ph_y_AtPMT_t,Ph_z_AtPMT_t);

	Ph_r_AtPMT_t = JUNORadius; 
	Ph_theta_AtPMT_t = Pbc_theta(theta_t); 
	Ph_phi_AtPMT_t = Pbc_phi(phi_t); 

	Closest_PMT_t = ClosestPMTIndex(Ph_x_AtPMT_t,Ph_y_AtPMT_t,Ph_z_AtPMT_t);

	TravelledDistance_t = Distance(Int_Vertex_x_t,Int_Vertex_y_t,Int_Vertex_z_t,Ph_x_AtPMT_t,Ph_y_AtPMT_t,Ph_z_AtPMT_t);

}

void GenerateCherenkovHit (double El_theta, double El_phi, double theta_Cher) {

	double ph_Direction_x,ph_Direction_y,ph_Direction_z,trash;

	std::tuple <double,double> Cher_phot_dir = Generate_Cone(El_theta,El_phi,theta_Cher);

	double theta_vers = get<0>(Cher_phot_dir);
	double phi_vers = get<1>(Cher_phot_dir);

	SphericalToCartesian(ph_Direction_x,ph_Direction_y,ph_Direction_z,1,theta_vers,phi_vers);

	Position_on_PMTsphere(Ph_x_AtPMT_t,Ph_y_AtPMT_t,Ph_z_AtPMT_t,Int_Vertex_x_t,Int_Vertex_y_t,Int_Vertex_z_t,ph_Direction_x,ph_Direction_y,ph_Direction_z);

	CartesianToSpherical(trash,theta_t,phi_t,Ph_x_AtPMT_t,Ph_y_AtPMT_t,Ph_z_AtPMT_t);

	Ph_r_AtPMT_t = JUNORadius; 
	Ph_theta_AtPMT_t = Pbc_theta(theta_t); 
	Ph_phi_AtPMT_t = Pbc_phi(phi_t); 

	Closest_PMT_t = ClosestPMTIndex(Ph_x_AtPMT_t,Ph_y_AtPMT_t,Ph_z_AtPMT_t);

	TravelledDistance_t = Distance(Int_Vertex_x_t,Int_Vertex_y_t,Int_Vertex_z_t,Ph_x_AtPMT_t,Ph_y_AtPMT_t,Ph_z_AtPMT_t);

}

void GenerateDNHit() {

	Closest_PMT_t = gRandom -> TRandom::Uniform(PMTNumber);
	Ph_theta_AtPMT_t = PMT_Position_Spherical[Closest_PMT_t][1];
	Ph_phi_AtPMT_t = PMT_Position_Spherical[Closest_PMT_t][2];

	SphericalToCartesian(Ph_x_AtPMT_t,Ph_y_AtPMT_t,Ph_z_AtPMT_t,JUNORadius,Ph_theta_AtPMT_t,Ph_phi_AtPMT_t);
	TravelledDistance_t = Distance(Int_Vertex_x_t,Int_Vertex_y_t,Int_Vertex_z_t,Ph_x_AtPMT_t,Ph_y_AtPMT_t,Ph_z_AtPMT_t);

}

//Generator for all the photons
int GenerateEvents (TTree* t, bool RandomPos, int NEvent, ifstream& ReadSolarPosition) {

	double x_Int,y_Int,z_Int,r_Int,theta_Int,phi_Int;
	double theta_vers,phi_vers,trash;
	int SeenPhotons = 0;
	double ph_Direction_x,ph_Direction_y,ph_Direction_z;
	//Generate a random interaction vertex

	if (RandomPos == true ) {
		r_Int = gRandom -> TRandom::Uniform(FV);
		theta_Int = gRandom -> TRandom::Uniform(2*PI);
		phi_Int = TMath::ACos(-1.+2.*gRandom->TRandom::Uniform(0,1));
		SphericalToCartesian(x_Int,y_Int,z_Int,r_Int,theta_Int,phi_Int);
	} else {
		x_Int = 0;
		y_Int = 0;
		z_Int = 0;
	}

	//Randomize the energy of the event
	double nu_energy = CalculateNuEnergy();
	double Event_Energy = CalculateEventEnergy(nu_energy);
	
	double theta_e = acos((1+m_e/nu_energy)*pow(Event_Energy/(Event_Energy+2*m_e),0.5)); //angle between the solar-nu and the electron scattered (assuming 7Be-nu)
	double beta_el = pow(1-(pow(m_e/(Event_Energy+m_e),2)),0.5) ; //beta of the electron generated
	double theta_Cher = acos(1/(beta_el*n)); //Cherenkov angle

	int ScintPhotons = gRandom -> Poisson(PEatMeV*Event_Energy);
	int CherenkovPhotons = gRandom -> Poisson(ChScRatio*ScintPhotons);
	int DNPhotons = gRandom -> Poisson(DN_per_event);
	int EffectivePhotonsGenerated;

	int EventPhotons = ScintPhotons + CherenkovPhotons + DNPhotons;

	if (fastmode) EffectivePhotonsGenerated = EventLimit;
	else EffectivePhotonsGenerated = EventPhotons;

	if (NEvent < 10) {
		cout << "Event #" << NEvent << " :  " << ScintPhotons << " scintillation and " << CherenkovPhotons << " Cherenkov photons emitted" << endl;
		cout << "        Dark Noise events = " << DNPhotons << endl;
		cout << "        Nu energy = " << nu_energy << "    Electron energy = " << Event_Energy << endl;  
		cout << "        Theta_e = " << theta_e << "    Theta_Cher = " << theta_Cher << endl << endl; 
	} else if (NEvent == 10) {
		cout<<endl;
	}

	bool IsFirstFlag = true;

	double Provv_StartTime;

	//Reading the solar positions
	double solar_nu_theta, solar_nu_phi;
	double solar_theta, solar_phi; //Direction opposite to the sun's position

	if (FixedSun) {
		
		solar_nu_theta = 0.;
		solar_nu_phi = 0.;
		solar_theta = 0.;
		solar_phi = 0.;

	} else {
		double solar_az,solar_alt;
		double theta_a,phi_a;
		
		if (ReadSolarPosition.peek() != EOF) ReadSolarPosition >> solar_az >> solar_alt >> theta_a >> phi_a ;
		else {
			cerr << "ERROR: Solar positions terminated, you generated more than 10yr of events" << endl;
			exit(1);
		}

		solar_phi = Pbc_phi(M_PI_2 - DegToRad( solar_alt ) );
		solar_theta = Pbc_theta (M_PI_2 - 0.989601686 - DegToRad (solar_az) ); //0.9... is the angle between the N pole and JUNO coordinates

		if (IsBackgrounds) {
			solar_nu_theta = gRandom -> TRandom::Uniform(2*PI);
			solar_nu_phi = TMath::ACos(-1.+2.*gRandom->TRandom::Uniform(0,1));
		} else {

			solar_nu_phi =  Pbc_phi(M_PI - solar_phi);
			solar_nu_theta = Pbc_theta( M_PI + solar_theta);

		}
		
	}

	//Generate Electron direction
	std::tuple<double,double> El_dir,Cher_phot_dir;

	El_dir = Generate_Cone(solar_nu_theta,solar_nu_phi,theta_e);

	double El_theta = get<0>(El_dir);
	double El_phi = get<1>(El_dir);

	double El_x , El_y, El_z;

	SphericalToCartesian(El_x,El_y,El_z,1.,El_theta,El_phi);

	//Set all event-level variables 
	Cher_angle_t = theta_Cher;
	Elec_angle_t = theta_e;
	Electron_Energy_t = Event_Energy;
	Neutrino_Energy_t = nu_energy;

	Int_Vertex_x_t = x_Int;
	Int_Vertex_y_t = y_Int;
	Int_Vertex_z_t = z_Int;

	El_Direction_x_t = El_x;
	El_Direction_y_t = El_y;
	El_Direction_z_t = El_z;

	Neutrino_theta_t = solar_nu_theta;
	Neutrino_phi_t = solar_nu_phi;
	Solar_phi_t = solar_phi;
	Solar_theta_t = solar_theta;

	//initialize vectors 

	Ph_x_AtPMT_v = new vector<double> ();
	Ph_y_AtPMT_v = new vector<double> ();
	Ph_z_AtPMT_v = new vector<double> ();
	TravelledDistance_v = new vector<double> ();
	Ph_theta_AtPMT_v = new vector<double> ();
	Ph_phi_AtPMT_v = new vector<double> ();
	Closest_PMT_v = new vector<int> ();
	Start_Time_v = new vector<double> ();
	Type_v = new vector<int> ();
	Hit_v = new vector<bool> ();

	//Generate SCINTILLATION Photons

	std::vector<tuple<double,int>> Times = GenerateTimes(ScintPhotons,CherenkovPhotons,DNPhotons); //generates and sorts the photons times

	for(int iPh=0;; iPh++){ //continues until it breaks because it generated enough photons

		//Generate unit vector over a sphere
		Start_Time_t = get<0>(Times[iPh]);
		Type_t = get<1>(Times[iPh]);
	
		if (Type_t == 0) {
			GenerateScintHit();
		} else if (Type_t == 1) {
			GenerateCherenkovHit(El_theta,El_phi,theta_Cher);
		} else if (Type_t == 2) {
			GenerateDNHit();
		} else {
			cerr << "ERROR: Unexpected photon type" << endl;
			exit(1); 
		}

		//cout << "Event generated: Type = "<< Type_t << "     Emission Time  = " << Start_Time_t << endl;  
		
		//Type_t = 0; 

		//Check if the photon hits a PMT or not, SeenPhotons counts the number of photon seen in the event
		SeenPhotons = CheckHit(SeenPhotons); //it also sets the value of Hit_t

		TotalPhotons++; 

		//fill vectors 

		Ph_x_AtPMT_v -> push_back(Ph_x_AtPMT_t);
		Ph_y_AtPMT_v -> push_back(Ph_y_AtPMT_t);
		Ph_z_AtPMT_v -> push_back(Ph_z_AtPMT_t);
		TravelledDistance_v -> push_back(TravelledDistance_t);
		Ph_theta_AtPMT_v -> push_back(Ph_theta_AtPMT_t);
		Ph_phi_AtPMT_v -> push_back(Ph_phi_AtPMT_t);
		Closest_PMT_v -> push_back(Closest_PMT_t);
		Start_Time_v -> push_back(Start_Time_t);
		Type_v -> push_back(Type_t);
		Hit_v -> push_back(Hit_t);

		if (fastmode && SeenPhotons == EffectivePhotonsGenerated) {
			break;
		} else if (iPh == EventPhotons - 1) {
			if (fastmode) {
				std::cerr << "ERROR in event #" << NEvent << ", not enough hits generated (" << SeenPhotons << "/" << EffectivePhotonsGenerated << ")" << std::endl;
			}
			break;
		}
						
	}

	t -> Fill();

	Times.clear();
	Ph_x_AtPMT_v -> clear();
	Ph_y_AtPMT_v -> clear();
	Ph_z_AtPMT_v -> clear();
	TravelledDistance_v -> clear();
	Ph_theta_AtPMT_v -> clear();
	Ph_phi_AtPMT_v -> clear();
	Closest_PMT_v -> clear();
	Start_Time_v -> clear();
	Type_v -> clear();
	Hit_v -> clear();
	
	return SeenPhotons;
	
}


double Directionality_ToyMC(string Configuration_Text, string Output_Rootfile) {

	ifstream file(Configuration_Text);
	//ofstream temp_out("Theta_e.txt");
	vector<string> col1;
	vector<string> col2;
	string line;
	string cher_times, scint_times, Boron_PDFs_file, Boron_PDFs_name, DN_times;
	string origin_rootfile, PMTPositions, SolarPositions;
	bool RandomIntVertex;

	// ### Parsing
	cout << "######### Configuration #########" << endl;
	cout << "Cfg file: " << Configuration_Text.c_str() << endl;
	while (getline(file, line)) {
		string s1, s2;
		istringstream iss(line);
		if (!(iss >> s1 >> s2)) { break; } // error
			col1.push_back(s1);
			col2.push_back(s2);
			cout << setw(11) << s1 <<  "\t\t" << s2 << endl;
	}
	cout << "################################" << endl;

	istringstream iss(col2[0]);
	iss >> PEatMeV;
	istringstream iss1(col2[1]);
	iss1 >> ChScRatio;
	istringstream iss2(col2[2]);	
	iss2 >> NEvents;
	istringstream iss3(col2[3]);	
	iss3 >> FV;	
	istringstream iss4(col2[4]);	
	iss4 >> origin_rootfile;
	istringstream iss5(col2[5]);	
	iss5 >> cher_times;
	istringstream iss6(col2[6]);	
	iss6 >> scint_times;
	istringstream iss7(col2[7]);	
	iss7 >> typenu;
	istringstream iss8(col2[8]);	
	iss8 >> fastmode;
	istringstream iss9(col2[9]);	
	iss9 >> min_eEnergy;
	istringstream iss10(col2[10]);	
	iss10 >> RandomIntVertex;
	istringstream iss11(col2[11]);	
	iss11 >> TimeCut;
	istringstream iss12(col2[12]);	
	iss12 >> PMTPositions;
	istringstream iss13(col2[13]);	
	iss13 >> PMTNumber;
	istringstream iss14(col2[14]);
	iss14 >> SolarPositions;
	istringstream iss15(col2[15]);
	iss15 >> FixedSun;
	istringstream iss16(col2[16]);
	iss16 >> IsBackgrounds;
	istringstream iss17(col2[17]);
	iss17 >> EventLimit;
	istringstream iss18(col2[18]);
	iss18 >> WindowEndpoint;
	istringstream iss19(col2[19]);
	iss19 >> DN_times;
	istringstream iss20(col2[20]);
	iss20 >> DN_per_event;
	istringstream iss21(col2[21]);
	iss21 >> Boron_PDFs_file;
	istringstream iss22(col2[22]);
	iss22 >> Boron_PDFs_name;

	// ### End parsing	

	TFile *rootfile = new TFile(origin_rootfile.c_str());
	
	Time_PDFs[0] = (TH1D*)rootfile->Get(scint_times.c_str());
	Time_PDFs[1] = (TH1D*)rootfile->Get(cher_times.c_str());
	Time_PDFs[2] = (TH1D*)rootfile->Get(DN_times.c_str());

	if (typenu == "B8") {

		TFile *B8_rootfile = new TFile(Boron_PDFs_file.c_str());
		B8_PDFs[0] = (TH1D*)B8_rootfile->Get(Boron_PDFs_name.c_str());

	}
	
	gRandom = new TRandom3(0);
	gRandom->SetSeed(0);

	double nu_energy;

	if (typenu == "B8") {
		nu_energy = B8_PDFs[0]->GetBinLowEdge(B8_PDFs[0]->FindLastBinAbove(0.) + 1);
	} else {
		nu_energy = CalculateNuEnergy();
	}

	double max_eEnergy = 2*pow(nu_energy,2)/(m_e+2*nu_energy); //MeV

	cout << "Nu_e survival Probability = " << survival_probability(nu_energy) << endl;
	cout << "Nu_e cross section @ 1 MeV = " << cross_section(1.0,nu_energy) << endl;
	cout << "Nu_x cross section @ 1 MeV = " << cross_section_nuX(1.0,nu_energy) << endl;
	cout << "Total cross section @ 1 MeV = " << total_cross_section(1.0,nu_energy) << endl;

	cout <<"Scintillation photons generated @ 1 MeV = " << PEatMeV << endl;
	cout <<"Number of events generated = " << NEvents << endl;
	cout <<"Maximum energy of the incoming neutrino = " << nu_energy << " MeV" << endl;
	cout <<"Maximum electron energy = "<<max_eEnergy<<" MeV" <<endl<<endl;

	//foutput must be defined before the tree to avoid errors 
	TFile *foutput = new TFile (Output_Rootfile.c_str(), "RECREATE");

	//Making the tree
	TTree *t = new TTree("t","New Tree");
	t->Branch("El_Direction_x", &El_Direction_x_t, "El_Direction_x/D");
	t->Branch("El_Direction_y", &El_Direction_y_t, "El_Direction_y/D");
	t->Branch("El_Direction_z", &El_Direction_z_t, "El_Direction_z/D");
	t->Branch("Ph_x_AtPMT","std::vector<double>", &Ph_x_AtPMT_v, 32000, 99);
	t->Branch("Ph_y_AtPMT","std::vector<double>", &Ph_y_AtPMT_v, 32000, 99);
	t->Branch("Ph_z_AtPMT","std::vector<double>", &Ph_z_AtPMT_v, 32000, 99);
	t->Branch("TravelledDistance","std::vector<double>", &TravelledDistance_v, 32000, 99);
	t->Branch("Ph_theta_AtPMT","std::vector<double>", &Ph_theta_AtPMT_v, 32000,99);
	t->Branch("Ph_phi_AtPMT","std::vector<double>", &Ph_phi_AtPMT_v,32000, 99);
	t->Branch("Closest_PMT","std::vector<int>", &Closest_PMT_v, 32000,99);
	t->Branch("Start_Time","std::vector<double>", &Start_Time_v, 32000,99);
	t->Branch("Electron_Energy", &Electron_Energy_t, "Electron_Energy/D");
	t->Branch("Int_Vertex_x", &Int_Vertex_x_t, "Int_Vertex_x/D");
	t->Branch("Int_Vertex_y", &Int_Vertex_y_t, "Int_Vertex_y/D");
	t->Branch("Int_Vertex_z", &Int_Vertex_z_t, "Int_Vertex_z/D");
	t->Branch("Neutrino_Energy", &Neutrino_Energy_t, "Neutrino_Energy/D");
	t->Branch("Type","std::vector<int>", &Type_v, 32000,99);
	t->Branch("Cherenkov_angle",&Cher_angle_t,"Cherenkov_angle/D");
	t->Branch("Electron_angle",&Elec_angle_t,"Electron_angle/D");

	t->Branch("Neutrino_theta",&Neutrino_theta_t,"Neutrino_theta/D");
	t->Branch("Neutrino_phi",&Neutrino_phi_t,"Neutrino_phi/D");
	t->Branch("Solar_theta",&Solar_theta_t,"Solar_theta/D");
	t->Branch("Solar_phi",&Solar_phi_t,"Solar_phi/D");

	t->Branch("Hit","std::vector<bool>" ,&Hit_v, 32000,99);


	// LOAD PMTs position
	ifstream ReadPMTPosition;
	ReadPMTPosition.open(PMTPositions.c_str());
	double blank;
	int Index;
	double x_PMT,y_PMT,z_PMT, r_PMT, theta_PMT, phi_PMT;

	for(int PMT=0;PMT<PMTNumber;PMT++){		
		ReadPMTPosition >> Index;
		ReadPMTPosition >> x_PMT;
		ReadPMTPosition >> y_PMT;
		ReadPMTPosition >> z_PMT;
		ReadPMTPosition >> blank >> blank;
		CartesianToSpherical(r_PMT, theta_PMT, phi_PMT,x_PMT,y_PMT,z_PMT);
		PMT_Position_Spherical.push_back({r_PMT,theta_PMT,phi_PMT});			
		//cout << PMT << "   " <<r_PMT << "  " << theta_PMT << "   " << phi_PMT << endl; 
	}	

	double past_phi = PMT_Position_Spherical[0][2];
	std::vector<tuple<int,double,double,double>> Filler;

	for (int i=0; i<PMTNumber; i++) {

		Index = i;
		r_PMT = PMT_Position_Spherical[i][0];
		theta_PMT = PMT_Position_Spherical[i][1];
		phi_PMT = PMT_Position_Spherical[i][2];

		if (abs(phi_PMT - past_phi) > 1e-5) {
			PMT_Position_Fast.push_back(Filler);
			Filler.clear();
			Filler.push_back(make_tuple(Index,r_PMT,theta_PMT,phi_PMT));
			
		} else {
			Filler.push_back(make_tuple(Index,r_PMT,theta_PMT,phi_PMT));
		}
		past_phi = phi_PMT;
	}

	PMT_Position_Fast.push_back(Filler);
	Filler.clear();

/* If you want to print out the PMT positions and PMT rings

	cout<< "SIZE 0 =      " << PMT_Position_Fast[0].size() << endl; 
	cout<< "SIZE 1 =      " << PMT_Position_Fast[1].size() << endl; 
	cout<< "SIZE 2 =      " << PMT_Position_Fast[2].size() << endl; 
	cout << "SIZE 88 =     " << PMT_Position_Fast[88].size() << endl; 

	
	for (int j=0; j<PMT_Position_Fast.size(); j++) {

		cout << "Block number" << j << ":   " << endl;
		
		for (int itr = 0; itr <  PMT_Position_Fast[j].size(); itr++ ) {
			
			cout << "       " << get<0>(PMT_Position_Fast[j][itr]) << "   " << get<1>(PMT_Position_Fast[j][itr]) << "  " << get<2>(PMT_Position_Fast[j][itr]) << "   " << get<3>(PMT_Position_Fast[j][itr]) << endl;
		}
	
	}
	
	//exit(0);
*/
	ifstream ReadSolarPosition;
	ReadSolarPosition.open(SolarPositions.c_str());
	
	foutput->cd();

	int SeenTotalPhotons = 0;

	for (int i=0; i<NEvents; i++) {
		SeenTotalPhotons += GenerateEvents(t, RandomIntVertex, i, ReadSolarPosition);
		if (NEvents > 10) {  //to avoid floating point exceptions for NEvents < 10
			if (i % 100 == 0 && i != 0) { // check if the index is a multiple of tenth
			std::cout << i << "-th Event ; " << round ( (double)i / (double)NEvents * 10000 ) / 100 << "% of events simulated \n";
    		}
		}
		
	}

	t->Write();
		
	foutput->Close();

	rootfile -> Close();

	ReadSolarPosition.close();

	cout << endl << "Geometric coverage = " << double(SeenTotalPhotons)/double(TotalPhotons) <<endl;
	cout << "#############" << endl;
	
	return 0;

}

int main(int argc, char** argv) {
        string macro = argv[0];

        if(argc!=3) {
                cout << "\n     USAGE:  "<< macro << " Configuration_Text  Output_Rootfile \n" << endl;
                return 1;
        }

        string Configuration_Text = argv[1];
        string Output_Rootfile = argv[2];

		Directionality_ToyMC(Configuration_Text, Output_Rootfile);
        
        return 0;
}
