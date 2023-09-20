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

#ifndef ROOT_TH1D
#endif

#define PI 3.141592653589793
#include "TMath.h"
#include "Math/Util.h"

int PMTNumber;
double PMTRadius = 0.25;
double JUNORadius = 19.0;
double FV;
int LY;
double ChScRatio;
int NEvents;
int TotalPhotons = 0;
TH1D** Time_PDFs = new TH1D*[2]; //container for the time PDFs, pos 0 for Cherenkov, pos 1 for scintillation

//Useful values

double n = 1.5 ; //refraction index
double c = 299792458 ; // m/s
double m_e = 0.51099895; //MeV    electron mass
double Be7_energy = 0.876; //MeV    energy of a 7Be neutrino
double pep_energy = 1.44; //MeV  energy of a pep neutrino
double nu_energy;
//double Event_Energy = 0.5; //MeV
double G_F = 1.1663787*pow(10,-11); //MeV^-2  Fermi constant
double sin2_thetaW = 0.23121; //Weinberg angle

double max_eEnergy; //  maximum electron energy from a neutrino scattering
double min_eEnergy; // minimum deposited energy from a neutrino scattering

double RefractionIndex = 1.5;

// VARIABLES NO LONGER USEFUL
double TravelledDistance = 19.0;
double AbsorptionLength = 80.;
double QE = 0.3;
double AbsorptionProbability = (1-TMath::Exp(-TravelledDistance/AbsorptionLength));
double ReemissionProbability = 0.75;
double SurvivingProbability = (1-AbsorptionProbability*(1-ReemissionProbability))*QE;

int PEatMeV; //photoelectron emitted @ 1 MeV

double El_Direction_x_t,El_Direction_y_t,El_Direction_z_t,phi_t,theta_t,Closest_PMT_t,Start_Time_t,Arr_Time_t,Electron_Energy_t,Neutrino_Energy_t,type_t;
double Ph_x_AtPMT_t,Ph_y_AtPMT_t,Ph_z_AtPMT_t, Ph_r_AtPMT_t, Ph_theta_AtPMT_t, Ph_phi_AtPMT_t, TravelledDistance_t;
double Int_Vertex_x_t,Int_Vertex_y_t,Int_Vertex_z_t;
bool Type_t, IsFirst_t;
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
std::vector <double> *Arr_Time_v;
std::vector <bool> *Type_v;
std::vector <double> *Min_Distance_v;
std::vector <bool> *Hit_v;

double Min_Distance_t;
bool Hit_t;
bool fastmode; //does not save the events at more than 5 ns
double TimeCut; //cut photons emitted at times lower than timecut
bool FixedSun, IsBackgrounds;
double Fastmode_cut; //time cut in fastmode

double MaxPMTDistance; //maximum distance in phi between two PMTs

using namespace std;

struct Tuple {
	double x;
	double y;
};

//to keep phi and theta in their intended range
double Pbc_phi (double in) {
    if (in<0) return in + M_PI;
    else if (in>M_PI) return in - M_PI;
    else return in;
}

double Pbc_theta (double in) {
    if (in<0) return in + 2*M_PI;
    else if (in>2*M_PI) return in - 2*M_PI;
    else return in;
}

double GetFloatPrecision(double value, double precision){
    return (floor((value * pow(10, precision) + 0.5)) / pow(10, precision));
}

inline bool exists_test0 (const std::string& name) {
    ifstream f(name.c_str());
    return f.good();
}

double DegToRad (double deg) {
	return deg / 180 * M_PI;
}

double RadToDeg (double rad) {
	return rad * 180 / M_PI;
}

void CartesianToSpherical(double & r,double & theta,double & phi,double x,double y,double z ){
	
	r = sqrt(x*x+y*y+z*z);
	if ( r != 0.0 ){
		phi = TMath::ACos( z / r );
		theta = TMath::ATan2( y, x ) + PI;
		}
	else
	theta = phi = 0.0;
}

void SphericalToCartesian(double & x, double & y, double & z, double r, double theta, double phi){
	if ( r < 0.0 ){
		throw "Negative radius in sphericalToCartesian()";
		}
	x = r * TMath::Sin( phi ) * TMath::Cos( theta );
	y = r * TMath::Sin( phi ) * TMath::Sin( theta );
	z = r * TMath::Cos( phi );
}

double Distance(double x1, double y1, double z1, double x2, double y2, double z2){
	return sqrt( pow(x2-x1,2) + pow(y2-y1,2) + pow(z2-z1,2));
}
	
double DistanceOnASphere(double r, double theta1, double phi1, double theta2, double phi2){
	return TMath::ACos(TMath::Cos(phi1)*TMath::Cos(phi2) + (TMath::Sin(phi1)*TMath::Sin(phi2))*TMath::Cos(theta1-theta2));
}

//calculate 1st order cross section for a neutrino-electron elastic scattering
double cross_section (double T) {
	double gl = 0.5 - sin2_thetaW -1 ;
	double gr = - sin2_thetaW;

	return ((2*(pow(G_F,2))*m_e)/M_PI)*(pow(gl,2) + pow(gr,2)*pow(1-T/nu_energy,2) - gl*gr*m_e*T/pow(nu_energy,2));
}

double cross_section_nuX (double T) {
	double g_L = 0.5 - sin2_thetaW;
	double g_R = -sin2_thetaW;

	return 2*m_e*G_F*G_F/M_PI * (pow(g_L,2) + pow(g_R,2)*pow(1-T/nu_energy,2) - g_L*g_R*m_e/nu_energy*T/nu_energy);
}

double survival_probability () {

	double E_nu = nu_energy*pow(10,6);
	//datas from pep, see Xuefeng technote
	double V_x = 5.13e-12;
	double Delta_V_x2__V_x2 = 0.076;

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
double total_cross_section (double T) {

	return survival_probability()*cross_section(T) + (1-survival_probability())*cross_section_nuX(T);

}

//sample cross_section with an accept-reject method
double CalculateEventEnergy () {
	double EventEnergy, test;
	bool flag = false;

	if (IsBackgrounds == true) {
		return gRandom -> Uniform(min_eEnergy,max_eEnergy);
	}

	while (flag == false) {
		flag = false;
		EventEnergy = gRandom -> Uniform(min_eEnergy,max_eEnergy);
		test = gRandom -> Uniform(0.,total_cross_section(0.)); //the maximum cross section is at T=0
		if (test < total_cross_section(EventEnergy)) {
			flag = true;
		}
	}

	return EventEnergy;
}

double ClosestPMTIndex (double x_Event,double y_Event,double z_Event, vector<vector<double>> & PMT_Position_Spherical) {

	double MinDistance=50000;
	double Distance_Temp;
	int Index=50000;
	int Closest=Index;
	double r_Event, theta_Event, phi_Event;
	double x_PMT,y_PMT,z_PMT, r_PMT, theta_PMT, phi_PMT;

	CartesianToSpherical(r_Event, theta_Event, phi_Event,x_Event, y_Event, z_Event);
		
	for(int PMT=0;PMT<PMTNumber;PMT++){
		r_PMT = PMT_Position_Spherical[PMT][0];
		theta_PMT = PMT_Position_Spherical[PMT][1];
		phi_PMT = PMT_Position_Spherical[PMT][2];
		Distance_Temp = DistanceOnASphere(r_PMT, theta_PMT, phi_PMT,theta_Event, phi_Event);
		//cout << PMT <<  " r " << r_PMT << "  th " <<   theta_PMT << "  phi " << phi_PMT <<  endl;
		
		if(Distance_Temp<MinDistance){
			Closest = PMT;
			MinDistance = Distance_Temp;
		} 
		if (MinDistance*JUNORadius <= PMTRadius) {

			Min_Distance_t = JUNORadius * MinDistance;

			return Closest;
		}

		//Code to speed up the process
		
		if (Distance_Temp > MaxPMTDistance*JUNORadius + 0.25) {
			PMT += 9;
		}
	
		//cout << theta_PMT << "  "  << phi_PMT << "    "  << Index << "   " << Distance_Temp << "  /   " << Closest << "   " << MinDistance << endl;	// DO NOT REMOVE
	}

	Min_Distance_t = JUNORadius * MinDistance;
	return -1;
}

//outputs the position (x,y,z) on the sphere for a photon generated in (j,k,l) parallel to the vector (a,b,c)

void MovePhoton (double & x,double& y , double& z, double j, double k, double l, double a, double b , double c ) {

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

Tuple Generate_Cone (double theta_0, double phi_0, double angle) {

	Tuple out;
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

	out.x = theta_out;
	out.y = phi_out;

	return out;
}

double GenerateScintStartTime () {

	double StartTime;

	double Check = gRandom -> Uniform(0,1);

	if (Check < 0.707) {
		StartTime = gRandom -> Exp(4.6*pow(10,-9));
	} else if (Check < 0.912) {
		StartTime = gRandom -> Exp(15.1*pow(10,-9));
	} else if (Check < 0.972) {
		StartTime = gRandom -> Exp(76.1*pow(10,-9));
	} else {
		StartTime = gRandom -> Exp(397*pow(10,-9));
	}

	return StartTime;
}

int CheckHit (int SeenPhotons) {

	if (Closest_PMT_t != -1) {
		Hit_t = 1;
		SeenPhotons++;

	} else {Hit_t = 0;}


	return SeenPhotons;
}

//Generator for all the photons

int GeneratePhotons (TTree* t, vector<vector<double>> PMT_Position_Spherical, bool RandomPos, int NEvent, ifstream& ReadSolarPosition) {

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
	double Event_Energy = CalculateEventEnergy();

	double theta_e = acos((1+m_e/nu_energy)*pow(Event_Energy/(Event_Energy+2*m_e),0.5)); //angle between the solar-nu and the electron scattered (assuming 7Be-nu)
	double beta_el = pow(1-(pow(m_e/(Event_Energy+m_e),2)),0.5) ; //beta of the electron generated
	double theta_Cher = acos(1/(beta_el*n)); //Cherenkov angle

	//temp_out << theta_e << "  " << theta_Cher << endl;

	int Photons = gRandom -> Poisson(PEatMeV*Event_Energy);
	int CherenkovPhotons = gRandom -> Poisson(ChScRatio*Photons);

	//Correction for the decimal part of CherenkovPhotons


	if (NEvent < 10) {
		cout << "Event #" << NEvent << " :  " << Photons << " scintillation and " << CherenkovPhotons << " Cherenkov photons emitted";
		cout << " (Electron energy = " << Event_Energy << " MeV)" << endl;  
	} else if (NEvent == 10) {
		cout<<endl;
	}

	bool IsFirstFlag = true;

	double Provv_StartTime;

	//cout<<Photons<<"  "<<CherenkovPhotons<<"  "<<Event_Energy<<endl;
	//cout<<max_eEnergy<<endl;

	//Reading the solar positions
	double solar_nu_theta, solar_nu_phi;
	double solar_theta, solar_phi; //Direction opposite to the sun's position

	if (FixedSun) {
		
		solar_nu_theta = 0.;
		solar_nu_phi = 0.;
		solar_theta = 0.;
		solar_phi = 0.;

	} else {
		double blank;
		double solar_az,solar_alt;
		double theta_a,phi_a;
		
		if (ReadSolarPosition.peek() != EOF) ReadSolarPosition >> solar_az >> solar_alt >> theta_a >> phi_a ;
		else {
			cerr << "ERROR: Solar positions terminated, you generated more than 10yr of events" << endl;
			exit(1);
		}

		solar_phi = Pbc_phi(theta_a);
		solar_theta = Pbc_theta(phi_a); //TO CHECK, DEPENDS FROM THE SYSTEM OF COORDINATES OF THE PMTs

		if (IsBackgrounds) {
			solar_nu_theta = gRandom -> TRandom::Uniform(2*PI);
			solar_nu_phi = TMath::ACos(-1.+2.*gRandom->TRandom::Uniform(0,1));
		} else {

			solar_nu_phi = Pbc_phi ( M_PI - solar_phi);
			solar_nu_theta = Pbc_theta ( M_PI + solar_theta);

		}

		
	}

	

/*
	cout << "Altitude " << solar_alt << " Azimut " << solar_az << endl;
	cout << "Solar phi " << RadToDeg(solar_phi) << " Solar theta " << RadToDeg(solar_theta) << endl;
	cout << "Nu phi " << RadToDeg(solar_nu_phi) << " Nu theta " << RadToDeg(solar_nu_theta) << endl;

	if (NEvent == 10 ) exit(1);
*/

	//Generate Electron direction
	Tuple a;
	Tuple b;

	a = Generate_Cone(solar_nu_theta,solar_nu_phi,theta_e);

	double El_theta = a.x;
	double El_phi = a.y;

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
	Arr_Time_v = new vector<double> ();
	Type_v = new vector<bool> ();
	Min_Distance_v = new vector<double> ();
	Hit_v = new vector<bool> ();

	//Generate SCINTILLATION Photons

	for(int iPh=0; iPh<Photons; iPh++){

		//Generate unit vector over a sphere

		do {
			Provv_StartTime = Time_PDFs[1] -> GetRandom();
		} while (Provv_StartTime < TimeCut);
		
		Start_Time_t = Provv_StartTime;

		if (fastmode && Start_Time_t > Fastmode_cut) {
			continue;
		}

		theta_vers = gRandom->TRandom::Uniform(2*PI);
		phi_vers = TMath::ACos(-1.+2.*gRandom->TRandom::Uniform(0,1));
		SphericalToCartesian(ph_Direction_x,ph_Direction_y,ph_Direction_z,1,theta_vers,phi_vers);

		MovePhoton(Ph_x_AtPMT_t,Ph_y_AtPMT_t,Ph_z_AtPMT_t,x_Int,y_Int,z_Int,ph_Direction_x,ph_Direction_y,ph_Direction_z);

		CartesianToSpherical(trash,theta_t,phi_t,Ph_x_AtPMT_t,Ph_y_AtPMT_t,Ph_z_AtPMT_t);

		Ph_r_AtPMT_t = JUNORadius; // NOT REALLY TRUE
		Ph_theta_AtPMT_t = theta_t; // TRUE ONLY WITHOUT ABSORPTION
		Ph_phi_AtPMT_t = phi_t; // TRUE ONLY WITHOUT ABSORPTION

		Closest_PMT_t = ClosestPMTIndex(Ph_x_AtPMT_t,Ph_y_AtPMT_t,Ph_z_AtPMT_t,PMT_Position_Spherical);

		TravelledDistance_t = Distance(x_Int,y_Int,z_Int,Ph_x_AtPMT_t,Ph_y_AtPMT_t,Ph_z_AtPMT_t);
		Type_t = 0; 
		
		Arr_Time_t = Start_Time_t + (TravelledDistance_t/c*n) * pow(10,9);

		
		//cout << iPh << "\t" << xx_at_PMTs << "\t" << yy_at_PMTs << "\t" << zz_at_PMTs << endl;
			
		//cout << "PHOTON " << iPh << " : CLOSEST INDEX IS = " << IndexExample << endl;

		//Generate the photon only if it hits the PMT
		SeenPhotons = CheckHit(SeenPhotons);

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
		Arr_Time_v -> push_back(Arr_Time_t);
		Type_v -> push_back(Type_t);
		Min_Distance_v -> push_back(Min_Distance_t);
		Hit_v -> push_back(Hit_t);

		//cout << i << "    " << iPh % (Photons/10) << endl;
						
	}

	//Generate CHERENKOV Photons


	for (int iPh=0; iPh<CherenkovPhotons; iPh++) {

		do {
			Provv_StartTime = Time_PDFs[0] -> GetRandom();
		} while (Provv_StartTime < TimeCut);

		Start_Time_t = Provv_StartTime;

		if (fastmode && Start_Time_t > Fastmode_cut) {
			continue;
		}

		b = Generate_Cone(El_theta,El_phi,theta_Cher);

		theta_vers = b.x;
		phi_vers = b.y;
		SphericalToCartesian(ph_Direction_x,ph_Direction_y,ph_Direction_z,1,theta_vers,phi_vers);

		MovePhoton(Ph_x_AtPMT_t,Ph_y_AtPMT_t,Ph_z_AtPMT_t,x_Int,y_Int,z_Int,ph_Direction_x,ph_Direction_y,ph_Direction_z);

		CartesianToSpherical(trash,theta_t,phi_t,Ph_x_AtPMT_t,Ph_y_AtPMT_t,Ph_z_AtPMT_t);

		Ph_r_AtPMT_t = JUNORadius; // NOT REALLY TRUE
		Ph_theta_AtPMT_t = theta_t; // TRUE ONLY WITHOUT ABSORPTION
		Ph_phi_AtPMT_t = phi_t; // TRUE ONLY WITHOUT ABSORPTION

		Cher_angle_t = theta_Cher;
		Elec_angle_t = theta_e;

		TravelledDistance_t = Distance(x_Int,y_Int,z_Int,Ph_x_AtPMT_t,Ph_y_AtPMT_t,Ph_z_AtPMT_t);

		Closest_PMT_t = ClosestPMTIndex(Ph_x_AtPMT_t,Ph_y_AtPMT_t,Ph_z_AtPMT_t,PMT_Position_Spherical);

		Type_t = 1;

		Arr_Time_t = Start_Time_t + (TravelledDistance_t/c*n)* pow(10,9);


		//Generate the photon only if it hits the PMT
		SeenPhotons = CheckHit(SeenPhotons);

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
		Arr_Time_v -> push_back(Arr_Time_t);
		Type_v -> push_back(Type_t);
		Min_Distance_v -> push_back(Min_Distance_t);
		Hit_v -> push_back(Hit_t);

		
	}

	t -> Fill();

	Ph_x_AtPMT_v -> clear();
	Ph_y_AtPMT_v -> clear();
	Ph_z_AtPMT_v -> clear();
	TravelledDistance_v -> clear();
	Ph_theta_AtPMT_v -> clear();
	Ph_phi_AtPMT_v -> clear();
	Closest_PMT_v -> clear();
	Start_Time_v -> clear();
	Arr_Time_v -> clear();
	Type_v -> clear();
	Min_Distance_v -> clear();
	Hit_v -> clear();
	
	return SeenPhotons;
	
}


double Directionality_ToyMC(string Configuration_Text, string Output_Rootfile) {

	ifstream file(Configuration_Text);
	//ofstream temp_out("Theta_e.txt");
	vector<string> col1;
	vector<string> col2;
	string line;
	string cher_times, scint_times, typenu;
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

	// quite rough; to be changed
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
	iss17 >> Fastmode_cut;

	//cout << "PMT Number " << PMTNumber << endl;

	// ### End parsing	

	TFile *rootfile = new TFile(origin_rootfile.c_str());
	
	Time_PDFs[0] = (TH1D*)rootfile->Get(cher_times.c_str());
	Time_PDFs[1] = (TH1D*)rootfile->Get(scint_times.c_str());

	gRandom = new TRandom3(0);
	gRandom->SetSeed(0);

	if (typenu == "Be7") {
		nu_energy = Be7_energy;
	} else if (typenu == "pep") {
		nu_energy = pep_energy;
	} else {
		cout << "ERROR : invalid nu type in "<< Configuration_Text.c_str() << endl;
		exit(1);
	}

	max_eEnergy = 2*pow(nu_energy,2)/(m_e+2*nu_energy); //MeV

	cout << "Nu_e survival Probability = " << survival_probability() << endl;
	cout << "Nu_e cross section @ 1 MeV = " << cross_section(1.0) << endl;
	cout << "Nu_x cross section @ 1 MeV = " << cross_section_nuX(1.0) << endl;
	cout << "Total cross section @ 1 MeV = " << total_cross_section(1.0) << endl;

	std::vector<vector<double>> PMT_Position_Spherical;		
		
	std::cout << "AbsorptionProbability = " << AbsorptionProbability << endl;
	cout << "SurvivingProbability = " << SurvivingProbability << endl;
	//cout<<"Neutrino-electron angle = "<<theta_e*180./M_PI<<" deg"<<endl;
	//cout<<"Cherenkov angle = "<<theta_Cher*180./M_PI<<" deg"<<endl;

	cout <<"Scintillation photons generated @ 1 MeV = " << int(PEatMeV) << endl;
	cout <<"Number of events generated = " << NEvents << endl;
	cout <<"Energy of the incoming neutrino = " << nu_energy << " MeV" << endl;
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
	t->Branch("Arr_Time","std::vector<double>", &Arr_Time_v, 32000,99);
	t->Branch("Electron_Energy", &Electron_Energy_t, "Electron_Energy/D");
	t->Branch("Int_Vertex_x", &Int_Vertex_x_t, "Int_Vertex_x/D");
	t->Branch("Int_Vertex_y", &Int_Vertex_y_t, "Int_Vertex_y/D");
	t->Branch("Int_Vertex_z", &Int_Vertex_z_t, "Int_Vertex_z/D");
	t->Branch("Neutrino_Energy", &Neutrino_Energy_t, "Neutrino_Energy/D");
	t->Branch("Type","std::vector<bool>", &Type_v, 32000,99);
	t->Branch("Cherenkov_angle",&Cher_angle_t,"Cherenkov_angle/D");
	t->Branch("Electron_angle",&Elec_angle_t,"Electron_angle/D");

	t->Branch("Neutrino_theta",&Neutrino_theta_t,"Neutrino_theta/D");
	t->Branch("Neutrino_phi",&Neutrino_phi_t,"Neutrino_phi/D");
	t->Branch("Solar_theta",&Solar_theta_t,"Solar_theta/D");
	t->Branch("Solar_phi",&Solar_phi_t,"Solar_phi/D");

	t->Branch("Min_Distance","std::vector<double>", &Min_Distance_v, 32000,99);
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

	//exit(1);

	// #################### PMT_Position_Sperical is already sorted in phi ####################################
	//sorting PMT_Position_Spherical in phi
	//std::sort(PMT_Position_Spherical.begin(),PMT_Position_Spherical.end(),[] (const std::vector<int>& a, const std::vector<int>& b) {return a[2] < b[2];});

	/*ofstream WriteTry;
	WriteTry.open("Prova.txt");
	cout << "##########################" << endl;
	cout << "Sorted PMT_Position_Spherical" << endl;
	for(int PMT=0;PMT<PMTNumber;PMT++){	
		WriteTry<<PMT_Position_Spherical[PMT][2]<<endl;
	}
	cout << "##########################" << endl;*/
	

	//Find MaxPMTDistance
	MaxPMTDistance = 0;
	double ProvPMTDistance;

	for(int PMT=0;PMT<PMTNumber-1;PMT++){		
		ProvPMTDistance = PMT_Position_Spherical[PMT+1][2] - PMT_Position_Spherical[PMT][2];
		if (ProvPMTDistance > MaxPMTDistance) {
			MaxPMTDistance = ProvPMTDistance;
		}	
	}	
	//cout << "MaxPMTDistance = " << MaxPMTDistance<<endl;

	ifstream ReadSolarPosition;
	ReadSolarPosition.open(SolarPositions.c_str());

	
	foutput->cd();

	int SeenPhotons = 0;

	for (int i=0; i<NEvents; i++) {
		SeenPhotons += GeneratePhotons(t, PMT_Position_Spherical, RandomIntVertex, i, ReadSolarPosition);
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

	cout << endl << "Geometric coverage = " << double(SeenPhotons)/double(TotalPhotons) <<endl;
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
