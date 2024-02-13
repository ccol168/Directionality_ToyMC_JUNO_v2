#include <fstream>
#include <string>
#include <iostream>
#include <cstdlib>
#include <sstream>
#include <iomanip>
#include "TMath.h"

#include "Functions.h"

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
		theta = TMath::ATan2( y, x );
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
	
// distance between two points on a sphere
double DistanceOnASphere(double r, double theta1, double phi1, double theta2, double phi2){
	return TMath::ACos(TMath::Cos(phi1)*TMath::Cos(phi2) + (TMath::Sin(phi1)*TMath::Sin(phi2))*TMath::Cos(theta1-theta2));
}

double mod (double x, double y, double z) {
    return sqrt(x*x + y*y + z*z);
};
