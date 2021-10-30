#include <string>
#include <stdio.h>
#include <cmath>
#include <complex>
#include <fstream>

using namespace std;

extern double pi;
extern double val_e;
extern double val_epsilonO;
extern double val_muO;
extern double val_freq;
extern double val_omega;
extern double val_c;
extern double val_lambdaO;
extern double val_kO;
extern double val_etaO;

void genetare_constants(){
    //
    //
    // genetare_constants
    // Setting the physical constants
    val_omega = 2*pi*val_freq;
    val_c = 1/sqrt(val_epsilonO*val_muO);
    val_lambdaO = val_c/val_freq;
    val_kO = 2*pi/val_lambdaO;
    val_etaO = sqrt(val_muO/val_epsilonO);

    printf("\n");
    printf("Frequency: %f GHz\n", val_freq/1e9);
    printf("Omega: %f rad Hz \n", val_omega);
    printf("Epsilon: %1.4e F/m \n", val_epsilonO);
    printf("Mu: %1.4e H/m \n", val_muO);
    printf("Speed of Light: %1.4e m/s \n", val_c);
    printf("Wavelength %f m \n", val_lambdaO);
    printf("Wave number: %f rad/m \n", val_kO);
    printf("Wave Impedance: %f ohm \n", val_etaO);
    printf("\n");

}
