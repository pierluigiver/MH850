#include <stdio.h>
//#include <stdlib.h>
#include <math.h>
#include "../Interpolazione/Interpolazione.h"

void long_modi(double m, double rho, double S, double V, double alpha_trim){

    //Costanti
    const double pi = atan(1) * 4; // pi greco
    const double g=9.81; //costante di gravità [m/s^2]
    // const double k = 0.0382;
    const double k = 0.1027;
    const double CD0 = 0.0235; //coefficiente di resistenza (file propeller)
    const double Iy = 0.01192188547073001;
    //const double CLa = 3.984900; //pendenza della curva CL-alpha [rad^-1]
    double **C = MatrixCoeff(10,10);
    const double c=0.2958;


    //Parametri
    double CWe;
    //double **CDe;
    double CXe;
    double CMa;
    double CMq;
    double CMap;

    //Variabili
    double mu;
    double omega_PH; //velocita' angolare modo di Fugoide
    double omega_SP; //velocita' angolare modo di Corto Periodo
    double smorz_PH; //smorzamento modo di Fugoide
    double smorz_SP; //smorzamento modo di Corto Periodo
    double T_PH; //periodo modo di Fugoide
    double T_SP; //periodo modo di Corto Periodo
    double t12_PH; //tempo di dimezzamento modo di Fugoide
    double t12_SP; //tempo di dimezzamento modo di Corto Periodo
    double Iy_ad; //momento di inerzia lungo y adimensionale


    CWe=(2*m*g)/(rho*(V)*(V)*S);
    //CXe=-(CD0+k*CWe*CWe);
    CXe=-0.0127;
    CMa=C[5][0];
    CMq=C[5][3];
    CMap=C[5][1];

    alpha_trim=alpha_trim*pi/180;
    double CLa=C[1][0]*sin(alpha_trim)-C[3][0]*cos(alpha_trim);

    //Modo di Fugoide (lungo periodo)
    mu=(2*m)/(rho*S*c);
    omega_PH=(V*(CWe)*2)/(sqrt(2)*mu*c);
    smorz_PH=(-CXe*3)/(2*sqrt(2)*(CWe));
    T_PH=(2*pi)/(omega_PH*sqrt(fabs((smorz_PH*smorz_PH)-1)));
    t12_PH=0.6931/(smorz_PH*omega_PH);

    //Modo di Corto Periodo
    Iy_ad=(8*Iy)/(rho*S*c*c*c);
    omega_SP=sqrt(-CMa/Iy_ad)*2*(V)/c;
    smorz_SP=(Iy_ad*CLa-2*mu*(CMq+CMap))/(2* sqrt(-2*mu*Iy_ad*(2*mu*(CMa)+(CMq)*CLa)));
    T_SP=(2*pi)/(omega_SP*sqrt(fabs((smorz_SP*smorz_SP)-1)));
    t12_SP=0.6931/(smorz_SP*omega_SP);


    printf("\nMODI LONGITUDINALI:\n");

    printf("\nMODO DI FUGOIDE: \n");
    printf("Pulsazione: %lf\n", omega_PH);
    printf("Smorzamento: %lf\n", smorz_PH);
    printf("Periodo: %lf\n", T_PH);
    printf("Tempo di dimezzamento: %lf\n", t12_PH);
    printf("\n");
    printf("MODO DI CORTO PERIODO: \n");
    printf("Pulsazione: %lf\n", omega_SP);
    printf("Smorzamento: %lf\n", smorz_SP);
    printf("Periodo: %lf\n", T_SP);
    printf("Tempo di dimezzamento: %lf", t12_SP);

    printf("\n--------------------------------------------------------------------------------");


    
}
