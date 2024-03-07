#include <stdio.h>
#include <stdlib.h>
#include <unistd.h>
#include <math.h>
#include "../Interpolazione/Interpolazione.h"

void Routh( double m, double S, double rho, double alpha, double V) {
    double **C = MatrixCoeff(10, 7);
    double chord = 0.295856197079013;
    double Jy = 0.01192188547073001;
    double Jy_ad = 8 * m / (rho*S*chord*chord*chord);
    double b = 0.872000000000000;
    double g=9.81;
    double CD0=0.0235;
    double lambda = b * b / S;
    double Cla1 = 3.25; //rad-1
    double mi = 2 * m / (rho * S * chord);
    double Cde, Ctu, Cla, Cmq, Cma1, Cwe, Cma,Cle, k, Cda, A1, B1, C1, D1, E1, delta;
    const double pi = atan(1) * 4;
    Cma1 = 0;
    Cwe=(2*m*g)/(rho*V*V*S);
    Cde=(CD0+k*Cwe*Cwe);
    Ctu = -3*Cde;
    Cla = C[1][0] * sin(alpha) - C[3][0] * cos(alpha);
    Cmq = C[5][3];
    Cle = Cwe;
    k = 0.0382;
    Cda = 2 * k * Cla * Cde;
    Cma = C[5][0];
    A1 = 2 * mi * Jy_ad * (2 * mi + Cla1);
    B1 = 2 * mi * Jy_ad * (Cla + Cde - Ctu) - Jy_ad * Ctu * Cla1 - 2 * mi * Cmq * Cla1 - 4 * mi * mi * (Cmq + Cma1);
    C1 = 2 * mi * (Cmq * (Ctu - Cla - Cde) - 2 * mi * Cma + Cma1 * Ctu) +
         Jy_ad * (2 * Cwe * (Cwe - Cda) + Ctu * Cla + Cde * Cla) + Cmq * Cla1 * Ctu;
    D1 = -2 * Cwe * Cwe * Cma1 + 2 * mi * Ctu * Cma + Ctu * Cmq * Cla - 2 * Cwe * Cmq * (Cle - Cda) +
         2 * Cde * Cmq * Ctu;
    E1 = -2 * Cwe * Cwe * Cma;
    delta = B1 * C1 * D1 - A1 * D1 * D1 - B1 * B1 * E1;

    printf("STABILITA':\n\n");

    if (E1 > 0) {
        printf("Il velivolo e' staticamente stabile.\n");
    } else {
        printf("Il velivolo e' staticamente instabile.simulazione interrotta\n");
        exit(1);
    }

    if (B1 > 0 && D1 > 0 && E1 > 0) {
        printf("Il velivolo e' dinamicamente stabile\n");
    } else {
        printf("Il velivolo e' dinamicamente instabile, simulazione interrotta.\n");
        exit(1);
    }

    printf("--------------------------------------------------------------------------------");

}




