#include <stdio.h>
#include <math.h> // Per le funzioni trigonometriche
#include <stdlib.h> 
#include "../Interpolazione/Interpolazione.h"
#include "../Calcolo_spinta_RPM//Calcolo_spinta_RPM.h"

int Condizioni_trim(double trim_cond[],double*trim_control, double alpha_trim, double  V,
                     double gamma, double h, double delta_e_trim, double m, double S,double rho, double RPM_min,
                     double RPM_max, double datipropeller[],double BA[]){

    double p_trim = 0,q_trim = 0,r_trim = 0,phi_trim = 0,psi_trim = 0,v_trim = 0, Da_tr=0;
    int i = 0;
    double res = 1;
    int flag1 = 0;
    double P_max=160; //[W]
    double Delta_RPM = 100;
    double g = 9.81;
    double inc = 1.5;
    const double pi = atan(1) * 4; // pi greco
    double d2r = pi/180;
    double u_trim, w_trim, theta_trim, h_trim, Dth_tr;
    double Cx_tot, T_trim, RPM_trim;
    int flag = 0;
    double **C = MatrixCoeff(10,10);

    double CXss = C[0][0];
    double CXa = C[1][0];
    double CXde = C[7][0];
    double Czss=C[0][2];
    double CLmax;
    double Vmin;

    u_trim = V * cos(alpha_trim);
    w_trim = V * sin(alpha_trim);
    theta_trim = alpha_trim * 180 / pi + gamma;
    h_trim = h;

    printf("Il vettore di stato di trim e':\n");
    printf("u_trim:\t%lf\t[m/s]\n", u_trim);
    printf("v_trim:\t%lf\t[m/s]\n", v_trim);
    printf("w_trim:\t%lf\t[m/s]\n", w_trim);
    printf("p_trim:\t%lf\t[rad/s]\n", p_trim);
    printf("q_trim:\t%lf\t[rad/s]\n", q_trim);
    printf("r_trim:\t%lf\t[rad/s]\n", r_trim);
    printf("phi_trim:\t%lf\t[deg]\n", phi_trim);
    printf("theta_trim:\t%lf\t[deg]\n", theta_trim);
    printf("psi_trim:\t%lf\t[deg]\n", psi_trim);
    printf("h_trim:\t%lf\t[m]\n", h);
    printf("x_trim:\t%lf\t[m]\n", 0.00);
    printf("y_trim:\t%lf\t[m]\n", 0.00);


    Cx_tot = CXss + CXa * alpha_trim + CXde * delta_e_trim;
    T_trim = m * g * sin(theta_trim * pi / 180) - 0.5 * rho * Cx_tot * S * V * V;

    trim_cond[0]=u_trim; trim_cond[1]=v_trim; trim_cond[2]=w_trim; trim_cond[3]=p_trim; trim_cond[4]=q_trim; trim_cond[5]=r_trim;
    trim_cond[6]=phi_trim; trim_cond[7]=theta_trim; trim_cond[8]=psi_trim; trim_cond[9]=h_trim; trim_cond[10]=0; trim_cond[11]=0;

    RPM_trim = propel_trim(T_trim, alpha_trim, delta_e_trim, V, rho);
    Dth_tr = ((1-0.1)/(30000-3600))*(RPM_trim-3600)+0.1;

    trim_control[0]=delta_e_trim; trim_control[1]=Da_tr; trim_control[2]=Dth_tr;

    printf("\nIl vettore di controllo di trim e':\n");
    printf("De_trim:\t%lf\t[deg]\n", trim_control[0]/d2r);
    printf("Da_trim:\t%lf\t[deg]\n", trim_control[1]/d2r);
    printf("Dth_trim:\t%lf\t[percentuale]\n", trim_control[2]);
    printf("--------------------------------------------------------------------------------\n");


}
