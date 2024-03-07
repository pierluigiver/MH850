
#ifndef FLIGHTSIMULATOR_INTEGRAZIONE_H
#define FLIGHTSIMULATOR_INTEGRAZIONE_H
#ifndef SIMULATOR_INTEGRATION_H
#define SIMULATOR_INTEGRATION_H

void lat_dir_controller(int i, double dt, double psi_ref, double psi, double phi,double*I_psi, double*I_phi, double*D_psi, double*D_phi,
                        double*e_phip, double*e_psip, double*da_c);

void long_controller(int i, double dt, double V_ref, double h_ref, double theta, double V, double h, double*I_V, double*D_V, double*e_Vp,
                     double*I_th, double*D_th, double*e_thp,double*I_h, double*D_h, double*e_hp, double*dth_c, double*de_c);

void integration(double trim_cond[], double trim_control[],double alpha_trim,double rho, double S, double b, double c, double m, double Jx, double Jy, double Jz, double Jxz,
                 double V_ref, double h_ref, double Vmin, int inp);

void PID(int i, double Kp, double Ki, double Kd, double N, double dt, double e, double *I, double *D, double *ep, double *PID);

void saturation(double umax, double umin, double*u,int*sat_flag,int*flagp);

double search_T(double dth, double*RPM, double*Torque, double*SPINTA);

void ReadPsi(double psi[], double* tf, double* x, double* y);

#endif //SIMULATOR_INTEGRATION_H

#endif //FLIGHTSIMULATOR_INTEGRAZIONE_H
