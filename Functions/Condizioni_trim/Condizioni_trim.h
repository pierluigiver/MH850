#ifndef FLIGHTSIMULATOR_CONDIZIONI_TRIM_H
#define FLIGHTSIMULATOR_CONDIZIONI_TRIM_H

void Condizioni_trim(double trim_cond[],double*trim_control,double alpha_trim, double  V,
                     double gamma, double h, double delta_e_trim, double m, double S,double rho, double RPM_min,
                     double RPM_max, double datipropeller[],double BA[]);
double Spinta(double RPM_tmp, double RPM[], double T[], double C[]);
#endif // !FLIGHTSIMULATOR_CONDIZIONI_TRIM_H