#ifndef FLIGHTSIMULATOR_CALCOLO_SPINTA_RPM_H
#define FLIGHTSIMULATOR_CALCOLO_SPINTA_RPM_H

double propel_trim(double T_trim, double alpha_trim, double deltae_trim, double Vel, double rho1);
double propel(double n, double Vel, double rho1);

#endif // !FLIGHTSIMULATOR_CALCOLO_SPINTA_RPM_H
