#ifndef FLIGHTSIMULATOR_ALPHA_TRIM_H
#define FLIGHTSIMULATOR_ALPHA_TRIM_H

int calcolo_alpha_trim(double alpha[], double m, double S, double rho, double* V, double *gamma, double css[][6], double czd[][7],
                        double cuf[][6],double cmd[][7], double cum[][6], double* alpha_trim, double* delta_e_trim, double* Vmin) ;

#endif // !FLIGHTSIMULATOR_ALPHA_TRIM_H
