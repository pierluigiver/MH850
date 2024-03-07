#ifndef FLIGHTSIMULATOR_INTERPOLAZIONE_H
#define FLIGHTSIMULATOR_INTERPOLAZIONE_H

void interp_alpha( double AoA, double alpha_des, double *Cp, double Cs);

void interpol_coeffs(int*index, int nAoA, double*ALPHA , double alpha_tr, double css[][6],double cxd[][7], double cyd[][6], double czd[][7],
                     double cld[][6], double cmd[][7],double cnd[][6],double cuf[][6],double cum[][6],double crd[][6]);

double **MatrixCoeff(int r,int c);

#endif //FLIGHTSIMULATOR_INTERPOLAZIONE_H
