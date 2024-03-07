//
// Created by Pierluigi on 13/03/2023.
//

#ifndef FLIGHTSIMULATOR_ATMOSFERA_H
#define FLIGHTSIMULATOR_ATMOSFERA_H

int AtmosphereChoice (double *press0,double *temp0,double *rho0,double *vsuono0,
                      double *press_h,double *temp_h,double *rho_h,double *vsuono_h, double *CI, int *flagatm);
int AtmosphereCalc (double *CI, double *Pmax_h,double *press0,double *temp0,
                    double *press_h,double *temp_h,double *rho_h, int *flagatm);
void loadCI(double *CI, int* inp);

#endif //FLIGHTSIMULATOR_ATMOSFERA_H


