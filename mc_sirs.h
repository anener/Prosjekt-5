#include <stdio.h>
#include <iostream>
#include <math.h>
#include <fstream>
#include <random>
#include <armadillo>

using namespace std;
using namespace arma;


void mc_rk4func(double*, double*, double, double, double, double, int, int);
void mc_SIRfunc(double*, double*, double, double, double, int, int);
void mc_SIRS(double, double, double, double, double*, double, int);
