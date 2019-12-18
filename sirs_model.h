#include <stdio.h>
#include <iostream>
#include <math.h>
#include <fstream>
#include <armadillo>

using namespace std;
using namespace arma;

void SIRS(double, double, double, double, double*, double, int);
void SIRfunc(double*, double*, double, double, double, int);
void rk4func(double*, double*, double, double, double, double, int);
