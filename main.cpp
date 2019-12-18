#include <iostream>
#include <armadillo>

#include "sirs_model.h"
#include "mc_sirs.h"

using namespace std;
using namespace arma;

int main()
{
    double b;
    cout << "b rate = ";
    cin >> b;

    double S = 300;
    double I = 100;
    double R = 0;

    /*
    //task a
    double t0, tN, h;
    t0 = 0;
    h = 0.1;
    tN = 7;

    //tids array
    int n = static_cast<int>((tN-t0)/h);
    double tid[n];
    tid[0] = t0;
    for(int i=1; i< n; i++) {
        tid[i] = tid[i-1] + h;
    }

    SIRS(b, S, R, I, tid, h, n);
    */


    //task b
    double t[3];
    t[0] = 4.0/(4.0*400.0);
    t[1] = 1.0/(b*400.0);
    t[2] = 1.0/(0.5*400.0);
    double dt = 0;
    if(t[0] <= t[1]) {
        if(t[0] < t[2]) {
            dt = t[0];
        }
        else if (t[0] > t[2]) {
            dt = t[2];
        }
    }
    else if (t[0] > t[1]) {
        if(t[1] < t[2]) {
            dt = t[1];
        }
        else if (t[1] > t[2]) {
            dt = t[2];
        }
    }
    //finn den minste verdien
    int n = 7/dt;
    double tid[n]; tid[0] = 0;
    for(int i=1; i<n; i++) {
        tid[i] = tid[i-1] + dt;
    }
    mc_SIRS(b, S, R, I, tid, dt, n);


}


//test
//sjekk at N er den samme hele tiden
