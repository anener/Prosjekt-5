#include "sirs_model.h"

void SIRS(double b, double S, double R, double I, double* tid, double h, int n) {
    double y_out[3], y_h[3]; //to lister som skal holde S, I og R verdiene under utregningene
    y_h[0] = S; y_h[1] = I; y_h[2] = R; //setter startverdiene til S, I og R

    double a = 4; double c = 0.5; int N = 400;

    string k = "SIRS_season5_b" + to_string(static_cast<int>(b)) + ".txt"; //lager en tekstfil for aa skrive svarene i
    ofstream fout(k);
    fout.setf(ios::scientific);
    fout.precision(4);

    fout << n << "\n"; //legger inn som forste linje antall tid steg
    fout << tid[0] << " " << y_h[0] << " " << y_h[1] << " " << y_h[2] << " " << y_h[0]+y_h[1]+y_h[2] <<"\n"; //skriver inn verdiene ved t=0

    //en for-loop som looper over en tids periode
    for(int i=0; i<=n; i++) {
        //double A = 3; double omega = M_PI*2.0;
        //double a = A*cos(omega*tid[i]) + a0;
        rk4func(y_h, y_out, h, a, b, c, N); //regner ut diff. ved hjelpf fra rk4
        fout << tid[i+1] << " " << y_out[0] << " " << y_out[1] << " " << y_out[2] << " " << y_out[0]+y_out[1]+y_out[2] <<"\n"; //skriver inn svarene i tekstfilen
        y_h[0] = y_out[0]; y_h[1] = y_out[1]; y_h[2] = y_out[2]; //setter nye start verdier
    }
    fout.close(); //lukker filen

}

void SIRfunc(double* in, double* out, double a, double b, double c, int N) {
    double d = 0; double di = 0; double e = 0; double f = 0;
    out[0] = c*in[2] - (a*in[0]*in[1])/N - d*in[0] + e*N - f;
    out[1] = (a*in[0]*in[1])/N - b*in[1] - d*in[1] - di*in[1];
    out[2] = b*in[1] - c*in[2] - d*in[2] + f;
}

void rk4func(double* y_in, double* y_out, double h, double a, double b, double c, int N) {
    double k1[3], k2[3], k3[3], k4[3]; //lager fire tomme lister som skal holde k1, k2, k3 og k4 verdiene til S, I og R utregningene
    double y_k[3]; //lager en tom liste som skal holde paa midlertidig verdier for S, I og R

    //k1
    SIRfunc(y_in, y_out, a, b, c, N);
    for(int i=0; i<=2; i++) {
        k1[i] = h*y_out[i];
        y_k[i] = y_in[i] + k1[i]*0.5;
    }

    //k2
    SIRfunc(y_k, y_out, a, b, c, N);
    for(int i=0; i<=2; i++) {
        k2[i] = h*y_out[i];
        y_k[i] = y_in[i] + k2[i]*0.5;
    }

    //k3
    SIRfunc(y_k, y_out, a, b, c, N);
    for(int i=0; i<=2; i++) {
        k3[i] = h*y_out[i];
        y_k[i] = y_in[i] + k3[i];
    }

    //k4
    SIRfunc(y_k, y_out, a, b, c, N);
    for(int i=0; i<=2; i++) {
        k4[i] = y_out[i]*h;
    }

    //calculate the new values for S, I and R
    for(int i=0; i<=2; i++) {
        y_out[i] = y_in[i] + (1.0/6.0)*(k1[i] + 2.0*k2[i] + 2.0*k3[i] + k4[i]);
    }

}

