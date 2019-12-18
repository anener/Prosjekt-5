#include "mc_sirs.h"

void mc_SIRS(double b, double S, double R, double I, double* tid, double dt, int n) {
    double y_out[3], y_h[3];
    y_h[0] = S; y_h[1] = I; y_h[2] = R;

    string k = "mcSIRS_death2_b" + to_string(static_cast<int>(b)) + ".txt";
    ofstream fout(k);
    fout.setf(ios::scientific);
    fout.precision(4);
    fout << n << "\n";
    fout << tid[0] << " " << y_h[0] << " " << y_h[1] << " " << y_h[2] << " " << y_h[0]+y_h[1]+y_h[2] <<"\n";

    //random
    random_device rd;
    mt19937 gen(rd());
    uniform_real_distribution<double> rr(0.0, 1.0);

    double a = 4.0; double c = 0.5; int N = 400;

    double E_S = 0; double E_I = 0; double E_R = 0;
    double E_S2 = 0; double E_I2 = 0; double E_R2 = 0;
    int X = 0;

    for(int x=0; x<=n; x++) {
        double P_SI = (a*y_h[0]*y_h[1])/N*dt;
        double P_IR = b*y_h[1]*dt;
        double P_RS = c*y_h[2]*dt;

        if(rr(gen) < P_SI) {
            if(rr(gen) < P_IR) {
                if(rr(gen) < P_RS) {
                    mc_rk4func(y_h, y_out, dt, a, b, c, N, 0);
                    //E_S = E_S + y_out[0]*P_RS*P_SI*P_IR;
                }
                else {
                    mc_rk4func(y_h, y_out, dt, a, b, c, N, 1);
                    //E_S = E_S + y_out[0]*P_SI*P_IR*(1-P_RS);
                }
            }
            else {
                if(rr(gen) < P_RS) {
                    mc_rk4func(y_h, y_out, dt, a, b, c, N, 2);
                    //E_S = E_S + y_out[0]*P_SI*(1-P_IR)*P_RS;
                }
                else {
                    mc_rk4func(y_h, y_out, dt, a, b, c, N, 3);
                    //E_S = E_S + y_out[0]*P_SI*(1-P_IR)*(1-P_RS);
                }
            }
        }
        else {
            if(rr(gen) < P_IR) {
                if(rr(gen) < P_RS) {
                    mc_rk4func(y_h, y_out, dt, a, b, c, N, 4);
                    //E_S = E_S + y_out[0]*(1-P_RS)*P_SI*P_IR;
                }
                else {
                    mc_rk4func(y_h, y_out, dt, a, b, c, N, 5);
                    //E_S = E_S + y_out[0]*(1-P_SI)*P_IR*(1-P_RS);
                }
            }
            else {
                if(rr(gen) < P_RS) {
                    mc_rk4func(y_h, y_out, dt, a, b, c, N, 6);
                    //E_S = E_S + y_out[0]*(1-P_SI)*(1-P_IR)*P_RS;
                }
                else {
                    mc_rk4func(y_h, y_out, dt, a, b, c, N, 7);
                    //E_S = E_S + y_out[0]*(1-P_SI)*(1-P_IR)*(1-P_RS);
                }
            }
        }
        fout << tid[x+1] << " " << y_out[0] << " " << y_out[1] << " " << y_out[2] << " " << y_out[0]+y_out[1]+y_out[2] <<"\n"; //skriver inn svarene i tekstfilen

        //finner forventningsverdiene ved/rundt equilibrium
        if(tid[x] > 4) {
            E_S = E_S + y_out[0]; E_S2 = E_S2 + y_out[0]*y_out[0];
            E_I = E_I + y_out[1]; E_I2 = E_I2 + y_out[1]*y_out[1];
            E_R = E_R + y_out[2]; E_R2 = E_R2 + y_out[2]*y_out[2];
            X = X + 1;
        }

        y_h[0] = y_out[0]; y_h[1] = y_out[1]; y_h[2] = y_out[2];
    }
    fout.close();
    //regner ut forventningsverdiene og standaravvikene
    E_S = E_S/X; E_I = E_I/X; E_R = E_R/X;
    E_S2 = E_S2/X; E_I2 = E_I2/X; E_R2 = E_R2/X;

    double std_S = sqrt(E_S2 - (E_S*E_S));
    double std_I = sqrt(E_I2 - (E_I*E_I));
    double std_R = sqrt(E_R2 - (E_R*E_R));
    //cout << E_S << " " << E_I << " " << E_R << "\n";
    //cout << std_S << " " << std_I << " " << std_R;

}

void mc_SIRfunc(double* in, double* out, double a, double b, double c, int N, int x) {
    double d = 0; double di = 0; double e = 0; double f = 0;
    if(x == 0) {
        out[0] = c*in[2] - (a*in[0]*in[1])/N - d*in[0] + e*N - f;
        out[1] = (a*in[0]*in[1])/N - b*in[1] - d*in[1] - di*in[1];
        out[2] = b*in[1] - c*in[2] - d*in[2] + f;
    }
    else if (x == 1) {
        out[0] = -(a*in[0]*in[1])/N - d*in[0] + e*N - f;
        out[1] = (a*in[0]*in[1])/N - b*in[1] - d*in[1] - di*in[1];
        out[2] = b*in[1] - d*in[2] + f;
    }
    else if (x == 2) {
        out[0] = c*in[2] - (a*in[0]*in[1])/N - d*in[0] + e*N - f;
        out[1] = (a*in[0]*in[1])/N - d*in[1] - di*in[1];
        out[2] = -c*in[2] - d*in[2] + f;
    }
    else if (x == 3) {
        out[0] = -(a*in[0]*in[1])/N - d*in[0] + e*N - f;
        out[1] = (a*in[0]*in[1])/N - d*in[1] - di*in[1];
        out[2] = -d*in[2] + f;
    }
    else if (x == 4) {
        out[0] = c*in[2] - d*in[0] + e*N - f;
        out[1] = -b*in[1] - d*in[1] - di*in[1];
        out[2] = b*in[1] - c*in[2] - d*in[2] + f;
    }
    else if (x == 5) {
        out[0] = -d*in[0] + e*N - f;
        out[1] = -b*in[1] - d*in[1] - di*in[1];
        out[2] = b*in[1] - d*in[2] + f;
    }
    else if (x == 6) {
        out[0] = c*in[2] - d*in[0] + e*N - f;
        out[1] = -d*in[1] - di*in[1];
        out[2] = -c*in[2] - d*in[2] + f;
    }
    else if (x == 7) {
        out[0] = -d*in[0] + e*N - f;
        out[1] = -d*in[1] - di*in[1];
        out[2] = -d*in[2] + f;
    }
}


void mc_rk4func(double* y_in, double* y_out, double dt, double a, double b, double c, int N, int x) {
    double k1[3], k2[3], k3[3], k4[3]; //lager fire tomme lister som skal holde k1, k2, k3 og k4 verdiene til S, I og R utregningene
    double y_k[3]; //lager en tom liste som skal holde paa midlertidig verdier for S, I og R

    //k1
    mc_SIRfunc(y_in, y_out, a, b, c, N, x);
    for(int i=0; i<=2; i++) {
        k1[i] = dt*y_out[i];
        y_k[i] = y_in[i] + k1[i]*0.5;
    }

    //k2
    mc_SIRfunc(y_k, y_out, a, b, c, N, x);
    for(int i=0; i<=2; i++) {
        k2[i] = dt*y_out[i];
        y_k[i] = y_in[i] + k2[i]*0.5;
    }

    //k3
    mc_SIRfunc(y_k, y_out, a, b, c, N, x);
    for(int i=0; i<=2; i++) {
        k3[i] = dt*y_out[i];
        y_k[i] = y_in[i] + k3[i];
    }

    //k4
    mc_SIRfunc(y_k, y_out, a, b, c, N, x);
    for(int i=0; i<=2; i++) {
        k4[i] = y_out[i]*dt;
    }
    //cout << x << " ";

    //calculate the new values for S, I and R
    for(int i=0; i<=2; i++) {
        y_out[i] = y_in[i] + (1.0/6.0)*(k1[i] + 2.0*k2[i] + 2.0*k3[i] + k4[i]);
    }

}
