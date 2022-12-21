#include <stdio.h>          /* printf */
#include <string.h>         /* strcspn */
#include <stdlib.h>         /* strtol */
#include <math.h>
#include <time.h>

#define M (500)
#define N (1000)
#define epsilon (0.000001)
#define A1 (0.0)
#define A2 (4.0)
#define B1 (0.0)
#define B2 (3.0)
#define h1 (4.0/M)
#define h2 (3.0/N)

#pragma dvm array distribute [block][block]
double B[M + 1][N + 1];
#pragma dvm array align([i][j] with B[i][j])
double omega_next[M + 1][N + 1];
#pragma dvm array align([i][j] with B[i][j]), shadow[1:1][1:1]
double r[M + 1][N + 1];
#pragma dvm array align([i][j] with B[i][j])
double A_r[M + 1][N + 1];
#pragma dvm array align([i][j] with B[i][j])
double A_omega[M + 1][N + 1];
#pragma dvm array align([i][j] with B[i][j]), shadow[1:1][1:1]
double omega[M + 1][N + 1];



double q(double x, double y) {
    return x + y;
}

double F(double x, double y) {
    double u = sqrt(4+x*y);
    return 1/(4*u*u*u)*(x*x + y*y) + (x + y)*u;
}

double u(double x, double y) {
    return sqrt(4+x*y);
}




double psi(double x, double y) {
    
    if (x == A2 && B1<y && y<B2) {
        // 1
        double u_ = sqrt(4.0+x*y);
        return y/(2*u_)+u_;
    } else if (x == A1 && B1<y && y<B2) {
        // 2
        return -y/4+2;
    } else if (y==B2 && A1<x && x<A2) {
        // 3
        double u_ = sqrt(4.0+x*y);
        return x/(2*u_)+u_;
    } else if (y==B1 && A1<x && x<A2) {
        // 4
        return -x/4+2;

    } else if (x==A1 && y==B1 ) {
        
        return (h1*(-x/4+2) + h2*(-y/4+2)) / (h1 + h2);
    } else if (x==A1 && y==B2 ) {
        
        double u = sqrt(4.0+x*y);
        return (h1*(x/(2*u)+u) + h2*(-y/4+2)) / (h1 + h2);
    } else if (x==A2 && y==B1 ) {
        double u = sqrt(4.0+x*y);
        return (h1*(-x/4+2) + h2*(y/(2*u)+u)) / (h1 + h2);
    } else if (x==A2 && y==B2 ) {
        double u = sqrt(4.0+x*y);
        return (h1*(x/(2*u)+u) + h2*(y/(2*u)+u)) / (h1 + h2);
    } else {
        printf("ERROR:(%.10f, %.10f)", x, y);
    }  
}

double ro(int index, int M_or_N) {
    if (index == 0 || index == M_or_N) {
        return 0.5;
    } else {
        return 1.0;
    }
}

// void applyA(double** whatApplyTo, double** whatWriteTo, int M, int N, double h1, double h2, double A1, double B1) {

// }

// void getB(double** whatWriteTo, int M, int N, double h1, double h2, double A1, double A2, double B1, double B2) {

// }

// void minus(double** first, double** second, double** whatWriteTo, double M, double N) {

// }



// double scalarProduct(double** first, double** second, double M, double N, double h1, double h2) {

// }

// void multiplyByNum(double** items, double num, double** whatWriteTo, double M, double N) {

// }


int main(/*int argc, char** argv*/) {


    double sum_ = 0.0;
    double tau = 0.0;
    double start = dvmh_wtime(); 
    double end = 0.0;
    double sq_eps = epsilon * epsilon;
    double squared_difference = sq_eps;
    int i, j;

    #pragma dvm region 
    {
        #pragma dvm parallel([i][j] on omega_next[i][j])
        for (i = 0; i <= M; ++i) {
            for (j = 0; j <= N; ++j) {
                omega_next[i][j] = 2.0;
            }
        }
    }
    // getB(B, M, N, h1, h2, A1, A2, B1, B2);
    #pragma dvm region
    {
        #pragma dvm parallel([i][j] on B[i][j])
        for (i = 1; i < M; ++i) {
            for (int j = 1; j < N; ++j) {
                // here is (7) equation works
                B[i][j] = F(A1+ i*h1,B1+ j*h2);
            }
        }

        #pragma dvm parallel([i] on B[i][N])
        for (i = 1; i <= M; ++i) {
            if (i != M) {
                B[i][N] = psi(A1+ i*h1, B1+ N*h2) * 2/h2 + F(A1 + i*h1, B1 + N*h2); 
            } else {
                B[M][N] = psi(A1+ M*h1, B1+ N*h2) * (2/h1 + 2/h2) + F(A1 + M*h1, B1 + N*h2);
            }
            
        }

        #pragma dvm parallel([i] on B[i][0])
        for (i = 0; i < M; ++i) {   
            if (i == 0) {
                B[0][0] = psi(A1+ 0*h1, B1+ 0*h2) * (2/h1 + 2/h2) + F(A1 + 0*h1, B1 + 0*h2);
            } else {
                B[i][0] = psi(A1+ i*h1, B1+ 0*h2) * 2/h2 + F(A1 + i*h1, B1 + 0*h2);
            }
        }

        #pragma dvm parallel([j] on B[0][j])
        for (j = 1; j <= N; ++j) {
            if (j == N) {
                B[0][N] = psi(A1+ 0*h1, B1+ N*h2) * (2/h1 + 2/h2) + F(A1 + 0*h1, B1 + N*h2);
            } else {
                B[0][j] = psi(A1+ 0*h1, B1+ j*h2) * 2/h1 + F(A1 + 0*h1, B1 + j*h2);
            }
        }

        #pragma dvm parallel([j] on B[M][j])
        for (j = 0; j < N; ++j) {
            if (j == 0){
                B[M][0] = psi(A1+ M*h1, B1+ 0*h2) * (2/h1 + 2/h2) + F(A1 + M*h1, B1 + 0*h2);
            } else {
                B[M][j] = psi(A1+ M*h1, B1+ j*h2) * 2/h1 + F(A1 + M*h1, B1 + j*h2);
            }
        } 
    }

    int count = 0;
    while (squared_difference >= sq_eps && count < 5000)
    {
        #pragma dvm region
        {
            #pragma dvm parallel([i][j] on omega[i][j])
            for (i = 0; i <= M; ++i) {
                for (j = 0; j <= N; ++j) {
                    omega[i][j] = omega_next[i][j];
                }
            }
        }

        // applyA(omega, A_omega, M, N, h1, h2, A1, B1);
        //with padding, inside "picture"
        #pragma dvm region
        {
            #pragma dvm parallel([i][j] on A_omega[i][j]) shadow_renew(omega) 
            for (i = 1; i < M; ++i) {
                for (j = 1; j < N; ++j) {
                    A_omega[i][j] = omega[i][j] * (2/(h1*h1) + 2/(h2*h2) + q(A1+ i*h1, B1+ j*h2)) + 
                                        omega[i-1][j] * (-1/(h1*h1)) +
                                        omega[i+1][j] * (-1/(h1*h1)) +
                                        omega[i][j-1] * (-1/(h2*h2)) +
                                        omega[i][j+1] * (-1/(h2*h2));
                }
            }
            #pragma dvm parallel ([i] on A_omega[i][N]) shadow_renew(omega) 
            for (i = 0; i < M; ++i) {
                if (i != 0) {
                    A_omega[i][N] = 2/(h2*h2) * (omega[i][N] - omega[i][N-1]) +
                                    ( q(A1+ i*h1,B1+ N*h2) + 2/h2 ) * omega[i][N] -
                                    1/(h1*h1)*(omega[i+1][N] - omega[i][N] - omega[i][N]+ omega[i-1][N]); 
                } else {
                    A_omega[0][N] = -2/(h1*h1)*(omega[1][N]- omega[0][N])+
                                2/(h2*h2)*(omega[0][N]- omega[0][N-1])+
                                (q(A1+ 0*h1,B1+ N*h2) + 2/h1 + 2/h2) * omega[0][N];
                }
            }

            #pragma dvm parallel([i] on A_omega[i][0]) shadow_renew(omega) 
            for (i = 1; i <= M; ++i) {
                if (i < M)
                {
                    A_omega[i][0] = -2/(h2*h2) * (omega[i][1] - omega[i][0]) +
                                    ( q(A1+ i*h1,B1+ 0*h2) + 2/h2 ) * omega[i][0] -
                                    1/(h1*h1)*(omega[i+1][0] - omega[i][0] - omega[i][0]+ omega[i-1][0]);
                } else {
                    A_omega[i][0] = 2/(h1*h1)*(omega[M][0] - omega[i-1][0]) - 
                                2/(h2*h2)*(omega[i][1] - omega[i][0]) +
                                (q(A1+ i*h1,B1+ 0*h2) + 2/h1 + 2/h2) * omega[i][0];
                }
            }

            #pragma dvm parallel([j] on A_omega[0][j]) shadow_renew(omega) 
            for (j = 0; j < N; ++j) {
                if (j == 0) {
                    A_omega[0][j] = -2/(h1*h1)*(omega[1][j] - omega[0][j]) - 
                                2/(h2*h2)*(omega[0][j+1] - omega[0][j]) +
                                (q(A1+ 0*h1,B1+ j*h2) + 2/h1 + 2/h2) * omega[0][j];
                } else {
                    A_omega[0][j] = -2/(h1*h1) * (omega[1][j] - omega[0][j]) + 
                                    (q(A1+ 0*h1,B1+ j*h2) + 2/h1) * omega[0][j] - 
                                    1/(h2*h2)*(omega[0][j+1] - omega[0][j] - omega[0][j]+ omega[0][j-1]);
                }
            }
            #pragma dvm parallel([j] on A_omega[M][j]) shadow_renew(omega) 
            for (j = 1; j <= N; ++j) {
                if (j < N) {
                    A_omega[M][j] = 2/(h1*h1) * (omega[M][j] - omega[M-1][j]) + 
                                    (q(A1+ M*h1,B1+ j*h2) + 2/h1) * omega[M][j] - 
                                    1/(h2*h2)*(omega[M][j+1] - omega[M][j] - omega[M][j]+ omega[M][j-1]);
                } else {
                    A_omega[M][j] = 2/(h1*h1)*(omega[M][j] - omega[M-1][j]) +
                                2/(h2*h2)*(omega[M][j] - omega[M][j-1]) +
                                (q(A1+ M*h1,B1+ j*h2) + 2/h1 + 2/h2) * omega[M][j];
                }
            }
        }


        // minus(A_omega, B, r, M, N);
        #pragma dvm region
        {
            #pragma dvm parallel([i][j] on r[i][j])
            for (i = 0; i <= M; ++i) {
                    for (j = 0; j <= N; ++j) {
                        r[i][j] = A_omega[i][j] - B[i][j];
                }
            } 
        }  




        // applyA(r, A_r, M, N, h1, h2, A1, B1);
        #pragma dvm region
        {

            //with padding, inside "picture"
            #pragma dvm parallel([i][j] on A_r[i][j]) shadow_renew(r) 
            for (i = 1; i < M; ++i) {
                for (j = 1; j < N; ++j) {
                    // here is (7) equation works
                    A_r[i][j] = r[i][j] * (2/(h1*h1) + 2/(h2*h2) + q(A1+ i*h1, B1+ j*h2)) + 
                                        r[i-1][j] * (-1/(h1*h1)) +
                                        r[i+1][j] * (-1/(h1*h1)) +
                                        r[i][j-1] * (-1/(h2*h2)) +
                                        r[i][j+1] * (-1/(h2*h2));
                }
            }
            #pragma dvm parallel ([i] on A_r[i][N]) shadow_renew(r) 
            for (i = 0; i < M; ++i) {
                if (i != 0) {
                    A_r[i][N] = 2/(h2*h2) * (r[i][N] - r[i][N-1]) +
                                    ( q(A1+ i*h1,B1+ N*h2) + 2/h2 ) * r[i][N] -
                                    1/(h1*h1)*(r[i+1][N] - r[i][N] - r[i][N]+ r[i-1][N]); 
                } else {
                    A_r[0][N] = -2/(h1*h1)*(r[1][N]- r[0][N])+
                                2/(h2*h2)*(r[0][N]- r[0][N-1])+
                                (q(A1+ 0*h1,B1+ N*h2) + 2/h1 + 2/h2) * r[0][N];
                }
            }
            #pragma dvm parallel([i] on A_r[i][0]) shadow_renew(r) 
            for (i = 1; i <= M; ++i) {
                if (i < M)
                {
                    A_r[i][0] = -2/(h2*h2) * (r[i][1] - r[i][0]) +
                                    ( q(A1+ i*h1,B1+ 0*h2) + 2/h2 ) * r[i][0] -
                                    1/(h1*h1)*(r[i+1][0] - r[i][0] - r[i][0]+ r[i-1][0]);
                } else {
                    A_r[i][0] = 2/(h1*h1)*(r[M][0] - r[i-1][0]) - 
                                2/(h2*h2)*(r[i][1] - r[i][0]) +
                                (q(A1+ i*h1,B1+ 0*h2) + 2/h1 + 2/h2) * r[i][0];
                }
            }
            #pragma dvm parallel([j] on A_r[0][j]) shadow_renew(r) 
            for (j = 0; j < N; ++j) {
                if (j == 0) {
                    A_r[0][j] = -2/(h1*h1)*(r[1][j] - r[0][j]) - 
                                2/(h2*h2)*(r[0][j+1] - r[0][j]) +
                                (q(A1+ 0*h1,B1+ j*h2) + 2/h1 + 2/h2) * r[0][j];
                } else {
                    A_r[0][j] = -2/(h1*h1) * (r[1][j] - r[0][j]) + 
                                    (q(A1+ 0*h1,B1+ j*h2) + 2/h1) * r[0][j] - 
                                    1/(h2*h2)*(r[0][j+1] - r[0][j] - r[0][j]+ r[0][j-1]);
                }
            }

            #pragma dvm parallel([j] on A_r[M][j]) shadow_renew(r) 
            for (j = 1; j <= N; ++j) {
                if (j < N) {
                    A_r[M][j] = 2/(h1*h1) * (r[M][j] - r[M-1][j]) + 
                                    (q(A1+ M*h1,B1+ j*h2) + 2/h1) * r[M][j] - 
                                    1/(h2*h2)*(r[M][j+1] - r[M][j] - r[M][j]+ r[M][j-1]);
                } else {
                    A_r[M][j] = 2/(h1*h1)*(r[M][j] - r[M-1][j]) +
                                2/(h2*h2)*(r[M][j] - r[M][j-1]) +
                                (q(A1+ M*h1,B1+ j*h2) + 2/h1 + 2/h2) * r[M][j];
                }
            }
        }
        // tau = scalarProduct(A_r, r, M, N, h1, h2) / scalarProduct(A_r, A_r, M, N, h1, h2);
        sum_ = 0.0;
        #pragma dvm actual(sum_)
        #pragma dvm region
        {
            #pragma dvm parallel([i][j] on A_r[i][j]) reduction(sum(sum_))
            for (i = 0; i <= M; ++i) {
                for (j = 0; j <= N; ++j) {
                    sum_ = sum_ + h1*h2*ro(i, M)*ro(j, N)*A_r[i][j] * r[i][j];
                }
            }
        }
        tau = sum_;
        sum_ = 0.0;
        #pragma dvm actual(sum_)
        #pragma dvm region
        {
            #pragma dvm parallel([i][j] on A_r[i][j]) reduction(sum(sum_))
            for (i = 0; i <= M; ++i) {
                for (j = 0; j <= N; ++j) {
                    sum_ = sum_ + h1*h2*ro(i, M)*ro(j, N)*A_r[i][j] * A_r[i][j];
                }
            }
        }
        tau = tau / sum_;


        // multiplyByNum(r, tau, tau_r, M, N);
        // minus(omega, tau_r, omega_next, M, N);
        #pragma dvm region
        {
            #pragma dvm parallel([i][j] on omega_next[i][j])
            for (i = 0; i <= M; ++i) {
                for (j = 0; j <= N; ++j) {
                    omega_next[i][j] = omega[i][j] - tau*r[i][j];
                }
            }
        }

        // squared_difference = scalarProduct(tau_r, tau_r, M, N, h1, h2);
        sum_ = 0.0;
        #pragma dvm actual(sum_)
        #pragma dvm region
        {
            #pragma dvm parallel([i][j] on r[i][j]) reduction(sum(sum_))
            for (i = 0; i <= M; ++i) {
                for (j = 0; j <= N; ++j) {
                    sum_ = sum_ + h1*h2*ro(i, M)*ro(j, N)*tau*r[i][j] * tau*r[i][j];
                }
            }
        }

        double squared_difference = sum_;
    
        if (count % 500 == 0)
            printf("n:%d, diff:%.10f\n", count, sqrt(squared_difference));
        count++;
    }
    
    end = dvmh_wtime();
    double max_ = 0.0;
    #pragma dvm actual(max_)
    #pragma dvm region
    {
        #pragma dvm parallel([i][j] on omega_next[i][j]) reduction(max(max_))
        for (i = 0; i < M + 1; ++i) {
            for (j = 0; j < N + 1; ++j) {
                double item = fabs(omega_next[i][j] - u(h1*i, h2*j));
                if (item > max_) {
                    max_ = item;
                }
                
            }
        }
    }
    
    printf("time:%.10f, max_diff:%.10f\n", (end-start), max_);
    return 0;
}
