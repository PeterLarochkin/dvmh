#include <stdio.h>          /* printf */
#include <string.h>         /* strcspn */
#include <stdlib.h>         /* strtol */
#include <math.h>



#define M (500)
#define N (500)
const double A1 = 0.0;
const double A2 = 4.0;
const double B1 = 0.0; 
const double B2 = 3.0;
const double h1 = ((/*A2-A1*/ 4.0)/ (double) M);
const double h2 = (/*B2-B1*/3.0 / (double) N);
const double epsilon = 0.000001;
double tau = 0.0;

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
#pragma dvm array align([i][j] with B[i][j])
double tau_r[M + 1][N + 1];

double q(double x, double y) {
    return x + y;
}

double u(double x, double y) {
    return sqrt(4+x*y);
}
double F(double x, double y) {
    double u_ = sqrt(4+x*y);
    return 1/(4*u_*u_*u_)*(x*x + y*y) + (x + y)*u_;
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
    } else /* if (x==A2 && y==B2 )*/ {
        double u = sqrt(4.0+x*y);
        return (h1*(x/(2*u)+u) + h2*(y/(2*u)+u)) / (h1 + h2);
    }
    /* 
    else {
        printf("ERROR:(%.10f, %.10f)", x, y);
    } */ 
}


void applyA_to_omega() {
    size_t i, j;
    #pragma dvm region
    {
        //with padding, inside "picture"
        #pragma dvm parallel([i][j] on A_omega[i][j]) shadow_renew(omega)
        for (i = 1; i < M; ++i) {
            for (j = 1; j < N; ++j) {
                // here is (7) equation works
                A_omega[i][j] = omega[i][j] * (2/(h1*h1) + 2/(h2*h2) + q(A1+ i*h1, B1+ j*h2)) + 
                                    omega[i-1][j] * (-1/(h1*h1)) +
                                    omega[i+1][j] * (-1/(h1*h1)) +
                                    omega[i][j-1] * (-1/(h2*h2)) +
                                    omega[i][j+1] * (-1/(h2*h2));
            }
        }

        #pragma dvm parallel([i] on A_omega[i][0]) shadow_renew(omega)
        for (i = 1; i < M; ++i) {
            // it's (10) equations
            // i=1,M-1
            // bottom applying
            A_omega[i][0] = -2/(h2*h2) * (omega[i][1] - omega[i][0]) +
                                ( q(A1+ i*h1,B1+ 0*h2) + 2/h2 ) * omega[i][0] -
                                1/(h1*h1)*(omega[i+1][0] - omega[i][0] - omega[i][0]+ omega[i-1][0]);
        }

        #pragma dvm parallel([i] on A_omega[i][N]) shadow_renew(omega)
        for (i = 1; i < M; ++i) {
            // it's (10) equations
            // i=1,M-1
            // top applying
            A_omega[i][N] = 2/(h2*h2) * (omega[i][N] - omega[i][N-1]) +
                                ( q(A1+ i*h1,B1+ N*h2) + 2/h2 ) * omega[i][N] -
                                1/(h1*h1)*(omega[i+1][N] - omega[i][N] - omega[i][N]+ omega[i-1][N]); 
        }

        #pragma dvm parallel([j] on A_omega[0][j]) shadow_renew(omega)
        for (j = 1; j < N; ++j) {
            // it's (9) equations
            // j=1,N-1
            // left applying
            A_omega[0][j] = -2/(h1*h1) * (omega[1][j] - omega[0][j]) + 
                                (q(A1+ 0*h1,B1+ j*h2) + 2/h1) * omega[0][j] - 
                                1/(h2*h2)*(omega[0][j+1] - omega[0][j] - omega[0][j]+ omega[0][j-1]);
        }

        #pragma dvm parallel([j] on A_omega[M][j]) shadow_renew(omega)
        for (j = 1; j < N; ++j) {
            // it's (9) equations
            // j=1,N-1
            // right applying
            A_omega[M][j] = 2/(h1*h1) * (omega[M][j] - omega[M-1][j]) + 
                                (q(A1+ M*h1,B1+ j*h2) + 2/h1) * omega[M][j] - 
                                1/(h2*h2)*(omega[M][j+1] - omega[M][j] - omega[M][j]+ omega[M][j-1]);
        }

        // remaining corner points
        // bottom left
        // it's (11) equation
        A_omega[0][0] = -2/(h1*h1)*(omega[1][0] - omega[0][0]) - 
                            2/(h2*h2)*(omega[0][1] - omega[0][0]) +
                            (q(A1+ 0*h1,B1+ 0*h2) + 2/h1 + 2/h2) * omega[0][0];

        // it's (12) equation
        // bottom right
        A_omega[M][0] = 2/(h1*h1)*(omega[M][0] - omega[M-1][0]) - 
                            2/(h2*h2)*(omega[M][1] - omega[M][0]) +
                            (q(A1+ M*h1,B1+ 0*h2) + 2/h1 + 2/h2) * omega[M][0];
        // it's (13) equation
        // top right
        A_omega[M][N] = 2/(h1*h1)*(omega[M][N] - omega[M-1][N]) +
                            2/(h2*h2)*(omega[M][N] - omega[M][N-1]) +
                            (q(A1+ M*h1,B1+ N*h2) + 2/h1 + 2/h2) * omega[M][N];
        // it's (14) equation
        // top left
        A_omega[0][N] = -2/(h1*h1)*(omega[1][N]- omega[0][N])+
                            2/(h2*h2)*(omega[0][N]- omega[0][N-1])+
                            (q(A1+ 0*h1,B1+ N*h2) + 2/h1 + 2/h2) * omega[0][N];
    }
}

// void applyA_to_r() {
//     size_t i, j;
//     #pragma dvm region
//     {
//         //with padding, inside "picture"
//         #pragma dvm parallel([i][j] on A_r[i][j]) shadow_renew(r)
//         for (i = 1; i < M; ++i) {
//             for (j = 1; j < N; ++j) {
//                 // here is (7) equation works
//                 A_r[i][j] = r[i][j] * (2/(h1*h1) + 2/(h2*h2) + q(A1+ i*h1, B1+ j*h2)) + 
//                                     r[i-1][j] * (-1/(h1*h1)) +
//                                     r[i+1][j] * (-1/(h1*h1)) +
//                                     r[i][j-1] * (-1/(h2*h2)) +
//                                     r[i][j+1] * (-1/(h2*h2));
//             }
//         }

//         #pragma dvm parallel([i] on A_r[i][0]) shadow_renew(r)
//         for (i = 1; i < M; ++i) {
//             // it's (10) equations
//             // i=1,M-1
//             // bottom applying
//             A_r[i][0] = -2/(h2*h2) * (r[i][1] - r[i][0]) +
//                                 ( q(A1+ i*h1,B1+ 0*h2) + 2/h2 ) * r[i][0] -
//                                 1/(h1*h1)*(r[i+1][0] - r[i][0] - r[i][0]+ r[i-1][0]);
//         }

//         #pragma dvm parallel([i] on A_r[i][N]) shadow_renew(r)
//         for (i = 1; i < M; ++i) {
//             // it's (10) equations
//             // i=1,M-1
//             // top applying
//             A_r[i][N] = 2/(h2*h2) * (r[i][N] - r[i][N-1]) +
//                                 ( q(A1+ i*h1,B1+ N*h2) + 2/h2 ) * r[i][N] -
//                                 1/(h1*h1)*(r[i+1][N] - r[i][N] - r[i][N]+ r[i-1][N]); 
//         }

//         #pragma dvm parallel([j] on A_r[0][j]) shadow_renew(r)
//         for (j = 1; j < N; ++j) {
//             // it's (9) equations
//             // j=1,N-1
//             // left applying
//             A_r[0][j] = -2/(h1*h1) * (r[1][j] - r[0][j]) + 
//                                 (q(A1+ 0*h1,B1+ j*h2) + 2/h1) * r[0][j] - 
//                                 1/(h2*h2)*(r[0][j+1] - r[0][j] - r[0][j]+ r[0][j-1]);
//         }

//         #pragma dvm parallel([j] on A_r[M][j]) shadow_renew(r)
//         for (j = 1; j < N; ++j) {
//             // it's (9) equations
//             // j=1,N-1
//             // right applying
//             A_r[M][j] = 2/(h1*h1) * (r[M][j] - r[M-1][j]) + 
//                                 (q(A1+ M*h1,B1+ j*h2) + 2/h1) * r[M][j] - 
//                                 1/(h2*h2)*(r[M][j+1] - r[M][j] - r[M][j]+ r[M][j-1]);
//         }

//         // remaining corner points
//         // bottom left
//         // it's (11) equation
//         A_r[0][0] = -2/(h1*h1)*(r[1][0] - r[0][0]) - 
//                             2/(h2*h2)*(r[0][1] - r[0][0]) +
//                             (q(A1+ 0*h1,B1+ 0*h2) + 2/h1 + 2/h2) * r[0][0];

//         // it's (12) equation
//         // bottom right
//         A_r[M][0] = 2/(h1*h1)*(r[M][0] - r[M-1][0]) - 
//                             2/(h2*h2)*(r[M][1] - r[M][0]) +
//                             (q(A1+ M*h1,B1+ 0*h2) + 2/h1 + 2/h2) * r[M][0];
//         // it's (13) equation
//         // top right
//         A_r[M][N] = 2/(h1*h1)*(r[M][N] - r[M-1][N]) +
//                             2/(h2*h2)*(r[M][N] - r[M][N-1]) +
//                             (q(A1+ M*h1,B1+ N*h2) + 2/h1 + 2/h2) * r[M][N];
//         // it's (14) equation
//         // top left
//         A_r[0][N] = -2/(h1*h1)*(r[1][N]- r[0][N])+
//                             2/(h2*h2)*(r[0][N]- r[0][N-1])+
//                             (q(A1+ 0*h1,B1+ N*h2) + 2/h1 + 2/h2) * r[0][N];
//     }
// }

// void getB() {
//     size_t i, j;
//     #pragma dvm region
//     {
//         //with padding, inside "picture"
//         #pragma dvm parallel([i][j] on B[i][j])
//         for (int i = 1; i < M; ++i) {
//             for (int j = 1; j < N; ++j) {
//                 // here is (7) equation works
//                 B[i][j] = F(A1+ i*h1,B1+ j*h2);
//             }
//         }
//         #pragma dvm parallel([i][j] on B[i][j])
//         for (int i=1; i<M; ++i) {
//             // it's (10) equations
//             // i=1,M-1
//             // top applying
//             B[i][N] = psi(A1+ i*h1, B1+ N*h2) * 2/h2 + F(A1 + i*h1, B1 + N*h2);
//             // bottom applying
//             B[i][0] = psi(A1+ i*h1, B1+ 0*h2) * 2/h2 + F(A1 + i*h1, B1 + 0*h2);
//         }
//         #pragma dvm parallel([i][j] on B[i][j])
//         for (int j=1; j<N; ++j) {
//             // it's (9) equations
//             // j=1,N-1
//             // right applying
//             B[M][j] = psi(A1+ M*h1, B1+ j*h2) * 2/h1 + F(A1 + M*h1, B1 + j*h2);
//             // left applying
//             B[0][j] = psi(A1+ 0*h1, B1+ j*h2) * 2/h1 + F(A1 + 0*h1, B1 + j*h2);
//         }

//         // remaining corner points
//         // bottom left
//         // it's (11) equation
//         B[0][0] = psi(A1+ 0*h1, B1+ 0*h2) * (2/h1 + 2/h2) + F(A1 + 0*h1, B1 + 0*h2);
//         // it's (12) equation
//         // bottom right
//         B[M][0] = psi(A1+ M*h1, B1+ 0*h2) * (2/h1 + 2/h2) + F(A1 + M*h1, B1 + 0*h2);
//         // it's (13) equation
//         // top right
//         B[M][N] = psi(A1+ M*h1, B1+ N*h2) * (2/h1 + 2/h2) + F(A1 + M*h1, B1 + N*h2);
//         // it's (14) equation
//         // top left
//         B[0][N] = psi(A1+ 0*h1, B1+ N*h2) * (2/h1 + 2/h2) + F(A1 + 0*h1, B1 + N*h2);
//     }
// }

//// void minus(double** first, double** second, double** whatWriteTo, double M, double N) {
////     for (size_t i = 0; i <= M; ++i) {
////             for (size_t j = 0; j <= N; ++j) {
////             whatWriteTo[i][j] = first[i][j] - second[i][j];
////         }
////     }
//// }

double ro(int index, int M_or_N) {
    if (index == 0 || index == M_or_N) {
        return 0.5;
    } else {
        return 1.0;
    }
}

// double scalarProduct_A_r_to_r() {
//     size_t i, j;
//     double sum_ = 0.0;
//     #pragma dvm actual(s)
//     #pragma dvm region
//     {
//         #pragma dvm parallel([i][j] on A_r[i][j]) reduction(sum(sum_))
//         for (i = 0; i < M + 1; ++i) {
//             for (j = 0; j < N + 1; ++j) {
//                 sum_ += h1 * h2 * ro(i, M)*ro(j, N) * A_r[i][j] * r[i][j];
//             }
//         }
//     }
//     return sum_;
// }

// double scalarProduct_A_r_to_A_r() {
//     size_t i, j;
//     double sum_ = 0.0;
//     #pragma dvm actual(sum_)
//     #pragma dvm region
//     {
//         #pragma dvm parallel([i][j] on A_r[i][j]) reduction(sum(sum_))
//         for (i = 0; i < M + 1; ++i) {
//             for (j = 0; j < N + 1; ++j) {
//                 sum_ += h1 * h2 * ro(i, M)*ro(j, N) * A_r[i][j] * A_r[i][j];
//             }
//         }
//     }
//     return sum_;
// }

// double scalarProduct_tau_r_to_tau_r() {
//     size_t i, j;
//     double sum_ = 0.0;
//     #pragma dvm actual(sum_)
//     #pragma dvm region
//     {
//         #pragma dvm parallel([i][j] on A_r[i][j]) reduction(sum(sum_))
//         for (i = 0; i < M + 1; ++i) {
//             for (j = 0; j < N + 1; ++j) {
//                 sum_ += h1 * h2 * ro(i, M)*ro(j, N) * r[i][j]*tau * r[i][j]*tau;
//             }
//         }
//     }
//     return sum_;
// }

//// void multiplyByNum(double** items, double num, double** whatWriteTo, double M, double N) {
////     for (size_t i = 0; i <= M; ++i) {
////             for (size_t j = 0; j <= N; ++j) {
////             whatWriteTo[i][j] = items[i][j]*num;
////         }
////     }
//// }
////

int main(/*int argc, char** argv*/) {
    // double sq_eps = epsilon * epsilon;
    // double squared_difference = sq_eps;
    // double start = dvmh_wtime(); 
    // double end = 0.0;
    // size_t i, j;

    // #pragma dvm region
    // {
    //     #pragma dvm parallel([i][j] on omega_next[i][j])
    //     for (i = 0; i <= M; ++i) {
    //         for (j = 0; j <= N; ++j) {
    //             // omega[i][j] = 0.0;
    //             omega_next[i][j] = 2.0;
    //             // B[i][j] = 0.0;
    //             // A_omega[i][j] = 0.0;
    //             // r[i][j] = 0.0;
    //             // A_r[i][j] = 0.0;
    //             // tau_r[i][j] = 0.0;
    //         }
    //     }
    // }

    // getB();

    // int count = 0;
    // while (squared_difference >= sq_eps)
    // {
    //     #pragma dvm region
    //     {
    //         #pragma dvm parallel([i][j] on omega[i][j])
    //         for (size_t i = 0; i <= M; ++i) {
    //             for (size_t j = 0; j <= N; ++j) {
    //                 omega[i][j] = omega_next[i][j];
    //             }
    //         }
    //     }
        
    //     // applyA(omega, A_omega, M, N, h1, h2, A1, B1);
    //     applyA_to_omega();
    //     #pragma dvm region
    //     {
    //     // minus(A_omega, B, r, M, N);
    //         #pragma dvm parallel([i][j] on r[i][j])
    //         for (i = 0; i <= M; ++i) {
    //             for (j = 0; j <= N; ++j) {
    //                 r[i][j] = A_omega[i][j] - B[i][j];
    //             }
    //         }
    //     }
    //     // applyA(r, A_r, M, N, h1, h2, A1, B1);
    //     applyA_to_r();
    //     // tau = scalarProduct(A_r, r, M, N, h1, h2) / scalarProduct(A_r, A_r, M, N, h1, h2);
    //     tau = scalarProduct_A_r_to_r()/scalarProduct_A_r_to_A_r();
    //     // multiplyByNum(r, tau, tau_r, M, N); tau_r = tau * r
    //     // minus(omega, tau_r, omega_next, M, N);
    //     #pragma dvm region
    //     {
    //         for (i = 0; i <= M; ++i) {
    //             for (j = 0; j <= N; ++j) {
    //                 omega_next[i][j] = omega[i][j] - r[i][j]*tau;
    //             }
    //         }
    //     }
    //     // squared_difference = scalarProduct(tau_r, tau_r, M, N, h1, h2);
    //     squared_difference = scalarProduct_tau_r_to_tau_r();
    //     count++;
    // }
    
    // dvmh_barrier();
    // end = dvmh_wtime();
    
    
    // double max_ = 0.0;
    // #pragma dvm actual(max_)
    // #pragma dvm region
    // {
    //     #pragma dvm parallel([i][j] on omega_next[i][j]) reduction(max(max_))
    //     for (i = 0; i < M + 1; ++i) {
    //         for (j = 0; j < N + 1; ++j) {
    //             double item = fabs(omega_next[i][j] - u(h1*i, h2*j));
    //             if (item > max_) {
    //                 max_ = item;
    //             }
                
    //         }
    //     }
    // }
    
    // printf("time:%.10f, diff:%.10f\n", (end-start), max_);
    
    return 0;
}