#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include <mpi.h>
#include <string.h>


#define MIN(a,b) (((a)<(b))?(a):(b))
#define MAX(a,b) (((a)>(b))?(a):(b))


struct Info {
  int rank;
  int size;
  // Process topological coords
  int coords[2];
  int left_rank, right_rank, up_rank, down_rank;
  // Process local domain size
  size_t m, n;
  size_t a1, a2, b1, b2;

};
typedef struct Info Info_t;

int log2_(int num) {
  if (num <= 0)
    return -1;
  int power = 0;
  while ((num & 1) == 0) {
    ++power;
    num = num >> 1; 
  } 
  if ((num >> 1) != 0)
    return -1;
  return power;
}

int split(size_t M, size_t N, int power) {
  double m = (double)M;
  double n = (double)N;
  int py = 0;
  for (int i = 0; i < power; ++i) {
    if (m > n) {
      n /= 2.0;
      ++py;
    } else {
      m /= 2.0;
    }
  } 
  return py;
}

double u(double x, double y) {
    return sqrt(4+x*y);
}

void partitioningDomain(size_t M, size_t N, MPI_Comm *Comm, int rank, int size, Info_t* info) {

    // ProcNum = 2^(power), power = px + py
    int power, px, py;
    // dims[0] = 2^px, dims[1] = 2^py
    int dims[2];
    // 2D topology dimension
    const int ndims = 2;
    // Topologically grid is not closed
    int periods[2] = {0, 0};
    // Init MPI lib
    MPI_Comm_size(MPI_COMM_WORLD, &size);
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);

    if (size == 20) {
        dims[0] = 5;
        dims[1] = 4;
    } else if (size == 40) {
        dims[0] = 8;
        dims[1] = 5;
    } else {
        MPI_Finalize();
        exit(1);
    }
    
    
    int m = M / dims[0]; 
    int n = N / dims[1];
    int rx = M + 1 - dims[0] * m;
    int ry = N + 1 - dims[1] * n;
    
    int coords[2];

    // printf("dims:[%d, %d]", dims[0], dims[1]);
    MPI_Cart_create(MPI_COMM_WORLD, ndims, dims, periods, 1, Comm);
    
    MPI_Cart_coords(*Comm, rank, ndims, coords);
    
    int a1 = MIN(rx, coords[0]) * (m + 1) + MAX(0, (coords[0] - rx)) * m;
    int b1 = MIN(ry, coords[1]) * (n + 1) + MAX(0, (coords[1] - ry)) * n;
    int a2 = a1 + m + (coords[0] < rx ? 1 : 0);
    int b2 = b1 + n + (coords[1] < ry ? 1 : 0);
    m = a2 - a1;
    n = b2 - b1;
    int up, left, down, right;
    // Get process neighbors
    MPI_Cart_shift(*Comm, 1, -1, &up, &down);
    MPI_Cart_shift(*Comm, 0, 1, &left, &right);
    // printf("----------------- (%d,%d)\n|            (%d)\n|             |\n|       (%d)-(%d, %d,%d)-(%d)  \n|             |\n|           (%d) \n(%d,%d)---------------\n", a2,b2, up, left, rank, coords[0], coords[1], right, down, a1,b1);
    info->rank = rank;
    info->coords[0] = coords[0];
    info->coords[1] = coords[1];
    info->left_rank=left;
    info->right_rank=right;
    info->up_rank = up;
    info->down_rank = down;
    info->m = m; 
    info->n = n;
    info->a1 = a1;
    info->a2 = a2;
    info->b1 = b1;
    info->b2 = b2;
    info->size = dims[0] * dims[1];
}

double q(double x, double y) {
    return x + y;
}


double F(double x, double y) {
    double u_ = u(x, y);
    return 1/(4*u_*u_*u_)*(x*x + y*y) + q(x, y)*u_;
}

double psi(double x, double y, double A1, double A2, double B1, double B2, double h1, double h2) {
    
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
        // printf("ERROR:(%.10f, %.10f)", x, y);
        double u_ = sqrt(4.0+x*y);
        return y/(2*u_)+u_;
    }  
}

void applyA(double** whatApplyTo, double** whatWriteTo, int M, int N, double h1, double h2, double A1, double A2, double B1, double B2, Info_t* info) {

    int m = info->m;
    int n = info->n;
    int a1 = info->a1;
    int b1 = info->b1;
    int left = info->left_rank;
    int right = info->right_rank;
    int up = info->up_rank;
    int down = info->down_rank;
    //with padding, inside, inside "picture"
    for (int i = 2; i <= m-1; ++i) {
        for (int j = 2; j <= n-1; ++j) {
            // here is (7) equation works
            whatWriteTo[i][j] = whatApplyTo[i][j] * (2/(h1*h1) + 2/(h2*h2) + q((a1 + i - 1)*h1, (b1 + j - 1)*h2)) + 
                                whatApplyTo[i-1][j] * (-1/(h1*h1)) +
                                whatApplyTo[i+1][j] * (-1/(h1*h1)) +
                                whatApplyTo[i][j-1] * (-1/(h2*h2)) +
                                whatApplyTo[i][j+1] * (-1/(h2*h2));
        }
    }
    // let's look on the top border
    if (up <0) {
        // if I'm near top, I need to use top boarding equation (10)
        // check (b1 + n) == N ?
        // if ((b1 + n - 1) == N) {
        //     printf("Yes\n");
        // }
        for (int i = 2; i <= m-1; ++i) { 
            // whatWriteTo[i][n] = psi((a1 + i - 1)*h1, (b1 + n - 1)*h2, A1, A2, B1, B2) * 2/h2 + F((a1 + i - 1)*h1, (b1 + n - 1)*h2);
            whatWriteTo[i][n] = 2/(h2*h2) * (whatApplyTo[i][n] - whatApplyTo[i][n-1]) +
                            ( q((a1 + i - 1)*h1, (b1 + n - 1)*h2) + 2/h2 ) * whatApplyTo[i][n] -
                            1/(h1*h1)*(whatApplyTo[i+1][n] - whatApplyTo[i][n] - whatApplyTo[i][n]+ whatApplyTo[i-1][n]);
        }
    } else {
        // I'm not near the top border
        for (int i = 2; i <= m-1; ++i) { 
            whatWriteTo[i][n] = whatApplyTo[i][n] * (2/(h1*h1) + 2/(h2*h2) + q((a1 + i - 1)*h1, (b1 + n - 1)*h2)) + 
                                whatApplyTo[i-1][n] * (-1/(h1*h1)) +
                                whatApplyTo[i+1][n] * (-1/(h1*h1)) +
                                whatApplyTo[i][n-1] * (-1/(h2*h2)) +
                                whatApplyTo[i][n+1] * (-1/(h2*h2)); // <- here is using border from another cell
        }
    }



    // let's look on the bottom border
    if (down <0) {
    // if I'm near the bottom border
        for (int i = 2; i <= m-1; ++i) { 
            whatWriteTo[i][1] = -2/(h2*h2) * (whatApplyTo[i][2] - whatApplyTo[i][1]) +
                            ( q((a1 + i - 1)*h1, (b1 + 1 - 1)*h2) + 2/h2 ) * whatApplyTo[i][1] -
                            1/(h1*h1)*(whatApplyTo[i+1][1] - whatApplyTo[i][1] - whatApplyTo[i][1]+ whatApplyTo[i-1][1]);
        }
    } else {
    // if I'm not near the bottom border
        for (int i = 2; i <= m-1; ++i) { 
            whatWriteTo[i][1] = whatApplyTo[i][1] * (2/(h1*h1) + 2/(h2*h2) + q((a1 + i - 1)*h1, (b1 + 1 - 1)*h2)) + 
                                whatApplyTo[i-1][1] * (-1/(h1*h1)) +
                                whatApplyTo[i+1][1] * (-1/(h1*h1)) +
                                whatApplyTo[i][1-1] * (-1/(h2*h2)) + // <- here is using border from another cell
                                whatApplyTo[i][1+1] * (-1/(h2*h2));
        }
    }

    // let's look on the right board
    if (right <0) {
        // I'm near right border
        // if ((a1 + m - 1) == M) {
        //     printf("Yes\n");
        // }
        for (int j = 2; j <= n-1; ++j) { 
            whatWriteTo[m][j] = 2/(h1*h1) * (whatApplyTo[m][j] - whatApplyTo[m-1][j]) + 
                            (q((a1 + m - 1)*h1, (b1 + j - 1)*h2) + 2/h1) * whatApplyTo[m][j] - 
                            1/(h2*h2)*(whatApplyTo[m][j+1] - whatApplyTo[m][j] - whatApplyTo[m][j]+ whatApplyTo[m][j-1]);
        }
    } else {
        for (int j = 2; j <= n-1; ++j) { 
            whatWriteTo[m][j] = whatApplyTo[m][j] * (2/(h1*h1) + 2/(h2*h2) + q((a1 + m - 1)*h1, (b1 + j - 1)*h2)) + 
                                whatApplyTo[m-1][j] * (-1/(h1*h1)) +
                                whatApplyTo[m+1][j] * (-1/(h1*h1)) + // <- here is using border from another cell
                                whatApplyTo[m][j-1] * (-1/(h2*h2)) +
                                whatApplyTo[m][j+1] * (-1/(h2*h2));
        }
    }

    //let's look on the left border
    if (left <0) {
        // if ((a1 + 1 - 1) == 0) {
        //     printf("Yes\n");
        // }
        for (int j = 2; j <= n-1; ++j) { 
            whatWriteTo[1][j] = -2/(h1*h1) * (whatApplyTo[2][j] - whatApplyTo[1][j]) + 
                            (q((a1 + 1 - 1)*h1, (b1 + j - 1)*h2) + 2/h1) * whatApplyTo[1][j] - 
                            1/(h2*h2)*(whatApplyTo[1][j+1] - whatApplyTo[1][j] - whatApplyTo[1][j]+ whatApplyTo[1][j-1]);
        }
    } else {
        for (int j = 2; j <= n-1; ++j) { 
            whatWriteTo[1][j] = whatApplyTo[1][j] * (2/(h1*h1) + 2/(h2*h2) + q((a1 + 1 - 1)*h1, (b1 + j - 1)*h2)) + 
                                whatApplyTo[1-1][j] * (-1/(h1*h1)) + // <- here is using border from another cell
                                whatApplyTo[1+1][j] * (-1/(h1*h1)) +
                                whatApplyTo[1][j-1] * (-1/(h2*h2)) +
                                whatApplyTo[1][j+1] * (-1/(h2*h2));
        }
    }

    // look on the down-left point
    if (down<0 && left<0) {
        // printf("Yes\n");
        // it's (11) equation
        whatWriteTo[1][1] = -2/(h1*h1)*(whatApplyTo[2][1] - whatApplyTo[1][1]) - 
                         2/(h2*h2)*(whatApplyTo[1][2] - whatApplyTo[1][1]) +
                         (q((a1 + 1 - 1)*h1, (b1 + 1 - 1)*h2) + 2/h1 + 2/h2) * whatApplyTo[1][1];
    } else if (down<0 && left >=0 ) {
        whatWriteTo[1][1] = -2/(h2*h2) * (whatApplyTo[1][2] - whatApplyTo[1][1]) +
                            ( q((a1 + 1 - 1)*h1, (b1 + 1 - 1)*h2) + 2/h2 ) * whatApplyTo[1][1] -
                            1/(h1*h1)*(whatApplyTo[1+1][1] - whatApplyTo[1][1] - whatApplyTo[1][1]+ whatApplyTo[1-1][1]); // <- here is using border from another cell
    } else if (down>=0 && left <0) {
        whatWriteTo[1][1] = -2/(h1*h1) * (whatApplyTo[2][1] - whatApplyTo[1][1]) + 
                            (q((a1 + 1 - 1)*h1, (b1 + 1 - 1)*h2) + 2/h1) * whatApplyTo[1][1] - 
                            1/(h2*h2)*(whatApplyTo[1][1+1] - whatApplyTo[1][1] - whatApplyTo[1][1]+ whatApplyTo[1][1-1]);// <- here is using border from another cell
    } else if (down>=0 && left >=0) {
        whatWriteTo[1][1] = whatApplyTo[1][1] * (2/(h1*h1) + 2/(h2*h2) + q((a1 + 1 - 1)*h1, (b1 + 1 - 1)*h2)) + 
                                whatApplyTo[1-1][1] * (-1/(h1*h1)) + // <- here is using border from another cell
                                whatApplyTo[1+1][1] * (-1/(h1*h1)) +
                                whatApplyTo[1][1-1] * (-1/(h2*h2)) + // <- here is using border from another cell
                                whatApplyTo[1][1+1] * (-1/(h2*h2));
    }



    // look on the bottom-right points
    if (down<0 && right<0) {
        // printf("Yes\n");
        // it's (12) equation
        
        whatWriteTo[m][1] = 2/(h1*h1)*(whatApplyTo[m][1] - whatApplyTo[m-1][1]) - 
                            2/(h2*h2)*(whatApplyTo[m][2] - whatApplyTo[m][1]) +
                            (q((a1 + m - 1)*h1, (b1 + 1 - 1)*h2) + 2/h1 + 2/h2) * whatApplyTo[m][1];
    } else if (down<0 && right >=0 ) {
        whatWriteTo[m][1] = -2/(h2*h2) * (whatApplyTo[m][2] - whatApplyTo[m][1]) +
                            ( q((a1 + m - 1)*h1, (b1 + 1 - 1)*h2) + 2/h2 ) * whatApplyTo[m][1] -
                            1/(h1*h1)*(whatApplyTo[m+1][1] - whatApplyTo[m][1] - whatApplyTo[m][1]+ whatApplyTo[m-1][1]); // <- here is using border from another cell
    } else if (down>=0 && right <0) {
        whatWriteTo[m][1] = 2/(h1*h1) * (whatApplyTo[m][1] - whatApplyTo[m-1][1]) + 
                            (q((a1 + m - 1)*h1, (b1 + 1 - 1)*h2) + 2/h1) * whatApplyTo[m][1] - 
                            1/(h2*h2)*(whatApplyTo[m][1+1] - whatApplyTo[m][1] - whatApplyTo[m][1]+ whatApplyTo[m][1-1]);
    } else if (down>=0 && right >=0) {
        whatWriteTo[m][1] = whatApplyTo[m][1] * (2/(h1*h1) + 2/(h2*h2) + q((a1 + m - 1)*h1, (b1 + 1 - 1)*h2)) + 
                                whatApplyTo[m-1][1] * (-1/(h1*h1)) +
                                whatApplyTo[m+1][1] * (-1/(h1*h1)) +// <- here is using border from another cell
                                whatApplyTo[m][1-1] * (-1/(h2*h2)) +// <- here is using border from another cell
                                whatApplyTo[m][1+1] * (-1/(h2*h2));
    }

    // look on the top-right points
    if (up<0 && right<0) {
        // printf("Yes\n");
        // it's (13) equation
        whatWriteTo[m][n] = 2/(h1*h1)*(whatApplyTo[m][n] - whatApplyTo[m-1][n]) +
                        2/(h2*h2)*(whatApplyTo[m][n] - whatApplyTo[m][n-1]) +
                        (q((a1 + m - 1)*h1, (b1 + n - 1)*h2) + 2/h1 + 2/h2) * whatApplyTo[m][n];
    } else if (up<0 && right >=0 ) {
        whatWriteTo[m][n] = 2/(h2*h2) * (whatApplyTo[m][n] - whatApplyTo[m][n-1]) +
                            ( q((a1 + m - 1)*h1, (b1 + n - 1)*h2) + 2/h2 ) * whatApplyTo[m][n] -
                            1/(h1*h1)*(whatApplyTo[m+1][n] - whatApplyTo[m][n] - whatApplyTo[m][n]+ whatApplyTo[m-1][n]);
    } else if (up>=0 && right <0) {
        whatWriteTo[m][n] = 2/(h1*h1) * (whatApplyTo[m][n] - whatApplyTo[m-1][n]) + 
                            (q((a1 + m - 1)*h1, (b1 + n - 1)*h2) + 2/h1) * whatApplyTo[m][n] - 
                            1/(h2*h2)*(whatApplyTo[m][n+1] - whatApplyTo[m][n] - whatApplyTo[m][n]+ whatApplyTo[m][n-1]);
    } else if (up>=0 && right >=0) {
        whatWriteTo[m][n] = whatApplyTo[m][n] * (2/(h1*h1) + 2/(h2*h2) + q((a1 + m - 1)*h1, (b1 + n - 1)*h2)) + 
                                whatApplyTo[m-1][n] * (-1/(h1*h1)) +
                                whatApplyTo[m+1][n] * (-1/(h1*h1)) +
                                whatApplyTo[m][n-1] * (-1/(h2*h2)) +
                                whatApplyTo[m][n+1] * (-1/(h2*h2));
    }

    // look on the top-left points
    if (up<0 && left<0) {
        // printf("Yes\n");
        // it's (14) equation
        
        whatWriteTo[1][n] = -2/(h1*h1)*(whatApplyTo[2][n]- whatApplyTo[1][n])+
                            2/(h2*h2)*(whatApplyTo[1][n]- whatApplyTo[1][n-1])+
                            (q((a1 + 1 - 1)*h1, (b1 + n - 1)*h2) + 2/h1 + 2/h2) * whatApplyTo[1][n];
    } else if (up<0 && left >=0 ) {
        whatWriteTo[1][n] = 2/(h2*h2) * (whatApplyTo[1][n] - whatApplyTo[1][n-1]) +
                            ( q((a1 + 1 - 1)*h1, (b1 + n - 1)*h2) + 2/h2 ) * whatApplyTo[1][n] -
                            1/(h1*h1)*(whatApplyTo[1+1][n] - whatApplyTo[1][n] - whatApplyTo[1][n]+ whatApplyTo[1-1][n]);
    } else if (up>=0 && left <0) {
        whatWriteTo[1][n] = -2/(h1*h1) * (whatApplyTo[2][n] - whatApplyTo[1][n]) + 
                            (q((a1 + 1 - 1)*h1, (b1 + n - 1)*h2) + 2/h1) * whatApplyTo[1][n] - 
                            1/(h2*h2)*(whatApplyTo[1][n+1] - whatApplyTo[1][n] - whatApplyTo[1][n]+ whatApplyTo[1][n-1]);
    } else if (up>=0 && left >=0) {
        whatWriteTo[1][n] = whatApplyTo[1][n] * (2/(h1*h1) + 2/(h2*h2) + q((a1 + 1 - 1)*h1, (b1 + n - 1)*h2)) + 
                                whatApplyTo[1-1][n] * (-1/(h1*h1)) +
                                whatApplyTo[1+1][n] * (-1/(h1*h1)) +
                                whatApplyTo[1][n-1] * (-1/(h2*h2)) +
                                whatApplyTo[1][n+1] * (-1/(h2*h2));
    }
}

void getB(double** whatWriteTo, int M, int N, double h1, double h2, double A1, double A2, double B1, double B2, Info_t* info) {
    int m = info->m;
    int n = info->n;
    int a1 = info->a1;
    int b1 = info->b1;
    int left = info->left_rank;
    int right = info->right_rank;
    int up = info->up_rank;
    int down = info->down_rank;
    //with padding, inside, inside "picture"
    for (int i = 2; i <= m-1; ++i) {
        for (int j = 2; j <= n-1; ++j) {
            whatWriteTo[i][j] = F((a1 + i - 1)*h1, (b1 + j - 1)*h2);
        }
    }
    // let's look on the top border
    if (up <0) {
        // if I'm near top, I need to use top boarding equation (10)
        for (int i = 2; i <= m-1; ++i) { 
            whatWriteTo[i][n] = psi((a1 + i - 1)*h1, (b1 + n - 1)*h2, A1, A2, B1, B2, h1, h2) * 2/h2 + F((a1 + i - 1)*h1, (b1 + n - 1)*h2);
        }
    } else {
        // I'm not near the top border
        for (int i = 2; i <= m-1; ++i) { 
            whatWriteTo[i][n] = F((a1 + i - 1)*h1, (b1 + n - 1)*h2);
        }
    }

    // let's look on the bottom border
    if (down <0) {
    // if I'm near the bottom border
        for (int i = 2; i <= m-1; ++i) { 
            whatWriteTo[i][1] = psi((a1 + i - 1)*h1, (b1 + 1 - 1)*h2, A1, A2, B1, B2, h1, h2) * 2/h2 + F((a1 + i - 1)*h1, (b1 + 1 - 1)*h2);
        }
    } else {
    // if I'm not near the bottom border
        for (int i = 2; i <= m-1; ++i) { 
            whatWriteTo[i][1] = F((a1 + i - 1)*h1, (b1 + 1 - 1)*h2);
        }
    }

    // let's look on the right board
    if (right <0) {
        for (int j = 2; j <= n-1; ++j) { 
            whatWriteTo[m][j] = psi((a1 + m - 1)*h1, (b1 + j - 1)*h2, A1, A2, B1, B2, h1, h2) * 2/h1 + F((a1 + m - 1)*h1, (b1 + j - 1)*h2);
        }
    } else {
        for (int j = 2; j <= n-1; ++j) { 
            whatWriteTo[m][j] = F((a1 + m - 1)*h1, (b1 + j - 1)*h2);
        }
    }

    //let's look on the left border
    if (left <0) {
        for (int j = 2; j <= n-1; ++j) { 
            whatWriteTo[1][j] = psi((a1 + 1 - 1)*h1, (b1 + j - 1)*h2, A1, A2, B1, B2, h1, h2) * 2/h1 + F((a1 + 1 - 1)*h1, (b1 + j - 1)*h2);
        }
    } else {
        for (int j = 2; j <= n-1; ++j) { 
            whatWriteTo[1][j] = F((a1 + 1 - 1)*h1, (b1 + j - 1)*h2);
        }
    }

    // look on the down-left point
    if (down<0 && left<0) {
        // printf("Yes\n");
        // it's (11) equation
        whatWriteTo[1][1] = psi((a1 + 1 - 1)*h1, (b1 + 1 - 1)*h2, A1, A2, B1, B2, h1, h2) * (2/h1 + 2/h2) + F((a1 + 1 - 1)*h1, (b1 + 1 - 1)*h2);
    } else if (down<0 && left >=0 ) {
        whatWriteTo[1][1] = psi((a1 + 1 - 1)*h1, (b1 + 1 - 1)*h2, A1, A2, B1, B2, h1, h2) * 2/h2 + F((a1 + 1 - 1)*h1, (b1 + 1 - 1)*h2);
    } else if (down>=0 && left <0) {
        whatWriteTo[1][1] = psi((a1 + 1 - 1)*h1, (b1 + 1 - 1)*h2, A1, A2, B1, B2, h1, h2) * 2/h1 + F((a1 + 1 - 1)*h1, (b1 + 1 - 1)*h2);
    } else if (down>=0 && left >=0) {
        whatWriteTo[1][1] = F((a1 + 1 - 1)*h1, (b1 + 1 - 1)*h2);
    }


    // look on the bottom-right points
    if (down<0 && right<0) {
        // printf("Yes\n");
        // it's (12) equation
        whatWriteTo[m][1] = psi((a1 + m - 1)*h1, (b1 + 1 - 1)*h2, A1, A2, B1, B2, h1, h2) * (2/h1 + 2/h2) + F((a1 + m - 1)*h1, (b1 + 1 - 1)*h2);
    } else if (down<0 && right >=0 ) {
        whatWriteTo[m][1] = psi((a1 + m - 1)*h1, (b1 + 1 - 1)*h2, A1, A2, B1, B2, h1, h2) * 2/h2 + F((a1 + m - 1)*h1, (b1 + 1 - 1)*h2);
    } else if (down>=0 && right <0) {
        whatWriteTo[m][1] = psi((a1 + m - 1)*h1, (b1 + 1 - 1)*h2, A1, A2, B1, B2, h1, h2) * 2/h1 + F((a1 + m - 1)*h1, (b1 + 1 - 1)*h2);
    } else if (down>=0 && right >=0) {
        whatWriteTo[m][1] = F((a1 + m - 1)*h1, (b1 + 1 - 1)*h2);
    }

    // look on the top-right points
    if (up<0 && right<0) {
        // printf("Yes\n");
        // it's (13) equation
        whatWriteTo[m][n] = psi((a1 + m - 1)*h1, (b1 + n - 1)*h2, A1, A2, B1, B2, h1, h2) * (2/h1 + 2/h2) + F((a1 + m - 1)*h1, (b1 + n - 1)*h2);
    } else if (up<0 && right >=0 ) {
        whatWriteTo[m][n] = psi((a1 + m - 1)*h1, (b1 + n - 1)*h2, A1, A2, B1, B2, h1, h2) * 2/h2 + F((a1 + m - 1)*h1, (b1 + n - 1)*h2);
    } else if (up>=0 && right <0) {
        whatWriteTo[m][n] = psi((a1 + m - 1)*h1, (b1 + n - 1)*h2, A1, A2, B1, B2, h1, h2) * 2/h1 + F((a1 + m - 1)*h1, (b1 + n - 1)*h2);
    } else if (up>=0 && right >=0) {
        whatWriteTo[m][n] = F((a1 + m - 1)*h1, (b1 + n - 1)*h2);
    }

    // look on the top-left points
    if (up<0 && left<0) {
        // printf("Yes\n");
        // it's (14) equation
        whatWriteTo[1][n] = psi((a1 + 1 - 1)*h1, (b1 + n - 1)*h2, A1, A2, B1, B2, h1, h2) * (2/h1 + 2/h2) + F((a1 + 1 - 1)*h1, (b1 + n - 1)*h2);
    } else if (up<0 && left >=0 ) {
        whatWriteTo[1][n] = psi((a1 + 1 - 1)*h1, (b1 + n - 1)*h2, A1, A2, B1, B2, h1, h2) * 2/h2 + F((a1 + 1 - 1)*h1, (b1 + n - 1)*h2);
    } else if (up>=0 && left <0) {
        whatWriteTo[1][n] = psi((a1 + 1 - 1)*h1, (b1 + n - 1)*h2, A1, A2, B1, B2, h1, h2) * 2/h1 + F((a1 + 1 - 1)*h1, (b1 + n - 1)*h2);
    } else if (up>=0 && left >=0) {
        whatWriteTo[1][n] = F((a1 + 1 - 1)*h1, (b1 + n - 1)*h2);
    }
}

void minus(double** first, double** second, double** whatWriteTo, double M, double N, Info_t* info) {
    int m = info->m;
    int n = info->n;
    for (size_t i = 1; i <= m; ++i) {
            for (size_t j = 1; j <= n; ++j) {
            whatWriteTo[i][j] = first[i][j] - second[i][j];
        }
    }
}

double ro(int index, int M_or_N) {
    if (index == 0 || index == M_or_N) {
        return 0.5;
    } else {
        return 1.0;
    }
}

double scalarProduct(double** first, double** second, double M, double N, double h1, double h2, Info_t* info, MPI_Comm* Comm) {
    int m = info->m;
    int n = info->n;
    int a1 = info->a1;
    int b1 = info->b1;
    double local_sum = 0.0;
    double reduced_sum = 0.0;
    for (size_t i = 1; i <= m; ++i) {
        for (size_t j = 1; j <= n; ++j) {
            local_sum = local_sum + h1*h2*ro((a1 + i - 1), M)*ro((b1 + j - 1), N) * first[i][j] * second[i][j];
        }
    }
    MPI_Allreduce(&local_sum, &reduced_sum, 1, MPI_DOUBLE, MPI_SUM, *Comm); 
    return reduced_sum;
}

double getMaxNorm(double** items, double M, double N, double h1, double h2, Info_t* info, MPI_Comm* Comm) {
    int m = info->m;
    int n = info->n;
    int a1 = info->a1;
    int b1 = info->b1;
    double local_max = 0.0;
    double reduced_max = 0.0;
    for (size_t i = 1; i <= m; ++i) {
        for (size_t j = 1; j <= n; ++j) {
            double item = fabs(items[i][j]);
            if (item > local_max) {
                local_max = item;
            }
        }
    }
    // printf("local_max: %.15f \n", local_max);
    MPI_Allreduce(&local_max, &reduced_max, 1, MPI_DOUBLE, MPI_MAX, *Comm); 
    return reduced_max;
}

void multiplyByNum(double** items, double num, double** whatWriteTo, double M, double N, Info_t* info) {
    int m = info->m;
    int n = info->n;
    for (size_t i = 1; i <= m; ++i) {
            for (size_t j = 1; j <= n; ++j) {
            whatWriteTo[i][j] = items[i][j]*num;
        }
    }
}

void sendrecv(double **domain, 
              double *send_up_row, double *recv_up_row, 
              double *send_down_row, double *recv_down_row,
              double *send_left_column, double *recv_left_column,
              double *send_right_column, double *recv_right_column,
              MPI_Comm* Comm, Info_t *info) {
    
    int up = info->up_rank;
    int down = info->down_rank;
    int left = info->left_rank;
    int right = info->right_rank;
    int m = info->m;
    int n = info->n;
    int rank = info->rank;
    for (size_t i = 0; i < info->m; ++i) {
      send_down_row[i] = domain[i + 1][1];
      send_up_row[i] = domain[i + 1][n]; 
    }
    for (size_t j = 0; j < n; ++j) {  
      send_left_column[j] = domain[1][j + 1]; 
      send_right_column[j] = domain[m][j + 1];
    }
    MPI_Status Status;

    
    

	if ((up < 0) && (down >= 0)) {
        
        // printf("Wait1");
		MPI_Sendrecv(send_down_row, m,MPI_DOUBLE, down, 0, recv_down_row, m, MPI_DOUBLE, down, 0, *Comm, &Status);
        // printf("End1");
	}
    else if ((up >= 0) && (down < 0)) {
        // printf("Wait3");
		MPI_Sendrecv(send_up_row, m, MPI_DOUBLE, up, 0, recv_up_row, m, MPI_DOUBLE, up, 0, *Comm, &Status);
        // printf("End3");
	}
	else if ((up >= 0) && (down >= 0)) {
        // printf("Wait2");
		MPI_Sendrecv(send_up_row, m, MPI_DOUBLE, up, 0, recv_up_row, m, MPI_DOUBLE, up, 0, *Comm, &Status);
		MPI_Sendrecv(send_down_row, m, MPI_DOUBLE, down, 0, recv_down_row, m, MPI_DOUBLE, down, 0, *Comm, &Status);
        // printf("End2");
	}
	


	if ((left < 0) && (right >= 0)) {
		MPI_Sendrecv(send_right_column, n, MPI_DOUBLE, right, 0, recv_right_column, n, MPI_DOUBLE, right, 0, *Comm, &Status);
	}
	else if ((left >= 0) && (right >= 0)) {
		MPI_Sendrecv(send_left_column, n, MPI_DOUBLE, left, 0, recv_left_column, n, MPI_DOUBLE, left, 0, *Comm, &Status);
		MPI_Sendrecv(send_right_column, n, MPI_DOUBLE, right, 0, recv_right_column, n, MPI_DOUBLE, right, 0, *Comm, &Status);

	}
	else if ((left >= 0) && (right < 0)) {
		MPI_Sendrecv(send_left_column, n, MPI_DOUBLE, left, 0, recv_left_column, n, MPI_DOUBLE, left, 0, *Comm, &Status);
	}
    
    // printf("I'm okay %d\n", info->rank);

    for (size_t i = 0; i < m; ++i) {
        domain[i + 1][0] = recv_down_row[i];
        domain[i + 1][n + 1] = recv_up_row[i];
    }
    for (size_t j = 0; j < n; ++j) {
        domain[0][j + 1] = recv_left_column[j];
        domain[m + 1][j + 1] = recv_right_column[j];
    }
}

void getAnalyticalSolution(double** whatWriteTo, double h1, double h2, Info_t* info) {
    int m = info->m;
    int n = info->n;
    int a1 = info->a1;
    int b1 = info->b1;
    for (int i=1; i <= m; ++i) {
        for (int j=1; j <= n; ++j) {
            whatWriteTo[i][j] = u((a1 + i - 1)*h1, (b1 + j - 1)*h2);
        }
    }
}

void solving (double h1, double h2, double epsilon, double A1, double A2, double B1, double B2, int M, int N, Info_t* info, MPI_Comm *Comm, int size) {
    double start_time = MPI_Wtime();
    int m = info->m;
    int n = info->n;
    int rank = info->rank;
    

    double** omega              = (double**)malloc((m + 2) * sizeof(double*)); 
    double** omega_next         = (double**)malloc((m + 2) * sizeof(double*)); 
    double** B                  = (double**)malloc((m + 2) * sizeof(double*)); 
    double** A_omega            = (double**)malloc((m + 2) * sizeof(double*)); 
    double** r                  = (double**)malloc((m + 2) * sizeof(double*)); 
    double** A_r                = (double**)malloc((m + 2) * sizeof(double*)); 
    double** tau_r              = (double**)malloc((m + 2) * sizeof(double*)); 
    double** difference_omega   = (double**)malloc((m + 2) * sizeof(double*));
    double** solution           = (double**)malloc((m + 2) * sizeof(double*)); 

    for (size_t i = 0; i <= m + 1; ++i) {
        omega[i]            = (double*)malloc((n + 2) * sizeof(double));
        omega_next[i]       = (double*)malloc((n + 2) * sizeof(double));
        B[i]                = (double*)malloc((n + 2) * sizeof(double));
        A_omega[i]          = (double*)malloc((n + 2) * sizeof(double));
        r[i]                = (double*)malloc((n + 2) * sizeof(double));
        A_r[i]              = (double*)malloc((n + 2) * sizeof(double));
        tau_r[i]            = (double*)malloc((n + 2) * sizeof(double));
        difference_omega[i] = (double*)malloc((n + 2) * sizeof(double));
        solution[i]         = (double*)malloc((n + 2) * sizeof(double));
    }

    double tau = 0.0;
    double difference_local = epsilon;
    double difference_global = epsilon;
    for (size_t i = 0; i <= m + 1; ++i) {
        for (size_t j = 0; j <= n + 1; ++j) {
            omega[i][j] = 0.0;
            omega_next[i][j] = 2.0;
            B[i][j] = 0.0;
            A_omega[i][j] = 0.0;
            r[i][j] = 0.0;
            A_r[i][j] = 0.0;
            tau_r[i][j] = 0.0;
            difference_omega[i][j] = 0.0;
            solution[i][j] = 0.0;
        }
    }
    
    
    getB(B, M, N, h1, h2, A1, A2, B1, B2, info);
    double *send_up_row =       (double*) malloc(m * sizeof(double));
    double *recv_up_row =       (double*) malloc(m * sizeof(double));
    double *send_down_row =     (double*) malloc(m * sizeof(double));
    double *recv_down_row =     (double*) malloc(m * sizeof(double));
    double *send_left_column =  (double*) malloc(n * sizeof(double));
    double *recv_left_column =  (double*) malloc(n * sizeof(double));
    double *send_right_column = (double*) malloc(n * sizeof(double));
    double *recv_right_column = (double*) malloc(n * sizeof(double));
    
    int count = 0;
    while (difference_global >= epsilon && count < 50)
    {
        
        for (size_t i = 1; i <= m; ++i) {
            for (size_t j = 1; j <= n; ++j) {
                omega[i][j] = omega_next[i][j];
            }
        }

        sendrecv(omega, 
             send_up_row, recv_up_row, 
             send_down_row, recv_down_row, 
             send_left_column, recv_left_column, 
             send_right_column, recv_right_column,
             Comm, info);
        
        applyA(omega, A_omega, M, N, h1, h2, A1, A2, B1, B2, info);
        minus(A_omega, B, r, M, N, info);
        sendrecv(r, 
             send_up_row, recv_up_row, 
             send_down_row, recv_down_row, 
             send_left_column, recv_left_column, 
             send_right_column, recv_right_column,
             Comm, info);
        applyA(r, A_r, M, N, h1, h2, A1, A2, B1, B2, info);
        tau = scalarProduct(A_r, r, M, N, h1, h2, info, Comm) / scalarProduct(A_r, A_r, M, N, h1, h2, info, Comm);
        multiplyByNum(r, tau, tau_r, M, N, info);
        minus(omega, tau_r, omega_next, M, N, info);
        difference_local = sqrt(scalarProduct(tau_r, tau_r, M, N, h1, h2, info, Comm));
        MPI_Allreduce(&difference_local, &difference_global, 1, MPI_DOUBLE, MPI_MAX, *Comm);
        if (rank == 0 && count % 1 == 0) {
            printf("n: %d, rank: %d, norm: %.10f, eps:%.8f\n", count, rank, difference_global, epsilon);
        }
        count++;

    }

    double local_time_diff = MPI_Wtime() - start_time;
    
    double global_time_diff = 0.0;
    getAnalyticalSolution(solution, h1, h2, info);
    minus(solution, omega_next, solution, M, N, info);
    double norm = getMaxNorm(solution, M, N, h1, h2, info, Comm);
    MPI_Allreduce(&local_time_diff, &global_time_diff, 1, MPI_DOUBLE, MPI_MAX, *Comm);
    
    // return;
    if (info->rank == 0) {
        printf("size ,  M , N   , time    , diff_u_w,  eps\n");
        printf("%d, %d, %d, %.5f, %.10f, %.8f\n", info->size, M, N, global_time_diff, norm, epsilon);
    }
    
}


int main(int argc, char** argv) {

    const size_t M = 16000;
    const size_t N = 16000;
    
    double epsilon = 0.00000001;    
    
    double A1 = 0.0;
    double A2 = 4.0;
    double B1 = 0.0;
    double B2 = 3.0;
    double h1 = (A2 - A1)/(double)M;
    double h2 = (B2 - B1)/(double)N;
    Info_t info;
    MPI_Comm Comm;
    MPI_Init(&argc, &argv);
    
    
    int rank, size;
    partitioningDomain(M, N, &Comm, rank, size, &info);
    
    solving(h1, h2, epsilon, A1, A2, B1, B2, M, N, &info, &Comm, size);
    
    MPI_Finalize();
}


// mpicc parallelMPI.c -o parallel && mpiexec -np 2 ./parallel
// mpicc parallelMPI.cpp -o parallel.o && mpiexec -np 4 ./parallel.o 20 20 73.51