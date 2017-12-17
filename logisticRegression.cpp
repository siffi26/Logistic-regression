/*Program to perform Logistic regression*/
/*Author: Siffi Singh */
/*Dated: 09/10/2017 */

/* standard Header */
#include <bits/stdc++.h>
#include <random>
#include <cmath>
using namespace std;
# define two_pi 2.0*3.14159265358979323846
/*Generating guassian distributed data*/
double guassian(double mu, double sigma) {
        double z0, z1;
        double epsilon = std::numeric_limits < double > ::min();
        bool generate;
        generate = !generate;
        if (!generate)
            return z1 * sigma + mu;
        double u1, u2;
        do {
            u1 = rand() * (1.0 / RAND_MAX);
            u2 = rand() * (1.0 / RAND_MAX);
        } while (u1 <= epsilon);
        z0 = sqrt(-2.0 * log(u1)) * cos(two_pi * u2);
        z1 = sqrt(-2.0 * log(u1)) * sin(two_pi * u2);
        return z0 * sigma + mu;
    }
    /*Matrix multiplication*/
void matrixmultiply(int n, int m, int p, int q, double A[10][10], double B[10][10], double res[10][10]) {
        for (int i = 0; i < n; i++) {
            for (int j = 0; j < q; j++) {
                res[i][j] = 0;
            }
        }
        for (int i = 0; i < n; i++) //no. of rows in the first matrix
        {
            for (int j = 0; j < q; j++) //no. of columns in the second matrix
            {
                for (int k = 0; k < p; k++) //q=m the dimension that is equal
                {
                    res[i][j] += A[i][k] * B[k][j];
                }
            }
        }
    }
    /*Matrix Inverse*/
void inverse(int n, double b[10][10], double a[10][10]) {
        for (int i = 0; i < n; i++) {
            for (int j = 0; j < n; j++) {
                a[i][j] = b[i][j];
            }
        }
        for (int i = 0; i < n; i++) {
            for (int j = n; j < 2 * n; j++) {
                if (i == j - n)
                    a[i][j] = 1;
                else
                    a[i][j] = 0;
            }
        }
        for (int i = 0; i < n; i++) {
            double t = a[i][i];
            for (int j = i; j < 2 * n; j++)
                a[i][j] = a[i][j] / t;
            for (int j = 0; j < n; j++) {
                if (i != j) {
                    t = a[j][i];
                    for (int k = 0; k < 2 * n; k++)
                        a[j][k] = a[j][k] - t * a[i][k];
                }
            }
        }
    }
    /*Single matrix multiply*/
double matrixmul(int n, double X[][10], double W[][10]) {
        double ans = 0;
        for (int i = 0; i < 1; i++) {
            for (int j = 0; j < 1; j++) {
                for (int k = 0; k < n; k++) {
                    ans = ans + X[i][k] * W[k][j];
                }
            }
        }
        return ans;
    }
    /*Matrix Transpose*/
void transpose(int n, int m, double a[10][10], double b[10][10]) {
        for (int i = 0; i < n; i++) {
            for (int j = 0; j < m; j++) {
                b[i][j] = a[j][i];
            }
        }
    }
    /*Determinant of Matrix*/
double d = 0;
double det(int n, double mat[10][10]) {
        int c, subi, i, j, subj;
        double submat[10][10];
        if (n == 2) {
            return ((mat[0][0] * mat[1][1]) - (mat[1][0] * mat[0][1]));
        } else {
            for (c = 0; c < n; c++) {
                subi = 0;
                for (i = 1; i < n; i++) {
                    subj = 0;
                    for (j = 0; j < n; j++) {
                        if (j == c) {
                            continue;
                        }
                        submat[subi][subj] = mat[i][j];
                        subj++;
                    }
                    subi++;
                }
                d = d + (pow(-1, c) * mat[0][c] * det(n - 1, submat));
            }
        }
        return d;
    }
    /*Driver Function*/
int main() {
    /*Variable Declarations*/
    double mx1, vx1, mx2, vx2, my1, vy1, my2, vy2;
    int n;
    double x1[10], y1[10], x2[10], y2[10], z1[10], z2[10];
    map < double, double > map1, map2;
    cout << "Enter mean and variance values for Data 1, X: ";
    cin >> mx1 >> vx1;
    cout << "Enter mean and variance values for Data 1, Y: ";
    cin >> my1 >> vy1;
    cout << "Enter mean and variance values for Data 2, X: ";
    cin >> mx2 >> vx2;
    cout << "Enter mean and variance values for Data 2, Y: ";
    cin >> my2 >> vy2;
    cout << "Enter the no. of values you want to generate in D1 and D2: ";
    cin >> n;
    for (int i = 0; i < n; i++) {
        double gen1, gen2, gen3, gen4;
        gen1 = guassian(mx1, vx1);
        gen2 = guassian(my1, vy1);
        map1[gen1] = gen2;
        x1[i] = gen1; //'x' value
        y1[i] = gen2; //'Y' value
        z1[i] = 0; //'label'
        gen3 = guassian(mx2, vx2);
        gen4 = guassian(mx2, vx2);
        map2[gen3] = gen4;
        x2[i] = gen3; //'x' value
        y2[i] = gen4; //'Y' value
        z2[i] = 1; //'label'
    }
    /* Shuffling random_shuffle(v.begin(), v.end(),myrandom);*/
    double data[10][3];
    for (int i = 0; i < 2 * n - 1; i++) {
        if (i % 2 == 0) {
            data[i][0] = x1[i];
            data[i][1] = y1[i];
            data[i][2] = z1[i];
        } else {
            data[i][0] = x2[i];
            data[i][1] = y2[i];
            data[i][2] = z2[i];
        }
    }
    int n1 = n;
    n = 2 * n - 1;
    int basis;
    /*Basis number*/
    cout << "Enter basis number: ";
    cin >> basis;
    /*Entering weight values*/
    cout << "Enter weight values: ";
    double W[10][10], WT[10][10];
    for (int i = 0; i < basis; i++) {
        cin >> W[i][0];
        WT[0][i] = W[i][0];
    }
    /*making design matrix*/
    for (int i = 0; i < n; i++) {
        double X[10][10], XT[10][10];
        for (int j = 0; j <= basis; j++) {
            X[j][0] = pow(data[i][0], j);
            XT[0][j] = pow(data[i][0], j);
        }
        double WtX = matrixmul(basis, WT, X);
        double R[10][10];
        double y = (1 / (1 + exp(WtX)));
        /*Forming 'R' vector*/
        R[0][0] = y * (1 - y);
        double XTR[10][10] = {
            0
        }, Hessian[10][10] = {
            0
        };
        /*Computing Hessian matrix*/
        matrixmultiply(basis, 1, 1, 1, XT, R, XTR);
        //		for(int j=0; j<n; j++)
        //		{
        //			for(int k=0; k<1; k++)
        //			cout<<XTR[i]<<" ";
        //		}
        matrixmultiply(basis, 1, 1, basis, XTR, X, Hessian);
        double detr;
        detr = det(basis, Hessian); //Computing determinant of Hessian matrix
        double w_new[10][10];
        /*if determinant of hessian matrix is not equal to zero, then we use the formula w_new = w_old + H* (transpose(x)) * (y-t) where y= (1/(1+exp(wT*x))) */
        if (detr != 0) {
            double Res[10][10];
            matrixmultiply(basis, basis, basis, 1, Hessian, X, Res);
            for (int m = 0; m < basis; m++) {
                for (int s = 0; s < 1; s++) {
                    Res[m][s] = Res[m][s] * (y - data[m][2]);
                    w_new[m][s] = W[m][s] + Res[m][s];
                }
            }
        }
        /*else if determinant of hessian matrix is equal to zero, then we use the formula w_new = w_old + (y-G(x,w))* where G(x,w)= (1/(1+exp(wT*x))) */
        else {
            for (int m = 0; m < basis; m++) {
                for (int s = 0; s < 1; s++) {
                    w_new[m][s] = W[m][s] + (y - data[m][2]);
                }
            }
        }
        /*Checking convergence*/
        int flag = 1;
        for (int j = 0; j < basis; j++) {
            if (abs((w_new[j] - W[j])) <= (10 ^ -5))
                flag = 0;
        }
        if (flag == 1) {
            cout << "Output converged at " << i + 1 << "th input!" << endl;	
            break;
        } else
            continue;
    }
    return 0;
}


