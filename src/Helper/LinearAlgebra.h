#ifndef CARIBOU_HELPER_LINEARALGEBRA_H
#define CARIBOU_HELPER_LINEARALGEBRA_H

#include <cmath>
#include <sofa/defaulttype/Mat.h>
#include <sofa/defaulttype/Vec.h>
#include <cfloat>

namespace sofa {

namespace caribou {

namespace helper {

namespace linear_algebra {

/// Jacobi rotate from Irving et al.
/// \tparam Matrix3 A 3x3 matrix type that define the accessor A[l][c] for line 'l' and column 'c'
/// \param A
/// \param R
/// \param p
/// \param q
template<class real>
inline void jacobiRotate(defaulttype::Mat<3,3,real> & A, defaulttype::Mat<3,3,real> & R, unsigned int p, unsigned int q) {
    if (A[p][q] == 0)
        return;

    double d = (A[p][p] - A[q][q]) / (2.0*A[p][q]);
    double t = 1.0 / (fabs(d) + sqrt(d*d + 1.0));
    if (d < 0.0)
        t = -t;
    double c = 1.0 / sqrt(t*t + 1);
    double s = t*c;
    A[p][p] += t*A[p][q];
    A[q][q] -= t*A[p][q];
    A[p][q] = A[q][p] = 0.0;
    unsigned int k;
    for (k = 0; k < 3; k++)
    {
        if (k != p && k != q)
        {
            double Akp = c*A[k][p] + s*A[k][q];
            double Akq = -s*A[k][p] + c*A[k][q];
            A[k][p] = A[p][k] = Akp;
            A[k][q] = A[q][k] = Akq;
        }
    }
    for (k = 0; k < 3; k++)
    {
        double Rkp = c*R[k][p] + s*R[k][q];
        double Rkq = -s*R[k][p] + c*R[k][q];
        R[k][p] = Rkp;
        R[k][q] = Rkq;
    }
}

template<class real>
void eigenDecomposition(const defaulttype::Mat<3,3,real> &A, defaulttype::Mat<3,3,real> &eigenVecs, defaulttype::Vec<3,real> &eigenVals)
{
    const int numJacobiIterations = 10;
    const double epsilon = 1e-15;
    defaulttype::Mat<3,3,real> D = A;
    eigenVecs = {{1, 0, 0}, {0, 1, 0}, {0, 0, 1}};
    int iter = 0;
    while (iter < numJacobiIterations)
    {
        unsigned int p, q;
        double a, max;
        max = fabs(D[0][1]);
        p = 0;
        q = 1;
        a = fabs(D[0][2]);
        if (a > max)
        {
            p = 0;
            q = 2;
            max = a;
        }
        a = fabs(D[1][2]);
        if (a > max)
        {
            p = 1;
            q = 2;
            max = a;
        }
        if (max < epsilon)
            break;
        jacobiRotate(D, eigenVecs, p, q);
        iter++;
    }
    eigenVals[0] = D[0][0];
    eigenVals[1] = D[1][1];
    eigenVals[2] = D[2][2];
}

template<class real>
inline defaulttype::Mat<3,3,real> rotationMatrixIrving(const defaulttype::Mat<3,3,real> &A) {
    defaulttype::Mat<3,3,real> AT_A, V;
    AT_A = A.transposed() * A;
    defaulttype::Vec<3,real> S;
    eigenDecomposition(AT_A, V, S);
    const double detV = determinant(V);
    if (detV < 0.0)
    {
        double minLambda = DBL_MAX;
        unsigned char pos = 0;
        for (unsigned char l = 0; l < 3; l++)
        {
            if (S[l] < minLambda)
            {
                pos = l;
                minLambda = S[l];
            }
        }
        V(0, pos) = -V(0, pos);
        V(1, pos) = -V(1, pos);
        V(2, pos) = -V(2, pos);
    }
    if (S[0] < 0.0f)
        S[0] = 0.0f;
    if (S[1] < 0.0f)
        S[1] = 0.0f;
    if (S[2] < 0.0f)
        S[2] = 0.0f;
    defaulttype::Vec<3,real>
        sigma;
    sigma[0]
        = sqrt(S[0]);
    sigma[1]
        = sqrt(S[1]);
    sigma[2]
        = sqrt(S[2]);
    unsigned char chk = 0;
    unsigned char pos = 0;
    defaulttype::Mat<3,3,real> U;
    for (unsigned char l = 0; l < 3; l++)
    {
        if (fabs(sigma[l]) < 1.0e-4)
        {
            pos = l;
            chk++;
        }
    }
    if (chk > 0)
    {
        if (chk > 1)
        {
            U.identity();
        }
        else
        {
            U = A * V;
            for (unsigned char l = 0; l < 3; l++)
            {
                if (l != pos)
                {
                    for (unsigned char m = 0; m < 3; m++)
                    {
                        U(m, l) *= 1.0f / sigma[l];
                    }
                }
            }
            defaulttype::Vec<3,real> v[2];
            unsigned char index = 0;
            for (unsigned char l = 0; l < 3; l++)
            {
                if (l != pos)
                {
                    v[index++] = defaulttype::Vec<3,real>(U(0, l), U(1, l), U
                        (2, l));
                }
            }
            defaulttype::Vec<3,real> vec = v[0].cross(v[1]);
            vec.normalize();
            U(0, pos) = vec[0];
            U(1, pos) = vec[1];
            U(2, pos) = vec[2];
        }
    }
    else
    {
        defaulttype::Vec<3,real> sigmaInv(1.0 / sigma[0], 1.0 / sigma[1],
                          1.0 / sigma[2]);
        U = A * V;
        for (unsigned char l = 0; l < 3; l++)
        {
            for (unsigned char m = 0; m < 3; m++)
            {
                U(m, l) *= sigmaInv[l];
            }
        }
    }
    const double detU = determinant(U);
    if (detU < 0.0)
    {
        double minLambda = DBL_MAX;
        pos = 0;
        for (unsigned char l = 0; l < 3; l++)
        {
            if (sigma[l] < minLambda)
            {
                pos = l;
                minLambda = sigma[l];
            }
        }
        sigma[pos] = -sigma[pos];
        U(0, pos) = -U(0, pos);
        U(1, pos) = -U(1, pos);
        U(2, pos) = -U(2, pos);
    }

    return U * V.transposed();
}



} // namespace triangle

} // namespace helper

} // namespace caribou

} // namespace sofa

#endif //CARIBOU_HELPER_LINEARALGEBRA_H
