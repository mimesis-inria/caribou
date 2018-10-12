#ifndef CARIBOU_MIRTICHINTEGRATION_H
#define CARIBOU_MIRTICHINTEGRATION_H
#include <iostream>
#include <cmath>
#include <vector>
#include <array>

void compProj( std::vector< std::vector< double > > &projInt, const std::vector< double > &face0, const std::vector< double > &face1,
               const std::vector < std::vector< unsigned int > > &exponentsAlphaBeta, const std::vector< unsigned int > &factorials)
{
    double delta;
    size_t idx;
    size_t p,q;

    for( size_t num = 0; num < exponentsAlphaBeta[0].size(); ++num )
    {

        if( exponentsAlphaBeta[1][num] > exponentsAlphaBeta[0][num] )
        {
            idx = 1;
            p = exponentsAlphaBeta[0][num];
            q = exponentsAlphaBeta[1][num] + 1;
        }
        else
        {
            idx = 0;
            p = exponentsAlphaBeta[0][num] + 1;
            q = exponentsAlphaBeta[1][num];
        }

        size_t k = face0.size();
        for( size_t e = 0; e < k; ++e )
        {
            if( idx == 0 )
            {
                delta = face1[(e+1)%k] - face1[e];
            }
            else
            {
                delta = face0[(e+1)%k] - face0[e];
            }

            for( size_t i = 0; i < p+1; ++i )
            {
                for( size_t j = 0; j < q+1; ++j )
                {
                    projInt[exponentsAlphaBeta[0][num]][exponentsAlphaBeta[1][num]] += delta*
                                                                                       ((factorials[p]/(double)(factorials[p-i]*factorials[i]))*(factorials[q]/(double)(factorials[q-j]*factorials[j]))/(factorials[p+q]/(double)(factorials[p+q-(i+j)]*factorials[i+j])))*
                                                                                       pow(face0[(e+1)%k],i) * pow(face0[e],p-i) * pow(face1[(e+1)%k],j) * pow(face1[e],q-j);
                }
            }
        }
        projInt[exponentsAlphaBeta[0][num]][exponentsAlphaBeta[1][num]] *= pow(-1.,idx)/(double)((exponentsAlphaBeta[idx][num]+1)*(p+q+1));
    }
}

void compFace( std::vector< double > &faceIntegrals, const std::vector< double > &face0, const std::vector< double > &face1, const std::vector< double > &face2, const double (&normal)[3],
               const std::vector< unsigned int > &exponents0, const std::vector< unsigned int > &exponents1, const std::vector< unsigned int > &exponents2, const std::vector < std::vector< unsigned int > > &exponentsAlphaBeta, const std::vector< unsigned int > &factorials )
{
    std::vector< std::vector< double > > projInt( 9, std::vector< double >(9, 0.));

    compProj( projInt, face0, face1, exponentsAlphaBeta, factorials);

    double w = -1.*(normal[0]*face0[1] + normal[1]*face1[1] + normal[2]*face2[1]); //fixes problems if normal not normed.

    for( size_t i = 0; i < exponents0.size(); ++i)
    {
        for( size_t a = 0; a < exponents2[i]+1; ++a)
        {
            for( size_t b = 0; b < exponents2[i]+1 - a; ++b)
            {
                faceIntegrals[i] += 1./(double)(factorials[a]*factorials[b]*factorials[exponents2[i]-a-b])*pow(normal[0],a)*pow(normal[1],b)*
                                    pow(w,exponents2[i]-a-b)*projInt[a + exponents0[i]][ b + exponents1[i]];
            }
        }
        faceIntegrals[i] *= pow(normal[2], -exponents2[i]-1)*pow(-1.,exponents2[i])*(double)factorials[exponents2[i]];
        //if( normal[2] < 0. )
        //faceIntegrals[i] *= (-1.); ///is this right?
    }
}

template<class Coord>
void compVol(
        std::vector< double > &volInt,
        std::vector< std::vector< double > > &bdryInt,
        const std::vector<  unsigned int > &cutElement,
        /*const vector< vector< int > > &bdryElements, */
        const std::vector< std::vector< unsigned int > > &faces,
        const std::vector<Coord> &coords,
        const std::vector< std::vector< unsigned int > > &exponents,
        /*const vector< vector< int > > &exponentsBdry, */
        const std::vector< std::vector< unsigned int > > &exponentsAlphaBeta,
        /*const vector< vector< int > > &exponentsAlphaBetaBdry, */
        const std::vector< unsigned int > &factorials,
        std::vector< double > &normal )
{


    std::vector< std::vector< unsigned int > > exponentsAfterDivThm = exponents;
    std::vector< int > indexMinExp (exponents[0].size(),0);

    //double normal[3];
    double permutedNormal[3];

    for( unsigned int j = 0; j < exponents[0].size(); ++j )
    {
        if( exponents[1][j] > 0 )
        {
            if( exponents[1][j] < exponents[indexMinExp[j]][j] )
                indexMinExp[j] = 1;
        }
        if( exponents[2][j] > 0 )
        {
            if( exponents[2][j] < exponents[indexMinExp[j]][j] )
                indexMinExp[j] = 2;
        }
        exponentsAfterDivThm[indexMinExp[j]][j] += 1;
    }


    //Permutation to  ensure reasonable projections onto coordinate planes
    int right_handed_permutation[3][3] = { {1,2,0}, {2,0,1}, {0,1,2} };

    size_t nFaces = cutElement.size();
    //cout << "nFaces = " << nFaces << "\n";

    std::fill(volInt.begin(), volInt.end(),0.);

    for( size_t fac = 0; fac < nFaces; ++fac )
    {
        // DEBUG;
        size_t face = cutElement[fac];
        size_t nPointsFaces = faces[face].size();
        //std::cout << nPointsFaces << "\n";

        // DEBUG;
        std::fill(bdryInt[fac].begin(), bdryInt[fac].end(),0.);

        //current Face coordinates
        std::vector< std::vector< double > > currentFace;
//        currentFace.assign(std::vector < double > (nPointsFaces));
        currentFace.resize(3);
        for( int dir = 0; dir < 3; ++dir )
            currentFace[dir].resize(nPointsFaces);
        for( size_t i = 0; i < nPointsFaces; ++i )
        {
            currentFace[0][i] = coords[faces[face][i]][0];          //transpose here to use the subprograms
            currentFace[1][i] = coords[faces[face][i]][1];
            currentFace[2][i] = coords[faces[face][i]][2];
        }

        //compute normal
        //Needs counterclockwise indexing in right-handed coordinate system to produce the correct normal

        /// CAREFUL! NORMAL MAY BE TOO SMALL. IF ALL PARTS SMALL -> ignore the element!
        double norm;

        //Compute the normal; not enough to take 3 neighboring points, face does not have to be convex!! (NEWELL's method)

        normal[0] = normal[1] = normal[2] = 0.;
        for( size_t h = 0; h < nPointsFaces; ++h )
        {
            normal[0] += (currentFace[1][h] - currentFace[1][(h+1)%nPointsFaces])*(currentFace[2][h] + currentFace[2][(h+1)%nPointsFaces]);
            normal[1] += (currentFace[2][h] - currentFace[2][(h+1)%nPointsFaces])*(currentFace[0][h] + currentFace[0][(h+1)%nPointsFaces]);
            normal[2] += (currentFace[0][h] - currentFace[0][(h+1)%nPointsFaces])*(currentFace[1][h] + currentFace[1][(h+1)%nPointsFaces]);
        }

        int normalBigEnough = 0;

        norm = sqrt(normal[0]*normal[0] + normal[1]*normal[1] + normal[2]*normal[2]);
        //SHOW(norm);
        if( norm > 1e-12) //was too big!
        {
            normalBigEnough = 1;
//           DEBUG;
        }
//        else
//            DEBUG;
        if( normalBigEnough )
        {
            normal[0] /= norm; normal[1] /= norm; normal[2] /= norm;//normal not normed! VERY IMPORTANT FOR THE BOUNDARY MATRICES!

            int idx = 0;
            if( std::abs(normal[1]) > std::abs(normal[idx]) )
                idx = 1;
            if( std::abs(normal[2]) > std::abs(normal[idx]) )
                idx = 2;

            permutedNormal[0] = normal[right_handed_permutation[idx][0]];
            permutedNormal[1] = normal[right_handed_permutation[idx][1]];
            permutedNormal[2] = normal[right_handed_permutation[idx][2]];


            std::vector< double > faceIntegrals(84,0.);
            std::vector< double > faceIntegrals2(84,0.);

            compFace( faceIntegrals,
                      currentFace[right_handed_permutation[idx][0]],
                      currentFace[right_handed_permutation[idx][1]],
                      currentFace[right_handed_permutation[idx][2]],
                      permutedNormal,
                      exponentsAfterDivThm[right_handed_permutation[idx][0]],
                      exponentsAfterDivThm[right_handed_permutation[idx][1]],
                      exponentsAfterDivThm[right_handed_permutation[idx][2]],
                      exponentsAlphaBeta, factorials );

            //the second compFace is only necessary for boundary faces, can be optimized!
            compFace( faceIntegrals2,
                      currentFace[right_handed_permutation[idx][0]],
                      currentFace[right_handed_permutation[idx][1]],
                      currentFace[right_handed_permutation[idx][2]],
                      permutedNormal,
                      exponents[right_handed_permutation[idx][0]],
                      exponents[right_handed_permutation[idx][1]],
                      exponents[right_handed_permutation[idx][2]],
                      exponentsAlphaBeta, factorials ); // for boundary faces

            for( int i = 0; i < 84; ++i )
            {
                volInt[i] += 1./(double)(exponentsAfterDivThm[indexMinExp[i]][i])*normal[indexMinExp[i]]*faceIntegrals[i];
                bdryInt[fac][i] = faceIntegrals2[i];
                //cout << faceIntegrals[i] << ", ";//to compare with MATLAB code
            }
        }
    }
}

template<class Coord>
std::vector<double> integrate(const std::vector<std::vector<Coord>> & faces,
                              const std::array<std::vector<unsigned int>, 3> & exponents)
{
    static const std::vector<unsigned int> factorials = {
            1,
            1,
            2,
            6,
            24,
            120,
            720,
            5040,
            40320,
            362880
    };

    std::vector<std::vector<unsigned int> > exponentsAlphaBeta(2, std::vector< unsigned int > (45));
    std::size_t counter = 0;
    for( unsigned int degA = 0; degA < 9; ++degA )
    {
        for( unsigned int degB = 0; degB < 9 - degA; ++degB )
        {
            exponentsAlphaBeta[0][counter] = degA;
            exponentsAlphaBeta[1][counter] = degB;
            ++counter;
        }
    }

    std::vector<  unsigned int > faces_indices(faces.size()); // Assign indices to each faces
    std::vector< std::vector< unsigned int > > faces_nodes_indices(faces.size()); // Assign indices to each node in a face
    std::vector<Coord> coords;
    std::vector< std::vector< unsigned int > > _exponents(3);
    for (unsigned int e : exponents[0]) _exponents[0].push_back((int) e);
    for (unsigned int e : exponents[1]) _exponents[1].push_back((int) e);
    for (unsigned int e : exponents[2]) _exponents[2].push_back((int) e);

    unsigned int node_counter = 0;
    for (unsigned int i = 0; i < faces.size(); ++i) {
        faces_indices[i] = i;
        const std::vector<Coord> & positions = faces[i];
        faces_nodes_indices[i].resize(positions.size());
        for (unsigned int j = 0; j < positions.size(); ++j) {
            faces_nodes_indices[i][j] = node_counter++;
            coords.push_back(positions[j]);
        }

    }


    std::vector< double > volInt(exponents[0].size(),0.);
    std::vector< std::vector< double > > bdryInt(faces.size(), std::vector< double > (exponents[0].size(), 0.)); //not needed (just for simplicity for now);
    std::vector< double > normal(3);

    compVol<Coord>( volInt, bdryInt, faces_indices, faces_nodes_indices, coords, _exponents, exponentsAlphaBeta, factorials, normal );

    return volInt;
}

template<class Coord>
std::vector<double> integrate(const std::vector<std::vector<Coord>> & faces)
{
    std::array<std::vector<unsigned int>, 3> exponents = {{
                                                                  std::vector<unsigned int>(84,0), // exponents of x
                                                                  std::vector<unsigned int>(84,0), // exponents of y
                                                                  std::vector<unsigned int>(84,0)  // exponents of z
                                                          }};


    std::size_t counter = 0;
    for( unsigned int degX = 0; degX < 7; ++degX )
    {
        for( unsigned int degY = 0; degY < 7 - degX; ++degY )
        {
            for( unsigned int degZ = 0; degZ < 7 - degX - degY; ++degZ )
            {
                exponents[0][counter] = degX;
                exponents[1][counter] = degY;
                exponents[2][counter] = degZ;

                ++counter;
            }
        }
    }

    std::vector<std::vector<int> > exponentsAlphaBeta(2, std::vector< int > (45));
    counter = 0;
    for( unsigned int degA = 0; degA < 9; ++degA )
    {
        for( unsigned int degB = 0; degB < 9 - degA; ++degB )
        {
            exponentsAlphaBeta[0][counter] = degA;
            exponentsAlphaBeta[1][counter] = degB;
            ++counter;
        }
    }


    return integrate<Coord>(faces, exponents);
}

#endif //CARIBOU_MIRTICHINTEGRATION_H
