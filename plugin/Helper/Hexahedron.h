#ifndef CARIBOU_HELPER_HEXAHEDRON_H
#define CARIBOU_HELPER_HEXAHEDRON_H

#include <array>
#include <vector>
#include <list>
#include <deque>
#include <stack>
#include <algorithm>
#include <cassert>
#include <sstream>
#include <iostream>
#include <numeric>

namespace sofa {

namespace caribou {

namespace helper {

namespace hexahedron {

static const int MarchingCubeEdgeTable[256] =
        {
                0x0  , 0x109, 0x203, 0x30a, 0x406, 0x50f, 0x605, 0x70c,
                0x80c, 0x905, 0xa0f, 0xb06, 0xc0a, 0xd03, 0xe09, 0xf00,
                0x190, 0x99 , 0x393, 0x29a, 0x596, 0x49f, 0x795, 0x69c,
                0x99c, 0x895, 0xb9f, 0xa96, 0xd9a, 0xc93, 0xf99, 0xe90,
                0x230, 0x339, 0x33 , 0x13a, 0x636, 0x73f, 0x435, 0x53c,
                0xa3c, 0xb35, 0x83f, 0x936, 0xe3a, 0xf33, 0xc39, 0xd30,
                0x3a0, 0x2a9, 0x1a3, 0xaa , 0x7a6, 0x6af, 0x5a5, 0x4ac,
                0xbac, 0xaa5, 0x9af, 0x8a6, 0xfaa, 0xea3, 0xda9, 0xca0,
                0x460, 0x569, 0x663, 0x76a, 0x66 , 0x16f, 0x265, 0x36c,
                0xc6c, 0xd65, 0xe6f, 0xf66, 0x86a, 0x963, 0xa69, 0xb60,
                0x5f0, 0x4f9, 0x7f3, 0x6fa, 0x1f6, 0xff , 0x3f5, 0x2fc,
                0xdfc, 0xcf5, 0xfff, 0xef6, 0x9fa, 0x8f3, 0xbf9, 0xaf0,
                0x650, 0x759, 0x453, 0x55a, 0x256, 0x35f, 0x55 , 0x15c,
                0xe5c, 0xf55, 0xc5f, 0xd56, 0xa5a, 0xb53, 0x859, 0x950,
                0x7c0, 0x6c9, 0x5c3, 0x4ca, 0x3c6, 0x2cf, 0x1c5, 0xcc ,
                0xfcc, 0xec5, 0xdcf, 0xcc6, 0xbca, 0xac3, 0x9c9, 0x8c0,
                0x8c0, 0x9c9, 0xac3, 0xbca, 0xcc6, 0xdcf, 0xec5, 0xfcc,
                0xcc , 0x1c5, 0x2cf, 0x3c6, 0x4ca, 0x5c3, 0x6c9, 0x7c0,
                0x950, 0x859, 0xb53, 0xa5a, 0xd56, 0xc5f, 0xf55, 0xe5c,
                0x15c, 0x55 , 0x35f, 0x256, 0x55a, 0x453, 0x759, 0x650,
                0xaf0, 0xbf9, 0x8f3, 0x9fa, 0xef6, 0xfff, 0xcf5, 0xdfc,
                0x2fc, 0x3f5, 0xff , 0x1f6, 0x6fa, 0x7f3, 0x4f9, 0x5f0,
                0xb60, 0xa69, 0x963, 0x86a, 0xf66, 0xe6f, 0xd65, 0xc6c,
                0x36c, 0x265, 0x16f, 0x66 , 0x76a, 0x663, 0x569, 0x460,
                0xca0, 0xda9, 0xea3, 0xfaa, 0x8a6, 0x9af, 0xaa5, 0xbac,
                0x4ac, 0x5a5, 0x6af, 0x7a6, 0xaa , 0x1a3, 0x2a9, 0x3a0,
                0xd30, 0xc39, 0xf33, 0xe3a, 0x936, 0x83f, 0xb35, 0xa3c,
                0x53c, 0x435, 0x73f, 0x636, 0x13a, 0x33 , 0x339, 0x230,
                0xe90, 0xf99, 0xc93, 0xd9a, 0xa96, 0xb9f, 0x895, 0x99c,
                0x69c, 0x795, 0x49f, 0x596, 0x29a, 0x393, 0x99 , 0x190,
                0xf00, 0xe09, 0xd03, 0xc0a, 0xb06, 0xa0f, 0x905, 0x80c,
                0x70c, 0x605, 0x50f, 0x406, 0x30a, 0x203, 0x109, 0x0
        };

static const int MarchingCubeFaceTable[256] =
        {
                0x0 , 0x19, 0x15, 0x1d, 0x25, 0x3d, 0x35, 0x3d, 0x29, 0x39,
                0x3d, 0x3d, 0x2d, 0x3d, 0x3d, 0x3c, 0x1a, 0x1b, 0x1f, 0x1f,
                0x3f, 0x3f, 0x3f, 0x3f, 0x3b, 0x3b, 0x3f, 0x3f, 0x3f, 0x3f,
                0x3f, 0x3e, 0x16, 0x1f, 0x17, 0x1f, 0x37, 0x3f, 0x37, 0x3f,
                0x3f, 0x3f, 0x3f, 0x3f, 0x3f, 0x3f, 0x3f, 0x3e, 0x1e, 0x1f,
                0x1f, 0xf , 0x3f, 0x3f, 0x3f, 0x2f, 0x3f, 0x3f, 0x3f, 0x2f,
                0x3f, 0x3f, 0x3f, 0x2e, 0x26, 0x3f, 0x37, 0x3f, 0x27, 0x3f,
                0x37, 0x3f, 0x2f, 0x3f, 0x3f, 0x3f, 0x2f, 0x3f, 0x3f, 0x3e,
                0x3e, 0x3f, 0x3f, 0x3f, 0x3f, 0x3f, 0x3f, 0x3f, 0x3f, 0x3f,
                0x3f, 0x3f, 0x3f, 0x3f, 0x3f, 0x3e, 0x36, 0x3f, 0x37, 0x3f,
                0x37, 0x3f, 0x33, 0x3b, 0x3f, 0x3f, 0x3f, 0x3f, 0x3f, 0x3f,
                0x3b, 0x3a, 0x3e, 0x3f, 0x3f, 0x2f, 0x3f, 0x3f, 0x3b, 0x2b,
                0x3f, 0x3f, 0x3f, 0x2f, 0x3f, 0x3f, 0x3b, 0x2a, 0x2a, 0x3b,
                0x3f, 0x3f, 0x2f, 0x3f, 0x3f, 0x3f, 0x2b, 0x3b, 0x3f, 0x3f,
                0x2f, 0x3f, 0x3f, 0x3e, 0x3a, 0x3b, 0x3f, 0x3f, 0x3f, 0x3f,
                0x3f, 0x3f, 0x3b, 0x33, 0x3f, 0x37, 0x3f, 0x37, 0x3f, 0x36,
                0x3e, 0x3f, 0x3f, 0x3f, 0x3f, 0x3f, 0x3f, 0x3f, 0x3f, 0x3f,
                0x3f, 0x3f, 0x3f, 0x3f, 0x3f, 0x3e, 0x3e, 0x3f, 0x3f, 0x2f,
                0x3f, 0x3f, 0x3f, 0x2f, 0x3f, 0x37, 0x3f, 0x27, 0x3f, 0x37,
                0x3f, 0x26, 0x2e, 0x3f, 0x3f, 0x3f, 0x2f, 0x3f, 0x3f, 0x3f,
                0x2f, 0x3f, 0x3f, 0x3f, 0xf , 0x1f, 0x1f, 0x1e, 0x3e, 0x3f,
                0x3f, 0x3f, 0x3f, 0x3f, 0x3f, 0x3f, 0x3f, 0x37, 0x3f, 0x37,
                0x1f, 0x17, 0x1f, 0x16, 0x3e, 0x3f, 0x3f, 0x3f, 0x3f, 0x3f,
                0x3b, 0x3b, 0x3f, 0x3f, 0x3f, 0x3f, 0x1f, 0x1f, 0x1b, 0x1a,
                0x3c, 0x3d, 0x3d, 0x2d, 0x3d, 0x3d, 0x39, 0x29, 0x3d, 0x35,
                0x3d, 0x25, 0x1d, 0x15, 0x19, 0x0
        };

static const int MarchingCubeTriTable[256][16] =
        {
                {-1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1},
                {0, 8, 3, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1},
                {0, 1, 9, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1},
                {1, 8, 3, 9, 8, 1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1},
                {1, 2, 10, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1},
                {0, 8, 3, 1, 2, 10, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1},
                {9, 2, 10, 0, 2, 9, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1},
                {2, 8, 3, 2, 10, 8, 10, 9, 8, -1, -1, -1, -1, -1, -1, -1},
                {3, 11, 2, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1},
                {0, 11, 2, 8, 11, 0, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1},
                {1, 9, 0, 2, 3, 11, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1},
                {1, 11, 2, 1, 9, 11, 9, 8, 11, -1, -1, -1, -1, -1, -1, -1},
                {3, 10, 1, 11, 10, 3, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1},
                {0, 10, 1, 0, 8, 10, 8, 11, 10, -1, -1, -1, -1, -1, -1, -1},
                {3, 9, 0, 3, 11, 9, 11, 10, 9, -1, -1, -1, -1, -1, -1, -1},
                {9, 8, 10, 10, 8, 11, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1},
                {4, 7, 8, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1},
                {4, 3, 0, 7, 3, 4, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1},
                {0, 1, 9, 8, 4, 7, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1},
                {4, 1, 9, 4, 7, 1, 7, 3, 1, -1, -1, -1, -1, -1, -1, -1},
                {1, 2, 10, 8, 4, 7, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1},
                {3, 4, 7, 3, 0, 4, 1, 2, 10, -1, -1, -1, -1, -1, -1, -1},
                {9, 2, 10, 9, 0, 2, 8, 4, 7, -1, -1, -1, -1, -1, -1, -1},
                {2, 10, 9, 2, 9, 7, 2, 7, 3, 7, 9, 4, -1, -1, -1, -1},
                {8, 4, 7, 3, 11, 2, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1},
                {11, 4, 7, 11, 2, 4, 2, 0, 4, -1, -1, -1, -1, -1, -1, -1},
                {9, 0, 1, 8, 4, 7, 2, 3, 11, -1, -1, -1, -1, -1, -1, -1},
                {4, 7, 11, 9, 4, 11, 9, 11, 2, 9, 2, 1, -1, -1, -1, -1},
                {3, 10, 1, 3, 11, 10, 7, 8, 4, -1, -1, -1, -1, -1, -1, -1},
                {1, 11, 10, 1, 4, 11, 1, 0, 4, 7, 11, 4, -1, -1, -1, -1},
                {4, 7, 8, 9, 0, 11, 9, 11, 10, 11, 0, 3, -1, -1, -1, -1},
                {4, 7, 11, 4, 11, 9, 9, 11, 10, -1, -1, -1, -1, -1, -1, -1},
                {9, 5, 4, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1},
                {9, 5, 4, 0, 8, 3, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1},
                {0, 5, 4, 1, 5, 0, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1},
                {8, 5, 4, 8, 3, 5, 3, 1, 5, -1, -1, -1, -1, -1, -1, -1},
                {1, 2, 10, 9, 5, 4, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1},
                {3, 0, 8, 1, 2, 10, 4, 9, 5, -1, -1, -1, -1, -1, -1, -1},
                {5, 2, 10, 5, 4, 2, 4, 0, 2, -1, -1, -1, -1, -1, -1, -1},
                {2, 10, 5, 3, 2, 5, 3, 5, 4, 3, 4, 8, -1, -1, -1, -1},
                {9, 5, 4, 2, 3, 11, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1},
                {0, 11, 2, 0, 8, 11, 4, 9, 5, -1, -1, -1, -1, -1, -1, -1},
                {0, 5, 4, 0, 1, 5, 2, 3, 11, -1, -1, -1, -1, -1, -1, -1},
                {2, 1, 5, 2, 5, 8, 2, 8, 11, 4, 8, 5, -1, -1, -1, -1},
                {10, 3, 11, 10, 1, 3, 9, 5, 4, -1, -1, -1, -1, -1, -1, -1},
                {4, 9, 5, 0, 8, 1, 8, 10, 1, 8, 11, 10, -1, -1, -1, -1},
                {5, 4, 0, 5, 0, 11, 5, 11, 10, 11, 0, 3, -1, -1, -1, -1},
                {5, 4, 8, 5, 8, 10, 10, 8, 11, -1, -1, -1, -1, -1, -1, -1},
                {9, 7, 8, 5, 7, 9, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1},
                {9, 3, 0, 9, 5, 3, 5, 7, 3, -1, -1, -1, -1, -1, -1, -1},
                {0, 7, 8, 0, 1, 7, 1, 5, 7, -1, -1, -1, -1, -1, -1, -1},
                {1, 5, 3, 3, 5, 7, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1},
                {9, 7, 8, 9, 5, 7, 10, 1, 2, -1, -1, -1, -1, -1, -1, -1},
                {10, 1, 2, 9, 5, 0, 5, 3, 0, 5, 7, 3, -1, -1, -1, -1},
                {8, 0, 2, 8, 2, 5, 8, 5, 7, 10, 5, 2, -1, -1, -1, -1},
                {2, 10, 5, 2, 5, 3, 3, 5, 7, -1, -1, -1, -1, -1, -1, -1},
                {7, 9, 5, 7, 8, 9, 3, 11, 2, -1, -1, -1, -1, -1, -1, -1},
                {9, 5, 7, 9, 7, 2, 9, 2, 0, 2, 7, 11, -1, -1, -1, -1},
                {2, 3, 11, 0, 1, 8, 1, 7, 8, 1, 5, 7, -1, -1, -1, -1},
                {11, 2, 1, 11, 1, 7, 7, 1, 5, -1, -1, -1, -1, -1, -1, -1},
                {9, 5, 8, 8, 5, 7, 10, 1, 3, 10, 3, 11, -1, -1, -1, -1},
                {5, 7, 0, 5, 0, 9, 7, 11, 0, 1, 0, 10, 11, 10, 0, -1},
                {11, 10, 0, 11, 0, 3, 10, 5, 0, 8, 0, 7, 5, 7, 0, -1},
                {11, 10, 5, 7, 11, 5, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1},
                {10, 6, 5, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1},
                {0, 8, 3, 5, 10, 6, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1},
                {9, 0, 1, 5, 10, 6, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1},
                {1, 8, 3, 1, 9, 8, 5, 10, 6, -1, -1, -1, -1, -1, -1, -1},
                {1, 6, 5, 2, 6, 1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1},
                {1, 6, 5, 1, 2, 6, 3, 0, 8, -1, -1, -1, -1, -1, -1, -1},
                {9, 6, 5, 9, 0, 6, 0, 2, 6, -1, -1, -1, -1, -1, -1, -1},
                {5, 9, 8, 5, 8, 2, 5, 2, 6, 3, 2, 8, -1, -1, -1, -1},
                {2, 3, 11, 10, 6, 5, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1},
                {11, 0, 8, 11, 2, 0, 10, 6, 5, -1, -1, -1, -1, -1, -1, -1},
                {0, 1, 9, 2, 3, 11, 5, 10, 6, -1, -1, -1, -1, -1, -1, -1},
                {5, 10, 6, 1, 9, 2, 9, 11, 2, 9, 8, 11, -1, -1, -1, -1},
                {6, 3, 11, 6, 5, 3, 5, 1, 3, -1, -1, -1, -1, -1, -1, -1},
                {0, 8, 11, 0, 11, 5, 0, 5, 1, 5, 11, 6, -1, -1, -1, -1},
                {3, 11, 6, 0, 3, 6, 0, 6, 5, 0, 5, 9, -1, -1, -1, -1},
                {6, 5, 9, 6, 9, 11, 11, 9, 8, -1, -1, -1, -1, -1, -1, -1},
                {5, 10, 6, 4, 7, 8, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1},
                {4, 3, 0, 4, 7, 3, 6, 5, 10, -1, -1, -1, -1, -1, -1, -1},
                {1, 9, 0, 5, 10, 6, 8, 4, 7, -1, -1, -1, -1, -1, -1, -1},
                {10, 6, 5, 1, 9, 7, 1, 7, 3, 7, 9, 4, -1, -1, -1, -1},
                {6, 1, 2, 6, 5, 1, 4, 7, 8, -1, -1, -1, -1, -1, -1, -1},
                {1, 2, 5, 5, 2, 6, 3, 0, 4, 3, 4, 7, -1, -1, -1, -1},
                {8, 4, 7, 9, 0, 5, 0, 6, 5, 0, 2, 6, -1, -1, -1, -1},
                {7, 3, 9, 7, 9, 4, 3, 2, 9, 5, 9, 6, 2, 6, 9, -1},
                {3, 11, 2, 7, 8, 4, 10, 6, 5, -1, -1, -1, -1, -1, -1, -1},
                {5, 10, 6, 4, 7, 2, 4, 2, 0, 2, 7, 11, -1, -1, -1, -1},
                {0, 1, 9, 4, 7, 8, 2, 3, 11, 5, 10, 6, -1, -1, -1, -1},
                {9, 2, 1, 9, 11, 2, 9, 4, 11, 7, 11, 4, 5, 10, 6, -1},
                {8, 4, 7, 3, 11, 5, 3, 5, 1, 5, 11, 6, -1, -1, -1, -1},
                {5, 1, 11, 5, 11, 6, 1, 0, 11, 7, 11, 4, 0, 4, 11, -1},
                {0, 5, 9, 0, 6, 5, 0, 3, 6, 11, 6, 3, 8, 4, 7, -1},
                {6, 5, 9, 6, 9, 11, 4, 7, 9, 7, 11, 9, -1, -1, -1, -1},
                {10, 4, 9, 6, 4, 10, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1},
                {4, 10, 6, 4, 9, 10, 0, 8, 3, -1, -1, -1, -1, -1, -1, -1},
                {10, 0, 1, 10, 6, 0, 6, 4, 0, -1, -1, -1, -1, -1, -1, -1},
                {8, 3, 1, 8, 1, 6, 8, 6, 4, 6, 1, 10, -1, -1, -1, -1},
                {1, 4, 9, 1, 2, 4, 2, 6, 4, -1, -1, -1, -1, -1, -1, -1},
                {3, 0, 8, 1, 2, 9, 2, 4, 9, 2, 6, 4, -1, -1, -1, -1},
                {0, 2, 4, 4, 2, 6, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1},
                {8, 3, 2, 8, 2, 4, 4, 2, 6, -1, -1, -1, -1, -1, -1, -1},
                {10, 4, 9, 10, 6, 4, 11, 2, 3, -1, -1, -1, -1, -1, -1, -1},
                {0, 8, 2, 2, 8, 11, 4, 9, 10, 4, 10, 6, -1, -1, -1, -1},
                {3, 11, 2, 0, 1, 6, 0, 6, 4, 6, 1, 10, -1, -1, -1, -1},
                {6, 4, 1, 6, 1, 10, 4, 8, 1, 2, 1, 11, 8, 11, 1, -1},
                {9, 6, 4, 9, 3, 6, 9, 1, 3, 11, 6, 3, -1, -1, -1, -1},
                {8, 11, 1, 8, 1, 0, 11, 6, 1, 9, 1, 4, 6, 4, 1, -1},
                {3, 11, 6, 3, 6, 0, 0, 6, 4, -1, -1, -1, -1, -1, -1, -1},
                {6, 4, 8, 11, 6, 8, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1},
                {7, 10, 6, 7, 8, 10, 8, 9, 10, -1, -1, -1, -1, -1, -1, -1},
                {0, 7, 3, 0, 10, 7, 0, 9, 10, 6, 7, 10, -1, -1, -1, -1},
                {10, 6, 7, 1, 10, 7, 1, 7, 8, 1, 8, 0, -1, -1, -1, -1},
                {10, 6, 7, 10, 7, 1, 1, 7, 3, -1, -1, -1, -1, -1, -1, -1},
                {1, 2, 6, 1, 6, 8, 1, 8, 9, 8, 6, 7, -1, -1, -1, -1},
                {2, 6, 9, 2, 9, 1, 6, 7, 9, 0, 9, 3, 7, 3, 9, -1},
                {7, 8, 0, 7, 0, 6, 6, 0, 2, -1, -1, -1, -1, -1, -1, -1},
                {7, 3, 2, 6, 7, 2, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1},
                {2, 3, 11, 10, 6, 8, 10, 8, 9, 8, 6, 7, -1, -1, -1, -1},
                {2, 0, 7, 2, 7, 11, 0, 9, 7, 6, 7, 10, 9, 10, 7, -1},
                {1, 8, 0, 1, 7, 8, 1, 10, 7, 6, 7, 10, 2, 3, 11, -1},
                {11, 2, 1, 11, 1, 7, 10, 6, 1, 6, 7, 1, -1, -1, -1, -1},
                {8, 9, 6, 8, 6, 7, 9, 1, 6, 11, 6, 3, 1, 3, 6, -1},
                {0, 9, 1, 11, 6, 7, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1},
                {7, 8, 0, 7, 0, 6, 3, 11, 0, 11, 6, 0, -1, -1, -1, -1},
                {7, 11, 6, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1},
                {7, 6, 11, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1},
                {3, 0, 8, 11, 7, 6, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1},
                {0, 1, 9, 11, 7, 6, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1},
                {8, 1, 9, 8, 3, 1, 11, 7, 6, -1, -1, -1, -1, -1, -1, -1},
                {10, 1, 2, 6, 11, 7, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1},
                {1, 2, 10, 3, 0, 8, 6, 11, 7, -1, -1, -1, -1, -1, -1, -1},
                {2, 9, 0, 2, 10, 9, 6, 11, 7, -1, -1, -1, -1, -1, -1, -1},
                {6, 11, 7, 2, 10, 3, 10, 8, 3, 10, 9, 8, -1, -1, -1, -1},
                {7, 2, 3, 6, 2, 7, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1},
                {7, 0, 8, 7, 6, 0, 6, 2, 0, -1, -1, -1, -1, -1, -1, -1},
                {2, 7, 6, 2, 3, 7, 0, 1, 9, -1, -1, -1, -1, -1, -1, -1},
                {1, 6, 2, 1, 8, 6, 1, 9, 8, 8, 7, 6, -1, -1, -1, -1},
                {10, 7, 6, 10, 1, 7, 1, 3, 7, -1, -1, -1, -1, -1, -1, -1},
                {10, 7, 6, 1, 7, 10, 1, 8, 7, 1, 0, 8, -1, -1, -1, -1},
                {0, 3, 7, 0, 7, 10, 0, 10, 9, 6, 10, 7, -1, -1, -1, -1},
                {7, 6, 10, 7, 10, 8, 8, 10, 9, -1, -1, -1, -1, -1, -1, -1},
                {6, 8, 4, 11, 8, 6, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1},
                {3, 6, 11, 3, 0, 6, 0, 4, 6, -1, -1, -1, -1, -1, -1, -1},
                {8, 6, 11, 8, 4, 6, 9, 0, 1, -1, -1, -1, -1, -1, -1, -1},
                {9, 4, 6, 9, 6, 3, 9, 3, 1, 11, 3, 6, -1, -1, -1, -1},
                {6, 8, 4, 6, 11, 8, 2, 10, 1, -1, -1, -1, -1, -1, -1, -1},
                {1, 2, 10, 3, 0, 11, 0, 6, 11, 0, 4, 6, -1, -1, -1, -1},
                {4, 11, 8, 4, 6, 11, 0, 2, 9, 2, 10, 9, -1, -1, -1, -1},
                {10, 9, 3, 10, 3, 2, 9, 4, 3, 11, 3, 6, 4, 6, 3, -1},
                {8, 2, 3, 8, 4, 2, 4, 6, 2, -1, -1, -1, -1, -1, -1, -1},
                {0, 4, 2, 4, 6, 2, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1},
                {1, 9, 0, 2, 3, 4, 2, 4, 6, 4, 3, 8, -1, -1, -1, -1},
                {1, 9, 4, 1, 4, 2, 2, 4, 6, -1, -1, -1, -1, -1, -1, -1},
                {8, 1, 3, 8, 6, 1, 8, 4, 6, 6, 10, 1, -1, -1, -1, -1},
                {10, 1, 0, 10, 0, 6, 6, 0, 4, -1, -1, -1, -1, -1, -1, -1},
                {4, 6, 3, 4, 3, 8, 6, 10, 3, 0, 3, 9, 10, 9, 3, -1},
                {10, 9, 4, 6, 10, 4, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1},
                {4, 9, 5, 7, 6, 11, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1},
                {0, 8, 3, 4, 9, 5, 11, 7, 6, -1, -1, -1, -1, -1, -1, -1},
                {5, 0, 1, 5, 4, 0, 7, 6, 11, -1, -1, -1, -1, -1, -1, -1},
                {11, 7, 6, 8, 3, 4, 3, 5, 4, 3, 1, 5, -1, -1, -1, -1},
                {9, 5, 4, 10, 1, 2, 7, 6, 11, -1, -1, -1, -1, -1, -1, -1},
                {6, 11, 7, 1, 2, 10, 0, 8, 3, 4, 9, 5, -1, -1, -1, -1},
                {7, 6, 11, 5, 4, 10, 4, 2, 10, 4, 0, 2, -1, -1, -1, -1},
                {3, 4, 8, 3, 5, 4, 3, 2, 5, 10, 5, 2, 11, 7, 6, -1},
                {7, 2, 3, 7, 6, 2, 5, 4, 9, -1, -1, -1, -1, -1, -1, -1},
                {9, 5, 4, 0, 8, 6, 0, 6, 2, 6, 8, 7, -1, -1, -1, -1},
                {3, 6, 2, 3, 7, 6, 1, 5, 0, 5, 4, 0, -1, -1, -1, -1},
                {6, 2, 8, 6, 8, 7, 2, 1, 8, 4, 8, 5, 1, 5, 8, -1},
                {9, 5, 4, 10, 1, 6, 1, 7, 6, 1, 3, 7, -1, -1, -1, -1},
                {1, 6, 10, 1, 7, 6, 1, 0, 7, 8, 7, 0, 9, 5, 4, -1},
                {4, 0, 10, 4, 10, 5, 0, 3, 10, 6, 10, 7, 3, 7, 10, -1},
                {7, 6, 10, 7, 10, 8, 5, 4, 10, 4, 8, 10, -1, -1, -1, -1},
                {6, 9, 5, 6, 11, 9, 11, 8, 9, -1, -1, -1, -1, -1, -1, -1},
                {3, 6, 11, 0, 6, 3, 0, 5, 6, 0, 9, 5, -1, -1, -1, -1},
                {0, 11, 8, 0, 5, 11, 0, 1, 5, 5, 6, 11, -1, -1, -1, -1},
                {6, 11, 3, 6, 3, 5, 5, 3, 1, -1, -1, -1, -1, -1, -1, -1},
                {1, 2, 10, 9, 5, 11, 9, 11, 8, 11, 5, 6, -1, -1, -1, -1},
                {0, 11, 3, 0, 6, 11, 0, 9, 6, 5, 6, 9, 1, 2, 10, -1},
                {11, 8, 5, 11, 5, 6, 8, 0, 5, 10, 5, 2, 0, 2, 5, -1},
                {6, 11, 3, 6, 3, 5, 2, 10, 3, 10, 5, 3, -1, -1, -1, -1},
                {5, 8, 9, 5, 2, 8, 5, 6, 2, 3, 8, 2, -1, -1, -1, -1},
                {9, 5, 6, 9, 6, 0, 0, 6, 2, -1, -1, -1, -1, -1, -1, -1},
                {1, 5, 8, 1, 8, 0, 5, 6, 8, 3, 8, 2, 6, 2, 8, -1},
                {1, 5, 6, 2, 1, 6, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1},
                {1, 3, 6, 1, 6, 10, 3, 8, 6, 5, 6, 9, 8, 9, 6, -1},
                {10, 1, 0, 10, 0, 6, 9, 5, 0, 5, 6, 0, -1, -1, -1, -1},
                {0, 3, 8, 5, 6, 10, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1},
                {10, 5, 6, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1},
                {11, 5, 10, 7, 5, 11, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1},
                {11, 5, 10, 11, 7, 5, 8, 3, 0, -1, -1, -1, -1, -1, -1, -1},
                {5, 11, 7, 5, 10, 11, 1, 9, 0, -1, -1, -1, -1, -1, -1, -1},
                {10, 7, 5, 10, 11, 7, 9, 8, 1, 8, 3, 1, -1, -1, -1, -1},
                {11, 1, 2, 11, 7, 1, 7, 5, 1, -1, -1, -1, -1, -1, -1, -1},
                {0, 8, 3, 1, 2, 7, 1, 7, 5, 7, 2, 11, -1, -1, -1, -1},
                {9, 7, 5, 9, 2, 7, 9, 0, 2, 2, 11, 7, -1, -1, -1, -1},
                {7, 5, 2, 7, 2, 11, 5, 9, 2, 3, 2, 8, 9, 8, 2, -1},
                {2, 5, 10, 2, 3, 5, 3, 7, 5, -1, -1, -1, -1, -1, -1, -1},
                {8, 2, 0, 8, 5, 2, 8, 7, 5, 10, 2, 5, -1, -1, -1, -1},
                {9, 0, 1, 5, 10, 3, 5, 3, 7, 3, 10, 2, -1, -1, -1, -1},
                {9, 8, 2, 9, 2, 1, 8, 7, 2, 10, 2, 5, 7, 5, 2, -1},
                {1, 3, 5, 3, 7, 5, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1},
                {0, 8, 7, 0, 7, 1, 1, 7, 5, -1, -1, -1, -1, -1, -1, -1},
                {9, 0, 3, 9, 3, 5, 5, 3, 7, -1, -1, -1, -1, -1, -1, -1},
                {9, 8, 7, 5, 9, 7, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1},
                {5, 8, 4, 5, 10, 8, 10, 11, 8, -1, -1, -1, -1, -1, -1, -1},
                {5, 0, 4, 5, 11, 0, 5, 10, 11, 11, 3, 0, -1, -1, -1, -1},
                {0, 1, 9, 8, 4, 10, 8, 10, 11, 10, 4, 5, -1, -1, -1, -1},
                {10, 11, 4, 10, 4, 5, 11, 3, 4, 9, 4, 1, 3, 1, 4, -1},
                {2, 5, 1, 2, 8, 5, 2, 11, 8, 4, 5, 8, -1, -1, -1, -1},
                {0, 4, 11, 0, 11, 3, 4, 5, 11, 2, 11, 1, 5, 1, 11, -1},
                {0, 2, 5, 0, 5, 9, 2, 11, 5, 4, 5, 8, 11, 8, 5, -1},
                {9, 4, 5, 2, 11, 3, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1},
                {2, 5, 10, 3, 5, 2, 3, 4, 5, 3, 8, 4, -1, -1, -1, -1},
                {5, 10, 2, 5, 2, 4, 4, 2, 0, -1, -1, -1, -1, -1, -1, -1},
                {3, 10, 2, 3, 5, 10, 3, 8, 5, 4, 5, 8, 0, 1, 9, -1},
                {5, 10, 2, 5, 2, 4, 1, 9, 2, 9, 4, 2, -1, -1, -1, -1},
                {8, 4, 5, 8, 5, 3, 3, 5, 1, -1, -1, -1, -1, -1, -1, -1},
                {0, 4, 5, 1, 0, 5, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1},
                {8, 4, 5, 8, 5, 3, 9, 0, 5, 0, 3, 5, -1, -1, -1, -1},
                {9, 4, 5, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1},
                {4, 11, 7, 4, 9, 11, 9, 10, 11, -1, -1, -1, -1, -1, -1, -1},
                {0, 8, 3, 4, 9, 7, 9, 11, 7, 9, 10, 11, -1, -1, -1, -1},
                {1, 10, 11, 1, 11, 4, 1, 4, 0, 7, 4, 11, -1, -1, -1, -1},
                {3, 1, 4, 3, 4, 8, 1, 10, 4, 7, 4, 11, 10, 11, 4, -1},
                {4, 11, 7, 9, 11, 4, 9, 2, 11, 9, 1, 2, -1, -1, -1, -1},
                {9, 7, 4, 9, 11, 7, 9, 1, 11, 2, 11, 1, 0, 8, 3, -1},
                {11, 7, 4, 11, 4, 2, 2, 4, 0, -1, -1, -1, -1, -1, -1, -1},
                {11, 7, 4, 11, 4, 2, 8, 3, 4, 3, 2, 4, -1, -1, -1, -1},
                {2, 9, 10, 2, 7, 9, 2, 3, 7, 7, 4, 9, -1, -1, -1, -1},
                {9, 10, 7, 9, 7, 4, 10, 2, 7, 8, 7, 0, 2, 0, 7, -1},
                {3, 7, 10, 3, 10, 2, 7, 4, 10, 1, 10, 0, 4, 0, 10, -1},
                {1, 10, 2, 8, 7, 4, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1},
                {4, 9, 1, 4, 1, 7, 7, 1, 3, -1, -1, -1, -1, -1, -1, -1},
                {4, 9, 1, 4, 1, 7, 0, 8, 1, 8, 7, 1, -1, -1, -1, -1},
                {4, 0, 3, 7, 4, 3, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1},
                {4, 8, 7, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1},
                {9, 10, 8, 10, 11, 8, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1},
                {3, 0, 9, 3, 9, 11, 11, 9, 10, -1, -1, -1, -1, -1, -1, -1},
                {0, 1, 10, 0, 10, 8, 8, 10, 11, -1, -1, -1, -1, -1, -1, -1},
                {3, 1, 10, 11, 3, 10, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1},
                {1, 2, 11, 1, 11, 9, 9, 11, 8, -1, -1, -1, -1, -1, -1, -1},
                {3, 0, 9, 3, 9, 11, 1, 2, 9, 2, 11, 9, -1, -1, -1, -1},
                {0, 2, 11, 8, 0, 11, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1},
                {3, 2, 11, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1},
                {2, 3, 8, 2, 8, 10, 10, 8, 9, -1, -1, -1, -1, -1, -1, -1},
                {9, 10, 2, 0, 9, 2, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1},
                {2, 3, 8, 2, 8, 10, 0, 1, 8, 1, 10, 8, -1, -1, -1, -1},
                {1, 10, 2, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1},
                {1, 3, 8, 9, 1, 8, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1},
                {0, 9, 1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1},
                {0, 3, 8, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1},
                {-1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1}
        };

using EdgeId = char;
using NodeId = char;
using AlphaValue = float;

/* Hexahedron node numbering
 *    7-------6
 *   /|      /|
 *  / |     / |
 * 3--|----2  |
 * |  4----|--5
 * | /     | /
 * 0-------1
 */

/* Cut hexahedron node numbering
 *        7-------------6
 *       /|            /|
 *      / |           / |
 *     /  |          /  |
 *    6---|--5------4   |
 *    |   |         |   |
 *    |   |         |   |
 *    7   |         3   |
 *    |   4---------|---5
 *    |  /          |  /
 *    | /           | /
 *    0------1------2
 */

/**
 * List of edges by node ID. For example, edges[5] = {5, 6} is represented by node 5 and node 6.
 */
static const std::array<std::array<NodeId, 2>, 12> edges {{
    {0, 1}, // 0
    {1, 2}, // 1
    {2, 3}, // 2
    {3, 0}, // 3
    {4, 5}, // 4
    {5, 6}, // 5
    {6, 7}, // 6
    {7, 4}, // 7
    {0, 4}, // 8
    {1, 5}, // 9
    {2, 6}, // 10
    {3, 7}, // 11
}};

/**
 * List of edge ID for each face. For example, edges_of_face[2] = {6, 5, 4, 7} is represented by
 * edge 6 (nodes 6 and 7), edge 5 (nodes 5 and 6), edge 4 (nodes 4 and 5), edge 7 (nodes 7 and 4) forming the face #2.
 */
static const std::array<std::array<EdgeId, 4>, 6> edges_of_face = {{
    {0, 1, 2, 3},   // 0 ; Nodes = {0 1 2 3}
    {11, 7, 8, 3},  // 1 ; Nodes = {3 7 4 0}
    {6, 5, 4, 7},   // 2 ; Nodes = {7 6 5 4}
    {5, 10, 1, 9},  // 3 ; Nodes = {5 6 2 1}
    {9, 0, 8, 4},   // 4 ; Nodes = {5 1 0 4}
    {10, 6, 11, 2}, // 5 ; Nodes = {2 6 7 3}
}};

/**
 * List of node ID for each face. For example, nodes_of_face[2] = {7, 6, 5, 4} is represented by
 * node 7, node 6, node 5 and node 4 forming the face #2.
 */
static const std::array<std::array<NodeId, 4>, 6> nodes_of_face = {{
    {0, 1, 2, 3}, // 0 ; Edges = {0, 1, 2, 3,}
    {3, 7, 4, 0}, // 1 ; Edges = {11, 7, 8, 3}
    {7, 6, 5, 4}, // 2 ; Edges = {6, 5, 4, 7}
    {5, 6, 2, 1}, // 3 ; Edges = {5, 10, 1, 9}
    {5, 1, 0, 4}, // 4 ; Edges = {9, 0, 8, 4}
    {2, 6, 7, 3}, // 5 ; Edges = {10, 6, 11, 2}
}};

/**
 * List of edges adjacent to a node. For example, edges_adjacent_to_node[2] = {2, 1, 10} is represented by the edges
 * 2, 1 and 10 that are adjacent to the node 2. Here, the order of the edges listed for a node allows for a normal
 * vector pointing outward of the cube interior.
 */
std::array<std::array<EdgeId,3>,8> edges_adjacent_to_node = {{
    {0, 3, 8},  // 0
    {1, 0, 9},  // 1
    {2, 1, 10}, // 2
    {3, 2, 11}, // 3
    {8, 7, 4},  // 4
    {4, 5, 9},  // 5
    {6, 10, 5}, // 6
    {6, 7, 11}, // 7
}};



/// An alpha position represent a relative position on a given edge.
/// The alpha value should be between 0 and 1, 0 being the first node and 1 being the second node of the edge.
using AlphaPosition = std::pair<EdgeId, AlphaValue>;

/**
 * Find the edge index from two node indices. The order of the nodes doesn't impact the search.
 * @param node_id_1 The first node.
 * @param node_id_2 The second node.
 * @return The index of the edge connecting the two nodes, or -1 if no such edge exists.
 */
char find_edge_id_from_two_nodes(const unsigned char node_id_1, const unsigned char node_id_2) {
    std::array<EdgeId,3> edge_indices = edges_adjacent_to_node[node_id_1];
    for (const EdgeId edge_id : edge_indices) {
        const std::array<NodeId, 2> & edge = edges[edge_id];
        if (edge[0] == node_id_2 || edge[1] == node_id_2)
            return edge_id;
    }

    return -1;
}

/**
 * Compute the alpha position of a given node and edge. The resulting alpha value will be 0 if the node is the first
 * node of the edge, or 1 if it is the second node of the edge.
 * @param node_id The node id for which we want the alpha value
 * @param edge_id The edge id for which we want the alpha value
 * @return The alpha position <edge_id, alpha-value>
 */
AlphaPosition alpha_position(const NodeId node_id, const EdgeId edge_id) {
    const std::array<NodeId, 2> & edge = edges[edge_id];
    if (edge[0] == node_id)
        return std::make_pair(edge_id, 0.);
    else
        return std::make_pair(edge_id, 1.);
}

#ifndef	NDEBUG
#define assert_trace(expr, message, trace) if (!(expr)) std::cout<<trace<<std::endl; assert(expr && message);
#else
#define assert_trace
#endif

/**
 * Triangulate an hexahedron cut by an iso-surface. This function will create a mesh of triangles that cover the iso-surface inside an hexahedron (the cut surface).
 * @param isInside [IN] Array of booleans that indicate whether or not one of the hexahedron's node is inside or outside of the iso-surface (cut-surface).
 * @param edgeIntersections [IN] Array of intersection coefficients that indicate where the iso-surface cut an edge.
 *                               Each coefficient should be between 0 and 1, 0 being the position of the first node of
 *                               the edge, 1 being the second node, and anywhere between is the exact location of the cut.
 * @return List of triangles where a triangle is represented by three edge intersections (an instersection is a pair of <edge_id, alpha-position> (see AlphaPosition)
 */
std::vector<std::array<AlphaPosition, 3>> triangulate(const std::array<bool, 8> &isInside, const std::array<AlphaValue, 12> &edgeIntersections)
{
    using AlphaTriangle = std::array<AlphaPosition, 3>;
    std::vector<AlphaTriangle> triangles;

    // Determine the index into the edge table which tells us which vertices are inside of the surface
    unsigned short flags = 0;
    for (unsigned char i = 0; i < 8; ++i)
        if (isInside[i])
            flags |= (1 << i);

    /* Cube is entirely in/out of the surface */
    if (MarchingCubeEdgeTable[flags] == 0)
        return triangles;

    // STEP 1: Meshing the triangles that cut the hexa
    for (unsigned short i = 0; MarchingCubeTriTable[flags][i] != -1; i += 3) {
        AlphaTriangle triangle;
        for (unsigned char j = 0; j < 3; ++j) {
            const EdgeId &edge_id = (const EdgeId) MarchingCubeTriTable[flags][i + j];
            AlphaValue alpha = edgeIntersections[edge_id];
            triangle[j] = std::make_pair(edge_id, alpha);
        }

        // Add the triangle to the output
        triangles.push_back(triangle);
    }

    // STEP 2: Meshing the hexahedron's faces that reside under the cut surface

    // Populate each edges with its connected triangles
    std::array<std::vector<AlphaTriangle>, 12> edge_triangles;
    for (const AlphaTriangle & t : triangles) {
        for (const AlphaPosition & p : t) {
            const EdgeId  edgeId = p.first;
            edge_triangles[edgeId].push_back(t);
        }
    }

    // Each edge has a maximum of 1 intersection represented by an
    // alpha position (between 0 and 1) : intersection = edge[0] + alpha * (edge[1] - edge[0]).
    // If no intersection is on an edge, the alpha value is -1.
    std::array<AlphaValue, 12> intersection_of_edge = {-1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1};
    for (const AlphaTriangle & triangle : triangles) {
        for (const AlphaPosition & intersection : triangle) {
            const EdgeId & edge_id = intersection.first;
            const AlphaValue & alpha = intersection.second;
            intersection_of_edge[edge_id] = alpha;
        }
    }

#ifndef	NDEBUG
    std::stringstream cube_debug;
    cube_debug << "Cube :";
    for (NodeId i=0;i<8;++i)
        if(isInside[i])
            cube_debug << " 1";
        else cube_debug << " 0";
    cube_debug << '\n';
    cube_debug << "  Triangles:" << '\n';
    for (const AlphaTriangle & triangle : triangles) {
        cube_debug << "    ";
        for (const AlphaPosition & intersection : triangle) {
            const EdgeId & edge_id = intersection.first;
            const AlphaValue & alpha = intersection.second;
            const auto & edge = edges[edge_id];
            cube_debug << "([" << std::to_string(edge[0]) << ", " << std::to_string(edge[1]) << "], " << std::to_string(alpha) << "), ";
        }
        cube_debug << '\n';
    }
#else
    std::stringstream cube_debug;
#endif

    // We loop through each of the hexa's faces. A face can contain one or more subfaces, or none if the face lies
    // completely outside of the cutting surface
    for (unsigned char face_id = 0; face_id < 6; ++face_id) {

        const auto & face_nodes = nodes_of_face[face_id];
        const auto & face_edges = edges_of_face[face_id];

#ifndef	NDEBUG
        std::stringstream face_debug;
        face_debug << "  Face #" << std::to_string(face_id) << " : (";
        for (const NodeId & nid : face_nodes)
            face_debug << std::to_string(nid)<<" ";
        face_debug << ")\n";
#else
        std::stringstream face_debug;
#endif

        // Create the sub-faces of this face that are under the cut surface. A sub-face is represented by a loop of
        // alpha-positions (edge ID + alpha (between 0 and 1)).
        std::vector<std::list<AlphaPosition>> subfaces;

        // Each of the face's node can be the starting point of a subface if it hasn't been already added to a subface yet.
        std::array<bool, 4> node_visited = {false, false, false, false};

        // Each of the face's node can be connected to an intersection found with the marching cube
        std::array<const AlphaPosition *, 4> connected_intersection = {nullptr, nullptr, nullptr, nullptr};

        // Populate the connected intersection to each face's nodes
        for (const AlphaTriangle & triangle : triangles) {
            const AlphaPosition * node = nullptr;
            const AlphaPosition * intersection = nullptr;

            for (const AlphaPosition & position : triangle) {
                // Check if the triangle's vertex is on the current face
                auto edge = std::find(std::begin(face_edges), std::end(face_edges), position.first);
                if (edge != std::end(face_edges)) {
                    // The vertex is on the current face
                    const AlphaValue & alpha = position.second;
                    if (alpha > 0 and alpha < 1) {
                        // The vertex is an intersect edge
                        intersection = &position;
                    } else {
                        // The vertex is one of the face's node
                        node = &position;
                    }
                };
            }

            // If the triangle was not connected to this face
            if (!node || !intersection)
                continue;

            assert_trace(node != nullptr && intersection != nullptr, "Each intersected edges should be connected to a node on the same face.", (cube_debug.str() + face_debug.str()));

            // Find the node id
            char node_id = edges[node->first][(unsigned char) node->second];
            for (unsigned char i = 0; i < 4; ++i) {
                if (face_nodes[i] == node_id) {
                    assert_trace(connected_intersection[i] == nullptr, "Only one intersected edge must be connected to a node for a given face", (cube_debug.str() + face_debug.str()));
                    connected_intersection[i] = intersection;
                    break;
                }
            }

        }

        for (unsigned char node_idx = 0; node_idx < 4; ++node_idx) {
            char node_id = face_nodes[node_idx];
            char next_node_id = face_nodes[(node_idx+1) % 4];

            if (not isInside[node_id])
                continue;

            if (node_visited[node_idx])
                continue;

            // A subface is a set (loop) of alpha-positions on a hexa's face.
            // (they can be either one of the hexahedron's node, or an intersection between two nodes)
            // Example:
            // A------------D
            // |            |
            // |            |
            // E            F
            // | \        / |
            // |  \      /  |
            // B---G----H---C
            //
            // Here (A, D, F, H, G, E), (E, B, G) and (F, H, C) are three subfaces of the quad face (A, B, C, D)

            std::list<AlphaPosition> subface;

            EdgeId currentEdgeId = find_edge_id_from_two_nodes(node_id, next_node_id);
            assert_trace (currentEdgeId>=0, "An edge must exists between two adjacent nodes.", (cube_debug.str() + face_debug.str()));

            subface.push_back(alpha_position(node_id, currentEdgeId));

            bool lastPositionIsAnIntersection = false;

            bool subfaceIsClosed = false; // Will be true once the contour of the subface is closed

            while(not subfaceIsClosed) {
                const AlphaPosition & currentPosition = subface.back();
                const EdgeId & edgeId = currentPosition.first;
                const AlphaValue & alpha = currentPosition.second;

#ifndef	NDEBUG
                const auto & edge = edges[edgeId];
                if (alpha > 0 && alpha < 1)
                    face_debug << "    Visiting "<< "[" << std::to_string(edge[0]) << ", " << std::to_string(edge[1]) << "], " << std::to_string(alpha) << "\n";
                else {
                    face_debug << "    Visiting " << std::to_string(edge[(unsigned char) alpha]) << '\n';
                }
#endif

                if (alpha > 0 and alpha < 1) {
                    // Current position is an intersection

                    if (lastPositionIsAnIntersection) {
                        // If the last position was also an intersection, we are coming from a triangle's edge and
                        // the next position should be a node

                        const auto & edge = edges[edgeId];
                        const NodeId & node1 = edge[0];
                        const NodeId & node2 = edge[1];

                        // We find the next node this way to make sure that our order of nodes in a face match the
                        // order to get an outward normal to the surface.
                        char currentFaceNode = -1;
                        for (unsigned char i = 0; i < 4; ++i) {
                            if ((face_nodes[i] == node1 && face_nodes[(i+1)%4] == node2)
                                || (face_nodes[i] == node2 && face_nodes[(i+1)%4] == node1)) {
                                currentFaceNode = (char) ((i+1)%4);
                            }
                        }

                        assert_trace ((currentFaceNode > -1), "Failed to match the two nodes of the current cut edge", (cube_debug.str() + face_debug.str()));
                        assert_trace (isInside[face_nodes[currentFaceNode]], "The next node adjacent to an intersection should be inside the surface.", (cube_debug.str() + face_debug.str()));

                        next_node_id = face_nodes[currentFaceNode];
                        subface.push_back(alpha_position(next_node_id, edgeId));

                    } else {
                        // If the last position was a node, we are following the triangle's edge to the next
                        // intersection on the same face

                        // First, check all triangles connected to this intersection and find the one that connects
                        // another intersection on the same face
                        std::vector<AlphaPosition> connectingIntersection;
                        for (const AlphaTriangle & triangle : edge_triangles[edgeId]) {
                            bool addTriangle = false;
                            for (const AlphaPosition & intersection : triangle) {
                                if (intersection.first != edgeId) {
                                    for (const EdgeId & e : face_edges) {
                                        if (e == intersection.first) {
                                            addTriangle = true;
                                            connectingIntersection.push_back(intersection);
                                            break;
                                        }
                                    }
                                }
                                if (addTriangle)
                                    break;
                            }
                        }

                        assert_trace(connectingIntersection.size() == 1, "There should be exactly one triangle per intersection that connects to the same face.", (cube_debug.str() + face_debug.str()));

                        subface.push_back(connectingIntersection[0]);
#ifndef	NDEBUG
                        face_debug << "      Pushing (" <<std::to_string(edges[subface.back().first][0]) <<", "<< std::to_string(edges[subface.back().first][1]) << ", " << std::to_string(subface.back().second) << ")" << '\n';
#endif
                    }

                    lastPositionIsAnIntersection = true;
                } else {
                    // Current position is a node, follow either the next adjacent node, or the intersection if there is one connected to it

                    node_id = edges[edgeId][(unsigned char) alpha];
                    char currentFaceNode = -1;
                    for (unsigned char i = 0; i < 4; ++i) {
                        if (face_nodes[i] == node_id)
                            currentFaceNode = i;
                    }

                    assert_trace(currentFaceNode > -1, "Cannot find the current face node", (cube_debug.str() + face_debug.str()));

                    if (node_visited[currentFaceNode]) {
                        subface.pop_back();
                        subfaceIsClosed = true;
                    } else {
                        lastPositionIsAnIntersection = false;
                        node_visited[currentFaceNode] = true;


                        next_node_id = face_nodes[(currentFaceNode + 1) % 4];
                        char next_next_node_id = face_nodes[(currentFaceNode + 2) % 4];

#ifndef    NDEBUG
                        face_debug << "      Node_id = " << std::to_string(node_id) << ", next_node_id = "
                                   << std::to_string(next_node_id) << ", next_next_node_id = "
                                   << std::to_string(next_next_node_id) << '\n';
#endif

                        currentEdgeId = find_edge_id_from_two_nodes(node_id, next_node_id);
                        const AlphaValue & intersection_alpha = intersection_of_edge[currentEdgeId];
                        if (intersection_alpha > 0 and intersection_alpha < 1) {
                            // We follow the intersection
                            assert_trace (not isInside[next_node_id], "The adjacent node of this node should be outside since there is an intersection between them.", (cube_debug.str() + face_debug.str()));
                            subface.push_back(std::make_pair(edgeId, intersection_alpha));
                        } else if (isInside[next_node_id]) {
                            // We follow the next node
                            currentEdgeId = find_edge_id_from_two_nodes(node_id, next_node_id);
                            assert_trace(
                                    (!connected_intersection[currentFaceNode] || currentEdgeId != connected_intersection[currentFaceNode]->first),
                                    "The adjacent node of this node should be outside since there it is connected to an intersected edge.",
                                    (cube_debug.str() + face_debug.str())
                            );

                            node_id = next_node_id;
                            next_node_id = next_next_node_id;
                            currentEdgeId = find_edge_id_from_two_nodes(node_id, next_node_id);
                            assert_trace (currentEdgeId >= 0, "An edge must exists between those nodes.",
                                          (cube_debug.str() + face_debug.str()));

                            subface.push_back(alpha_position(node_id, currentEdgeId));
#ifndef    NDEBUG
                            face_debug << "      Pushing (" << std::to_string(edges[subface.back().first][0])
                                       << ", " << std::to_string(edges[subface.back().first][1]) << ", "
                                       << std::to_string(subface.back().second) << ")" << '\n';
#endif
                        } else if (connected_intersection[currentFaceNode]) {
                            // We follow the intersection
                            subface.push_back(*connected_intersection[currentFaceNode]);
#ifndef    NDEBUG
                            face_debug << "      Pushing (" << std::to_string(edges[subface.back().first][0])
                                       << ", " << std::to_string(edges[subface.back().first][1]) << ", "
                                       << std::to_string(subface.back().second) << ")" << '\n';
#endif
                        } else {
                            subfaceIsClosed = true;
                        }
                    }
                }
            }

            if (subface.size() > 2) // Make sure that it is a complete surface and not only an edge
                subfaces.push_back(subface);
        }

        // Create triangles for each subfaces found
        std::vector<std::array<AlphaPosition, 3>> alpha_triangles;
        for (std::list<AlphaPosition> & subface : subfaces) {

            assert_trace(subface.size() <= 5, "A subface cannot contain more than 5 nodes", (cube_debug.str() + face_debug.str()));

            AlphaTriangle t1, t2, t3;
            t1[0] = subface.front(); subface.pop_front();
            t1[1] = subface.front(); subface.pop_front();
            t1[2] = subface.front(); subface.pop_front();
            alpha_triangles.push_back(t1);

            if (!subface.empty()) {
                t2[0] = t1[2];
                t2[1] = subface.front(); subface.pop_front();
                t2[2] = t1[0];
                alpha_triangles.push_back(t2);
            }

            if (!subface.empty()) {
                t3[0] = t2[2];
                t3[1] = subface.front(); subface.pop_front();
                t3[2] = t1[0];
                alpha_triangles.push_back(t3);
            }
        }

        triangles.insert(std::end(triangles), std::begin(alpha_triangles), std::end(alpha_triangles));
    }

    return triangles;
}

/**
 * Triangulate an hexahedron cut by an iso-surface. This function will create a mesh of triangles that cover the iso-surface inside an hexahedron.
 * @tparam IsInsideFct bool IsInsideFct(unsigned char node_id)
 * @tparam IntersectFct Real IntersectFct(unsigned char edge_id)
 * @param IsInside [IN] Callback function that returns true if the vertex passed by argument is inside or outside the iso-surface.
 * @param ComputeIntersection [IN] Callback function that returns the alpha position (between 0 and 1) of the intersection of the iso-surface and the edge passed by argument.
 * @param triangles [OUT] List of triangles where a triangle is represented by three edge intersections (an instersection is a pair of <edge_id, alpha-position> (see @param ComputeIntersection).
 * @param snap_treshold Distance threshold where two point should be snapped together into one point.
 * @return True if the hexahedron was triangulated, false otherwise (the hexahedron is completely inside or outside of the iso-surface).
 */
template<typename IntersectFct, typename IsInsideFct>
bool triangulate (const IsInsideFct & IsInside, const IntersectFct & ComputeIntersection, std::vector<std::array<AlphaPosition, 3>> & triangles, const double snap_treshold = 0.000000000001)
{

    /*
      Determine the index into the edge table which
      tells us which vertices are inside of the surface
      */
    std::array<bool, 8> isInside = {{false, false, false, false, false, false, false, false}};
    unsigned short flags = 0;
    for (unsigned char i = 0; i < 8; ++i)
        if (IsInside(i)) {
            flags |= (1 << i);
            isInside[i] = true;
        }

    /* Cube is entirely in/out of the surface */
    if ( MarchingCubeEdgeTable[flags] == 0 )
        return false;

    /* Find the positions where the surface intersects the edges */
    std::array<AlphaValue, 12> edgeIntersections = {{-1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1 ,-1}};
    for (unsigned char edge_id = 0; edge_id < 12; ++edge_id) {
        if (MarchingCubeEdgeTable[flags] & (1 << edge_id)) {
            // Linearly interpolate the position where an isosurface cuts an edge between two vertices, each with their
            // own scalar value
            edgeIntersections[edge_id] = ComputeIntersection(edge_id);
        }
    }

    for (std::array<AlphaPosition, 3> & triangle : triangles) {
        for (AlphaPosition node : triangle) {
            AlphaValue & alpha = node.second;
            if (alpha < snap_treshold)
                alpha = 0;
            else if (1. - alpha < snap_treshold || alpha - 1. < snap_treshold )
                alpha = 1;
        }
    }

    std::vector<std::array<AlphaPosition, 3>> t = triangulate(isInside, edgeIntersections);
    triangles.insert(std::end(triangles), std::begin(t), std::end(t));

    return true;
}

/**
 * Triangulate an hexahedron cut by an iso-surface. This function will create a mesh of triangles that cover the iso-surface inside an hexahedron.
 * @tparam Coord Coordinate type, should implement the + and * operators
 * @tparam IntersectFct Real IntersectFct(unsigned char edge_id)
 * @tparam IsInsideFct bool IsInsideFct(unsigned char node_id)
 * @param nodes [IN] Array of the hexahedron's nodes positions.
 * @param IsInside [IN] Callback function that returns true if the vertex passed by argument is inside or outside the iso-surface.
 * @param ComputeIntersection [IN] Callback function that returns the alpha position (between 0 and 1) of the intersection of the iso-surface and the edge passed by argument.
 * @param snap_treshold Distance threshold where two point should be snapped together into one point.
 * @return A list of triangles forming the mesh. A triangle is represented by an array of 3 coordinate.
 */
template<class Coord, typename IntersectFct, typename IsInsideFct>
std::vector<std::array<Coord, 3>> triangulate (const std::array<Coord, 8> & nodes, const IsInsideFct & IsInside, const IntersectFct & ComputeIntersection, const double snap_treshold = 0.000000000001)
{
    std::vector<std::array<Coord, 3>> triangles;

    std::vector<std::array<AlphaPosition, 3>> triangles_by_intersection;
    bool isboundary = triangulate(IsInside, ComputeIntersection, triangles_by_intersection, snap_treshold);

    // Cube is entirely in/out of the surface, returns an empty set of triangles
    if (isboundary) {
        for (const auto &triangle_by_intersection : triangles_by_intersection) {
            std::array<Coord, 3> triangle;
            for (unsigned char i = 0; i < 3; ++i) {
                const auto edge_id = triangle_by_intersection[i].first;
                const AlphaValue alpha = triangle_by_intersection[i].second;
                const auto &edge = edges[edge_id];

                const auto &p0 = nodes[edge[0]];
                const auto &p1 = nodes[edge[1]];

                triangle[i] = p0 + alpha * (p1 - p0);
            }
            triangles.push_back(triangle);
        }
    }

    return triangles;
}


template<class Coord, typename IntersectFct, typename IsInsideFct>
std::vector<std::array<Coord, 3>> triangulate_interior (const std::array<Coord, 8> & nodes, const IsInsideFct & IsInside, const IntersectFct & ComputeIntersection, const double snap_treshold = 0.000000000001)
{
    using AlphaTriangle = std::array<AlphaPosition, 3>;
    using Triangle = std::array<Coord, 3>;

    std::vector<Triangle> triangles;

    std::vector<AlphaTriangle> triangles_by_intersection;
    triangulate(IsInside, ComputeIntersection, triangles_by_intersection, snap_treshold);

    // Add the triangles while converting alpha-position back to coordinates
    for (const AlphaTriangle & alpha_triangle : triangles_by_intersection) {
        Triangle triangle;
        for (unsigned char i = 0; i < 3; ++i) {
            const Coord p0 = nodes[edges[alpha_triangle[i].first][0]];
            const Coord p1 = nodes[edges[alpha_triangle[i].first][1]];
            const AlphaValue alpha = alpha_triangle[i].second;
            triangle[i] = p0 + (p1-p0)*alpha;
        }
        triangles.push_back(triangle);
    }

    return triangles;
}

} // namespace hexahedron

} // namespace helper

} // namespace caribou

} // namespace sofa

#endif //CARIBOU_HELPER_HEXAHEDRON_H
