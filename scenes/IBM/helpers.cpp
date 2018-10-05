
#undef NDEBUG

#include "../../plugin/Helper/Hexahedron.h"
#include "../../plugin/Helper/MirtichIntegration.h"

#include <Python.h>

#include <array>

static PyObject *
helpers_triangulate(PyObject *self, PyObject *args)
{
    PyObject *pIsInside;
    PyObject *pEdgeIntersections;
    PyObject *pItem;
    Py_ssize_t n;
    int i;

    // Parse the parameters
    if (!PyArg_ParseTuple(args, "O!O!", &PyList_Type, &pIsInside, &PyList_Type, &pEdgeIntersections)) {
        PyErr_SetString(PyExc_TypeError, "parameter must be a list.");
        return nullptr;
    }

    // Check the size of the isInside array
    n = PyList_Size(pIsInside);
    if (n != 8) {
        PyErr_SetString(PyExc_TypeError, "The first parameter must be a list of 8 boolean identifying if the given node is inside or outside of the iso-surface");
    }

    // Parse the isInside array
    std::array<bool, 8> isInside;
    for (i=0; i<n; i++) {
        pItem = PyList_GetItem(pIsInside, i);
        if(!PyBool_Check(pItem) and !PyInt_Check(pItem)) {
            PyErr_SetString(PyExc_TypeError, "The IsInside list must contain only boolean.");
            return nullptr;
        }
        isInside[i] = PyObject_IsTrue(pItem);
    }

    // Check the size of the edgeIntersections array
    n = PyList_Size(pEdgeIntersections);
    if (n != 12) {
        PyErr_SetString(PyExc_TypeError, "The second parameter must be a list of 12 float identifying where the iso-surface cut an edge.\n"
                                         " Each coefficient should be between 0 and 1, 0 being the position of the first node of\n"
                                         " the edge, 1 being the second node, and anywhere between is the exact location of the cut.");
    }

    // Parse the edgeIntersections array
    std::array<sofa::caribou::helper::hexahedron::AlphaValue, 12> edgeIntersections;
    for (i=0; i<n; i++) {
        pItem = PyList_GetItem(pEdgeIntersections, i);
        if(!PyFloat_Check(pItem) and !PyInt_Check(pItem)) {
            PyErr_SetString(PyExc_TypeError, "The edgeIntersections argument must be contain only float values between 0 and 1 if an intersection is found, or -1 if no intersection is on the edge.");
            return nullptr;
        }

        if (PyInt_Check(pItem))
            edgeIntersections[i] = (float) PyInt_AsLong(pItem);
        else
            edgeIntersections[i] = PyFloat_AsDouble(pItem);
    }

    // Get the triangles
    using AlphaPosition = sofa::caribou::helper::hexahedron::AlphaPosition;
    using AlphaTriangle = std::array<AlphaPosition, 3>;
    std::array<std::vector<AlphaTriangle>, 6> faces_triangles;
    std::vector<AlphaTriangle> cut_triangles;

    sofa::caribou::helper::hexahedron::triangulate(isInside, edgeIntersections, faces_triangles, cut_triangles);

    // Serialize the result to python lists
    PyObject * pFaceTriangles = PyList_New(6);
    for (unsigned int face_id = 0; face_id < 6; ++ face_id) {
        const std::vector<AlphaTriangle> & triangles = faces_triangles[face_id];
        PyObject * pTriangles = PyList_New(triangles.size());
        for (unsigned int i = 0; i < triangles.size(); ++i) {
            PyObject *pTriangleNodes = PyList_New(3);
            for (unsigned int j = 0; j < 3; ++j) {
                sofa::caribou::helper::hexahedron::EdgeId edgeId = triangles[i][j].first;
                sofa::caribou::helper::hexahedron::AlphaValue alpha = triangles[i][j].second;
                PyObject *pAlphaPosition = PyTuple_Pack(
                        3,
                        PyInt_FromLong(sofa::caribou::helper::hexahedron::edges[edgeId][0]),
                        PyInt_FromLong(sofa::caribou::helper::hexahedron::edges[edgeId][1]),
                        PyFloat_FromDouble(alpha)
                );
                PyList_SetItem(pTriangleNodes, j, pAlphaPosition);
            }
            PyList_SetItem(pTriangles, i, pTriangleNodes);
        }
        PyList_SetItem(pFaceTriangles, face_id, pTriangles);
    }

    PyObject * pCutTriangles = PyList_New(cut_triangles.size());
    for (unsigned int i = 0; i < cut_triangles.size(); ++i) {
        PyObject *pTriangleNodes = PyList_New(3);
        for (unsigned int j = 0; j < 3; ++j) {
            sofa::caribou::helper::hexahedron::EdgeId edgeId = cut_triangles[i][j].first;
            sofa::caribou::helper::hexahedron::AlphaValue alpha = cut_triangles[i][j].second;
            PyObject *pAlphaPosition = PyTuple_Pack(
                    3,
                    PyInt_FromLong(sofa::caribou::helper::hexahedron::edges[edgeId][0]),
                    PyInt_FromLong(sofa::caribou::helper::hexahedron::edges[edgeId][1]),
                    PyFloat_FromDouble(alpha)
            );
            PyList_SetItem(pTriangleNodes, j, pAlphaPosition);
        }
        PyList_SetItem(pCutTriangles, i, pTriangleNodes);
    }

    return PyTuple_Pack(2, pFaceTriangles, pCutTriangles);
}

static PyObject *
helpers_integrate(PyObject *self, PyObject *args)
{
    PyObject *pFaces;
    // Parse the parameters
    if (!PyArg_ParseTuple(args, "O!", &PyList_Type, &pFaces)) {
        PyErr_SetString(PyExc_TypeError, "parameter must be a list.");
        return nullptr;
    }

    // get the faces and their subfaces
    std::vector<std::vector<std::array<double, 3>>> faces;
    Py_ssize_t nb_faces = PyList_Size(pFaces);
    for (unsigned int i = 0; i < nb_faces; ++i) {
        std::vector<std::array<double, 3>> subfaces;
        PyObject *pSubFaces = PyList_GetItem(pFaces, i);
        if (!PyList_Check(pSubFaces)) {
            PyErr_SetString(PyExc_TypeError, "parameter must be a list of faces where each faces is a list of subfaces.");
            return nullptr;
        }

        Py_ssize_t nb_subfaces = PyList_Size(pSubFaces);
        for (unsigned int j = 0; j < nb_subfaces; ++j) {
            std::array<double, 3> pos;
            PyObject * pPos = PyList_GetItem(pSubFaces, j);
            if (PyList_Check(pPos) && PyList_Size(pPos) == 3) {
                pos[0] = PyFloat_AsDouble(PyList_GetItem(pPos, 0));
                pos[1] = PyFloat_AsDouble(PyList_GetItem(pPos, 1));
                pos[2] = PyFloat_AsDouble(PyList_GetItem(pPos, 2));
            } else if(PyTuple_Check(pPos) && PyTuple_Size(pPos) == 3) {
                pos[0] = PyFloat_AsDouble(PyTuple_GetItem(pPos, 0));
                pos[1] = PyFloat_AsDouble(PyTuple_GetItem(pPos, 1));
                pos[2] = PyFloat_AsDouble(PyTuple_GetItem(pPos, 2));
            } else {
                PyErr_SetString(PyExc_TypeError, "Subfaces must be a list of positions where a position is a list of 3 floats or a tuple of 3 floats");
                return nullptr;
            }
            subfaces.push_back(pos);
        }
        faces.push_back(subfaces);
    }

    std::vector<double> result = integrate<std::array<double, 3>>(faces);

    PyObject * pResults = PyList_New(result.size());
    for (unsigned int i = 0; i < result.size(); ++i) {
        PyObject * pValue = PyFloat_FromDouble(result[i]);
        PyList_SetItem(pResults, i, pValue);
    }

    return pResults;

}

static PyMethodDef HexahedronMethods[] = {
        {"triangulate",  helpers_triangulate, METH_VARARGS, "Triangulate an hexahedron cut by an iso-surface. This function will create a mesh of triangles that cover the iso-surface inside an hexahedron (the cut surface)."},
        {NULL, NULL, 0, NULL}        /* Sentinel */
};

static PyMethodDef HelperMethods[] = {
        {"integrate",  helpers_integrate, METH_VARARGS, "Integrate an hexahedron cut by an iso-surface."},
        {NULL, NULL, 0, NULL}        /* Sentinel */
};

extern "C" void initHelpers(void)
{
    if (!Py_IsInitialized()) {
        Py_Initialize();
    }

    PyEval_InitThreads();

    PyObject *libraryModule = Py_InitModule("Helpers", HelperMethods);
    PyObject *hexahedronModule = Py_InitModule("Hexahedron", HexahedronMethods);

    PyModule_AddObject(libraryModule, "Hexahedron", hexahedronModule);
}