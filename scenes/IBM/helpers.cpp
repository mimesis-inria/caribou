
#undef NDEBUG

#include "../../plugin/Helper/Hexahedron.h"

#include <Python.h>

#include <array>

static PyObject *
triangulate(PyObject *self, PyObject *args)
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
    std::vector<std::array<sofa::caribou::helper::hexahedron::AlphaPosition, 3>> triangles = sofa::caribou::helper::hexahedron::triangulate(isInside, edgeIntersections);

    // Serialize the result to python lists
    PyObject * pTriangles = PyList_New(triangles.size());
    for (unsigned int i = 0; i < triangles.size(); ++i) {
        PyObject * pTriangleNodes = PyList_New(3);
        for (unsigned int j = 0; j < 3; ++j) {
            sofa::caribou::helper::hexahedron::EdgeId edgeId = triangles[i][j].first;
            sofa::caribou::helper::hexahedron::AlphaValue alpha = triangles[i][j].second;
            PyObject * pAlphaPosition = PyTuple_Pack(
                3,
                PyInt_FromLong(sofa::caribou::helper::hexahedron::edges[edgeId][0]),
                PyInt_FromLong(sofa::caribou::helper::hexahedron::edges[edgeId][1]),
                PyFloat_FromDouble(alpha)
            );
            PyList_SetItem(pTriangleNodes, j, pAlphaPosition);
        }
        PyList_SetItem(pTriangles, i, pTriangleNodes);
    }

    return pTriangles;
}

static PyMethodDef HexahedronMethods[] = {
        {"triangulate",  triangulate, METH_VARARGS, "Triangulate an hexahedron cut by an iso-surface. This function will create a mesh of triangles that cover the iso-surface inside an hexahedron (the cut surface)."},
        {NULL, NULL, 0, NULL}        /* Sentinel */
};

extern "C" void initHelpers(void)
{
    if (!Py_IsInitialized()) {
        Py_Initialize();
    }

    PyEval_InitThreads();

    PyObject *libraryModule = Py_InitModule("Helpers", nullptr);
    PyObject *hexahedronModule = Py_InitModule("Hexahedron", HexahedronMethods);

    PyModule_AddObject(libraryModule, "Hexahedron", hexahedronModule);
}