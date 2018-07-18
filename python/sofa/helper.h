#ifndef CARIBOU_PYTHON_SOFA_HELPER_H
#define CARIBOU_PYTHON_SOFA_HELPER_H

#include <Python.h>
#include <vector>
#include <string>

namespace caribou {
namespace python {
namespace sofa {
namespace helper {

// =============================================================================
// Ptr objects passed to python
// deletion can be made by Python IF the "deletable" flag is true,
// and if a FreeFunc is provided in the class definition
// =============================================================================
template <class T>
struct PyPtr
{
    PyObject_HEAD
    T* object;
    bool deletable;
};

template <class T>
class ClassCreator {
public:

    // Functions
    inline PyTypeObject* create_type() {

    }

    // Members
    std::string name = "";
    std::string doc = "";
    std::vector<PyMethodDef*> methods;
    cmpfunc* compare = nullptr;
    getattrfunc* getattribute = nullptr;
    setattrfunc* setattribute = nullptr;
    printfunc* print = nullptr;
    destructor* destructor = nullptr;
    initproc* constructor = nullptr;
    PyTypeObject* parent = nullptr;

};

} // namespace helper
} // namespace sofa
} // namespace python
} // namespace caribou

#endif //CARIBOU_PYTHON_SOFA_HELPER_H
