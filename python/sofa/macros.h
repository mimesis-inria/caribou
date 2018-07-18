#ifndef CARIBOU_PYTHON_SOFA_MACROS_H
#define CARIBOU_PYTHON_SOFA_MACROS_H

#define CP_MODULE_METHODS_BEGIN(MODULENAME) static PyMethodDef module##_##MODULENAME##_##methods[] = {
#define CP_MODULE_METHOD(MODULENAME,M) {#M, M, METH_VARARGS, ""},
#define CP_MODULE_METHOD_DOC(MODULENAME,M, D) {#M, M, METH_VARARGS, D},
#define CP_MODULE_METHOD_KW(MODULENAME,M) {#M, (PyCFunction)M, METH_KEYWORDS|METH_VARARGS, ""},
#define CP_MODULE_METHOD_KW_DOC(MODULENAME,M, D) {#M, (PyCFunction)M, METH_KEYWORDS|METH_VARARGS, D},
#define CP_MODULE_METHODS_END {NULL,NULL,0,NULL} };

#define CP_MODULE_DOC(MODULENAME, DOC) static const char * module##_##MODULENAME##_##doc = DOC;

#define CP_DECLARE_MODULE(MODULENAME) \
    extern PyMethodDef * module##_##MODULENAME##_##methods;\
    extern const char * module##_##MODULENAME##_##doc;

#define CP_CREATE_MODULE(MODULENAME) \
    Py_InitModule3(#MODULENAME, module##_##MODULENAME##_##methods, module##_##MODULENAME##_##doc);

#endif //CARIBOU_PYTHON_SOFA_MACROS_H
