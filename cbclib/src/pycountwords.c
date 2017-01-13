#include <Python.h>

#include "count_words.h"


typedef struct{
    PyObject_HEAD
    
} pcw_Counts;

static PyTypeObject pcw_CountsType = {
    PyObject_HEAD_INIT(NULL)
    0,                         /*ob_size*/
    "countwords.Counts",       /*tp_name*/
    sizeof(pcw_Counts),        /*tp_basicsize*/
    0,                         /*tp_itemsize*/
    0,                         /*tp_dealloc*/
    0,                         /*tp_print*/
    0,                         /*tp_getattr*/
    0,                         /*tp_setattr*/
    0,                         /*tp_compare*/
    0,                         /*tp_repr*/
    0,                         /*tp_as_number*/
    0,                         /*tp_as_sequence*/
    0,                         /*tp_as_mapping*/
    0,                         /*tp_hash */
    0,                         /*tp_call*/
    0,                         /*tp_str*/
    0,                         /*tp_getattro*/
    0,                         /*tp_setattro*/
    0,                         /*tp_as_buffer*/
    Py_TPFLAGS_DEFAULT,        /*tp_flags*/
    "Counts help",             /* tp_doc */
};

static PyMethodDef Counts_methods[] = {
    {NULL}
};

#ifndef PyMODINIT_FUNC
#define PyMODINIT_FUNC void
#endif

PyMODINIT_FUNC
initcounts(void)
{
    PyObject* counts;

    pcw_CountsType.tp_new = PyType_GenericNew;
    if (PyType_Ready(&pcw_CountsType) < 0)
        return;

    counts = Py_InitModule3("Counts", Counts_methods, "Counts help");

    Py_INCREF(&pcw_CountsType);
    PyModule_AddObject(counts, "Counts", (PyObject *)&pcw_CountsType);
}

static PyObject *py_count_words(PyObject *self, PyObject *args){
    char *filename;
    int len, pos = 0, gap = 0;
    if(!PyArg_Parse(args, "si|ii", &filename, len, pos, gap))
        return NULL;
    else{
        return NULL;
    }
}

static PyMethodDef pcw_methods[] = {
    {"count_words", py_count_words, METH_VARARGS, "count words help"},
    {NULL, NULL, 0, NULL}
};

static struct PyModuleDef pcw_module = {
    PyModuleDef_HEAD_INIT,
    "countwords", "module help", -1, pcw_methods
};

PyMODINIT_FUNC
PyInit_countwords(){
    return PyModule_Create(&pcw_module);
}
