#include <Python.h>

#include "countscalc.h"

typedef struct{
    PyObject_HEAD
    long **counts;
} counts_Counts;

static PyTypeObject counts_CountsType = {
    PyObject_HEAD_INIT(NULL)
    0,                         /*ob_size*/
    "counts.Counts",           /*tp_name*/
    sizeof(counts_Counts),     /*tp_basicsize*/
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
    0,                         /*tp_hash*/
    0,                         /*tp_call*/
    0,                         /*tp_str*/
    0,                         /*tp_getattro*/
    0,                         /*tp_setattro*/
    0,                         /*tp_as_buffer*/
    Py_TPFLAGS_DEFAULT,        /*tp_flags*/
    "Counts help",             /*tp_doc*/
};

static PyObject *Counts_count(PyObject *counts, PyObject *args){
    char *filename;
    int len, pos = 0, gap = 0;
    if(!PyArg_Parse(args, "si|ii", &filename, len, pos, gap))
        return NULL;
    else{
        return NULL;
    }
}

static PyMethodDef counts_methods[] = {
    {NULL}
};

#ifndef PyMODINIT_FUNC
#define PyMODINIT_FUNC void
#endif

PyMODINIT_FUNC
initCounts(void)
{
    PyObject* counts;

    counts_CountsType.tp_new = PyType_GenericNew;
    if (PyType_Ready(&counts_CountsType) < 0)
        return;

    counts = Py_InitModule3("counts", counts_methods, "module help");

    Py_INCREF(&counts_CountsType);
    PyModule_AddObject(counts, "Counts", (PyObject *)&counts_CountsType);
}

