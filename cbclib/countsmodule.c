#include <Python.h>

#include "countscalc.h"

const int MAX_LENGTH = 14;
const int MAX_GAP_LENGTH = 14;

typedef struct{
    int slen, pos, glen;
    int hash;
} SiteStruct;

static void fill_struct_hash(SiteStruct *s){
    int struct_hash = 0;
    struct_hash += s->slen;
    if(s->glen){
        struct_hash += MAX_LENGTH + 1;
        struct_hash += (s->pos - 1) * MAX_LENGTH;
        struct_hash += (s->glen - 1) * (MAX_LENGTH - 1) * MAX_LENGTH;
    }
    s->hash = struct_hash;
}

typedef struct{
    PyObject_HEAD
    long **countset;
    long *totals;
    SiteStruct *structs;
    int structs_num;
} Counts;

static PyObject *Counts_new(PyTypeObject *type, PyObject *args,
                            PyObject *kwds){
    Counts *self;
    self = (Counts *)type->tp_alloc(type, 0);
    return (PyObject *)self;
}

static int Counts_init(Counts *self, PyObject *args, PyObject *kwds){
    char *filename; // leak of filename; uncomment the last
    PyObject *structs_seq;
    if(!PyArg_ParseTuple(args, "sO", &filename, &structs_seq))
        return -1;
    else{
        int structs_num = PySequence_Size(structs_seq);
        int i;
        SiteStruct *structs;
        structs = (SiteStruct *)calloc(structs_num, sizeof(SiteStruct));
        int max_hash = 0;
        for(i=0; i<structs_num; i++){
            PyObject *item = PySequence_GetItem(structs_seq, i);
            SiteStruct *sp = &(structs[i]);
            if(!PyArg_ParseTuple(item, "iii", &sp->slen, &sp->pos,
                                 &sp->glen))
                return -1;
            fill_struct_hash(sp);
            if(sp->hash > max_hash)
                max_hash = sp->hash;
        }
        // reuse of __init__ will cause memmory leak!! #TOFIX
        self->countset = (long **)calloc(max_hash+1, sizeof(long *));
        self->totals = (long *)calloc(max_hash+1, sizeof(long));
        for(i=0; i<structs_num; i++){
            SiteStruct s = structs[i];
            long **counts;
            int j;
            if(s.glen){
                counts = count_bipart_words(filename, s.slen,
                                            s.pos, s.glen);
                long total_index = 1ul << (2 * s.slen);
                int num = s.slen-s.pos;
                for(j=0; j<num; j++){
                    self->countset[s.hash-j] = counts[j];
                    self->totals[s.hash-j] = *(counts[j] + total_index);
                    total_index >>= 2;
                }
            }else{
                counts = count_short_words(filename, s.slen);
                long total_index = 1ul << (2 * s.slen);
                for(j=0; j<s.slen; j++){
                    self->countset[s.slen-j] = counts[j];
                    self->totals[s.slen-j] = *(counts[j] + total_index);
                    total_index >>= 2;
                }
            }
            free(counts);
        }
        self->structs = structs;
        self->structs_num = structs_num;
        //free(filename);
    }
    return 0;
}

static void Counts_dealloc(Counts *self){
    SiteStruct *structs;
    structs = self->structs;
    int i;
    for(i=0; i<self->structs_num; i++){
        SiteStruct s = structs[i];
        long **ptr;
        ptr = self->countset + s.hash;
        while(s.slen > s.pos){
            free(*ptr);
            ptr --;
            s.slen --;
        }
    }
    free(structs);
    free(self->countset);
    free(self->totals);
    self->ob_type->tp_free((PyObject *)self);
}

static PyObject *Counts_get_count(Counts *self, PyObject *args,
                                  PyObject *kwds){
    PyObject *dsite_seq;
    int struct_hash = 1;
    if(!PyArg_ParseTuple(args, "O|i", &dsite_seq, &struct_hash))
        return NULL;
    else{
        long *counts;
        counts = self->countset[struct_hash];
        long count = 0;
        int dsite_num = PySequence_Size(dsite_seq);
        long dsite;
        int i;
        for(i=0; i<dsite_num; i++){
            dsite = PyInt_AsLong(PySequence_GetItem(dsite_seq, i));
            count += counts[dsite];
        }
        return Py_BuildValue("l", count);
    }
}

static PyObject *Counts_get_freq(Counts *self, PyObject *args,
                                 PyObject *kwds){
    PyObject *dsite_seq;
    int struct_hash = 1;
    if(!PyArg_ParseTuple(args, "O|i", &dsite_seq, &struct_hash))
        return NULL;
    else{
        long *counts;
        counts = self->countset[struct_hash];
        long count = 0;
        int dsite_num = PySequence_Size(dsite_seq);
        long dsite;
        int i;
        for(i=0; i<dsite_num; i++){
            dsite = PyInt_AsLong(PySequence_GetItem(dsite_seq, i));
            count += counts[dsite];
        }
        double freq = (double) count;
        long total = self->totals[struct_hash];
        // if(total == 0)
        return Py_BuildValue("d", freq/total);
    }
}

static PyObject *Counts_get_total(Counts *self, PyObject *args,
                                  PyObject *kwds){
    int struct_hash = 1;
    if(!PyArg_ParseTuple(args, "|i", &struct_hash))
        return NULL;
    else{
        long total = self->totals[struct_hash];
        return Py_BuildValue("l", total);
    }
}

static PyMethodDef Counts_methods[] = {
    {
        "get_count", (PyCFunction)Counts_get_count, METH_VARARGS,
        "get_count(dsite, [struct_hash]) -> int"
    }, {
        "get_freq", (PyCFunction)Counts_get_freq, METH_VARARGS,
        "get_freq(dsite, [struct_hash]) -> float"
    }, {
        "get_total", (PyCFunction)Counts_get_total, METH_VARARGS,
        "get_total([struct_hash]) -> int"
    }, {NULL}
};

static PyTypeObject CountsType = {
    PyObject_HEAD_INIT(NULL)
    0,                          /*ob_size*/
    "counts.Counts",            /*tp_name*/
    sizeof(Counts),             /*tp_basicsize*/
    0,                          /*tp_itemsize*/
    (destructor)Counts_dealloc, /*tp_dealloc*/
    0,                          /*tp_print*/
    0,                          /*tp_getattr*/
    0,                          /*tp_setattr*/
    0,                          /*tp_compare*/
    0,                          /*tp_repr*/
    0,                          /*tp_as_number*/
    0,                          /*tp_as_sequence*/
    0,                          /*tp_as_mapping*/
    0,                          /*tp_hash*/
    0,                          /*tp_call*/
    0,                          /*tp_str*/
    0,                          /*tp_getattro*/
    0,                          /*tp_setattro*/
    0,                          /*tp_as_buffer*/
    Py_TPFLAGS_DEFAULT | Py_TPFLAGS_BASETYPE, /*tp_flags*/
    "Counts(filename, struct_tuple)",         /*tp_doc*/
    0,                          /*tp_traverse*/
    0,                          /*tp_clear*/
    0,                          /*tp_richcompare*/
    0,                          /*tp_weaklistoffset*/
    0,                          /*tp_iter*/
    0,                          /*tp_iternext*/
    Counts_methods,             /*tp_methods*/
    0,                          /*tp_members*/
    0,                          /*tp_getset*/
    0,                          /*tp_base*/
    0,                          /*tp_dict*/
    0,                          /*tp_descr_get*/
    0,                          /*tp_descr_set*/
    0,                          /*tp_dictoffset*/
    (initproc)Counts_init,      /*tp_init*/
    0,                          /*tp_alloc*/
    Counts_new,                 /*tp_new*/
};

static PyMethodDef counts_methods[] = {
    {NULL}
};

#ifndef PyMODINIT_FUNC
#define PyMODINIT_FUNC void
#endif

PyMODINIT_FUNC
initcounts(void){
    PyObject* counts;
    CountsType.tp_new = PyType_GenericNew;
    if (PyType_Ready(&CountsType) < 0)
        return;
    counts = Py_InitModule3("counts", counts_methods, "module help");
    Py_INCREF(&CountsType);
    PyModule_AddObject(counts, "Counts", (PyObject *)&CountsType);
}
