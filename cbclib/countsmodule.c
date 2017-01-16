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
    long **totals;
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
    char *filename;
    PyObject *structs_tuple;
    if(!PyArg_ParseTuple(args, "sO", &filename, &structs_tuple))
        return NULL;
    else{
        int structs_num = PyTuple_Size(structs_tuple);
        int i;
        SiteStruct *structs;
        structs = (SiteStruct *)calloc(structs_num, sizeof(SiteStruct));
        int max_hash = 0;
        for(i=0; i<structs_num; i++){
            PyObject *item = PyTuple_GetItem(structs_tuple, i);
            SiteStruct s = structs[i];
            if(!PyArg_ParseTuple(item, "iii", &s.slen, &s.pos, &s.glen))
                continue;
            fill_struct_hash(&s);
            if(s.hash > max_hash)
                max_hash = s.hash;
        }
        self->countset = (long **)calloc(max_hash+1, sizeof(long *));
        self->totals = (long **)calloc(max_hash+1, sizeof(long *));
        for(i=0; i<structs_num; i++){
            SiteStruct s = structs[i];
            long **counts;
            int j;
            if(s.glen){
                counts = count_short_words(filename, s.slen);
                long total_index = 1ul << (2 * s.slen);
                for(j=s.slen-1; j>=0; j--){
                    self->countset[j] = counts[j];
                    self->totals[j] = counts[j] + total_index;
                    total_index >>= 2;
                }
            }else{
                counts = count_bipart_words(filename, s.slen,
                                            s.pos, s.glen);
                long total_index = 1ul << (2 * s.slen);
                int len = s.slen-s.pos;
                for(j=0; j<len; j++){
                    self->countset[s.hash-j] = counts[len-j-1];
                    self->totals[s.hash-j] = counts[len-j-1] + total_index;
                    total_index >>= 2;
                }
            }
            free(counts);
        }
        self->structs = structs;
        self->structs_num = structs_num;
    }
}

static void Counts_dealloc(Counts *self){
    SiteStruct *structs;
    structs = self->structs;
    int i;
    for(i=0; i<self->structs_num; i++){
        SiteStruct s = structs[i];
        long **ptr;
        ptr = self->countset + s.hash
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
    PyObject *dsite_list;
    int struct_hash = 1;
    if(!PyArg_ParseTuple(args, "O|i", &dsite_list, &struct_hash))
        return NULL;
    else{
        long *counts;
        counts = self->countset[struct_hash];
        long count = 0;
        int dsite_num = PyList_Size(dsite_list);
        long dsite;
        int i;
        for(i=0; i<dsite_num; i++){
            dsite = PyInt_AsLong(dsite_list[i])
            count += counts[dsite];
        }
        return Py_BuildValue("l", count);
    }
}

static PyObject *Counts_get_freq(Counts *self, PyObject *args,
                                 PyObject *kwds){
    PyObject *dsite_list;
    int struct_hash = 1;
    if(!PyArg_ParseTuple(args, "O|i", &dsite_list, &struct_hash))
        return NULL;
    else{
        long *counts;
        counts = self->countset[struct_hash];
        long count = 0;
        int dsite_num = PyList_Size(dsite_list);
        long dsite;
        int i;
        for(i=0; i<dsite_num; i++){
            dsite = PyInt_AsLong(dsite_list[i])
            count += counts[dsite];
        }
        double freq = (double) count;
        long total = *(self->totals[struct_hash])
        //if(total == 0)
        return Py_BuildValue("d", freq/total);
    }
}

static PyObject *Counts_get_total(Counts *self, PyObject *args,
                                  PyObject *kwds){
    int struct_hash = 1;
    if(!PyArg_ParseTuple(args, "|i", &struct_hash))
        return NULL;
    else{
        long total = *(self->totals[struct_hash])
        return Py_BuildValue("l", total);
    }
}

static PyTypeObject CountsType = {
    PyObject_HEAD_INIT(NULL)
    0,                         /*ob_size*/
    "counts.Counts",           /*tp_name*/
    sizeof(Counts),            /*tp_basicsize*/
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

static PyMethodDef counts_methods[] = {
    {NULL}
};

#ifndef PyMODINIT_FUNC
#define PyMODINIT_FUNC void
#endif

PyMODINIT_FUNC
initCounts(void){
    PyObject* counts;

    CountsType.tp_new = PyType_GenericNew;
    if (PyType_Ready(&CountsType) < 0)
        return;

    counts = Py_InitModule3("counts", counts_methods, "module help");

    Py_INCREF(&CountsType);
    PyModule_AddObject(counts, "Counts", (PyObject *)&CountsType);
}
