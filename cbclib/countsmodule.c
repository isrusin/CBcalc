#include <Python.h>

#include "countscalc.h"

const int MAX_LENGTH = 14;
const int MAX_GAP_LENGTH = 14;

typedef struct{
    PyObject_HEAD
    long **countset;
} Counts;

static void Counts_dealloc(Counts *self){
    //free counts
    self->ob_type->tp_free((PyObject *)self);
}

static PyObject *Counts_new(PyTypeObject *type, PyObject *args,
                            PyObject *kwds){
    Counts *self;
    self = (Counts *)type->tp_alloc(type, 0);
    return (PyObject *)self;
}

typedef struct{
    int slen, pos, glen;
    int hash;
} Struct;

static void fill_struct_hash(Struct *s){
    int struct_hash = 0;
    struct_hash += s->slen;
    if(s->glen){
        struct_hash += MAX_LENGTH + 1;
        struct_hash += (s->pos - 1) * MAX_LENGTH;
        struct_hash += (s->glen - 1) * (MAX_LENGTH - 1) * MAX_LENGTH;
    }
    s->hash = struct_hash;
}

static int Counts_init(Counts *self, PyObject *args, PyObject *kwds){
    char *filename;
    PyObject *structs_tuple;
    if(!PyArg_ParseTuple(args, "sO", &filename, &structs_tuple))
        return NULL;
    else{
        int structs_num = PyTuple_Size(structs_tuple);
        int i;
        Struct *structs;
        structs = (Struct *)calloc(structs_num, sizeof(Struct));
        int max_hash = 0;
        for(i=0; i<structs_num; i++){
            PyObject *item = PyTuple_GetItem(structs_tuple, i);
            Struct s = structs[i];
            if(!PyArg_ParseTuple(item, "iii", &s.slen, &s.pos, &s.glen))
                continue;
            fill_struct_hash(&s);
            if(s.hash > max_hash)
                max_hash = s.hash;
        }
        self->countset = (long **)calloc(max_hash+1, sizeof(long *));
        for(i=0; i<structs_num; i++){
            Struct s = structs[i];
            long **counts;
            long **pcountset;
            pcountset = self->countset;
            int j;
            if(s.glen)
                counts = count_short_words(filename, s.slen);
                for(j=0; j<s.slen; j++)
                    pcountset[j] = counts[j]
            else
                counts = count_bipart_words(filename, s.slen,
                                            s.pos, s.glen);
                int len = s.slen-s.pos;
                for(j=0; j<len; j++)
                    pcountset[s.hash-j] = counts[len-j-1];
            free(counts);
        }
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
