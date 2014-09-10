/** @file pygetkw.c
 * @brief C-interface to the Python getkw module
 *
 * Written by Jonas Juselius <jonas.juselius@chem.uit.no>
 * CTCC, University of Tromsø, April 2008
 *
 */

#define PYTHON_EXTMOD

#include <Python.h>
#include <getkw.h>
#include <assert.h>

#define ASSERT_VAL(val) if (val == NULL) { Py_DECREF(self->kw); \
	MSG_ERROR2_0("Value error:", path); return 0; }
#define ASSERT_MEM(mem) if (mem == NULL) { Py_DECREF(self->kw); \
		MSG_CRITICAL("malloc() failed!", MEM_ERROR); }

static int verbose=1;
static int strict=1;

/** Toggle the verbosity flag for explicit error messages
 *
 * @param v boolan
 * @return Original state of the verbosity flag
 */
int kw_set_verbose(int v)
{
	int old=verbose;
	verbose=v;
	return old;
}

/** Toggle the strict flag, making most errors fatal.
 *
 * @param v boolan
 * @return Original state of the strict flag
 */
int kw_set_strict(int v)
{
	int old=strict;
	strict=v;
	return old;
}

static int except_error(PyObject *obj)
{
	PyObject *a, *b, *c;

	if (obj != NULL) return 0;
	PyErr_Fetch(&a, &b, &c);

	if (verbose) {
		fprintf(stderr, "Error: ");
		if (b != NULL) PyObject_Print(b, stderr, Py_PRINT_RAW);
		fprintf(stderr, "\n");
	}

	Py_DECREF(a);
	Py_DECREF(b);
	Py_XDECREF(c);
	return 1;
}

static int check_kw_type(PyObject *kw, const char *type)
{
	PyObject *typ;
	int i;

	typ=PyObject_CallMethod(kw, "is_type", "s", type);
	if (except_error(typ) || typ == Py_None) {
		Py_DECREF(typ);
		return 0;
	}

	i=PyInt_AsLong(typ);
	Py_DECREF(typ);
	if (!i) return 0;
	return 1;
}

static int get_kw_obj(Getkw_t *self, const char *path, const char *type)
{
	PyObject *obj;

	obj=PyObject_CallMethod(self->getkw, "get_keyword", "s", path);
	if (except_error(obj)) return 0;
	if (obj == Py_None) {
		Py_DECREF(obj);
		return 0;
	}
	if (!check_kw_type(obj, type)) { 
		Py_DECREF(obj);
		if (verbose) MSG_ERROR2_0("Keyword is not of type", type);
		return 0;
	}

	self->kw=PyObject_GetAttrString(obj, "arg");
	if (except_error(self->kw)) {
		Py_DECREF(obj);
		return 0;
	}
	self->len=PyObject_Size(self->kw);
	Py_DECREF(obj);
	return 1;
}

/** Constructor for the Getkw "class"
 *
 * Allocates memory for the Getkw_t struct
 *
 * @param getkw Python Getkw object pointer
 * @return Getkw_t object
 */
Getkw_t *kw_InitGetkw(PyObject *getkw)
{
	Getkw_t *self;

	Py_INCREF(getkw);
	self=(Getkw_t *) malloc(sizeof(Getkw_t));
	if (self == NULL) return NULL;
	self->getkw=getkw;
	return self;
}

/** Destructor for the Getkw "class"
 *
 * @param self Python Getkw object pointer to be deallocated
 */
void kw_DelGetkw(Getkw_t *self)
{
	Py_DECREF(self->getkw);
	self->getkw=NULL;
	free(self);
}

/** Get 'long' keyword 
 *
 * @param self Getkw_t Object handle
 * @param path Path name of keyword to return
 * @param result Pointer to where to store result
 * @return 1 on succes, 0 on failure
 */
int kw_GetLong(Getkw_t *self, const char *path, long *result)
{	
	PyObject *val;

	if (!get_kw_obj(self, path, "INT")) 
		MSG_ERROR2_0("Invalid key:", path);
	if (self->len > 1) {
		Py_DECREF(self->kw);
		MSG_ERROR2_0("Keyword is not a scalar:", path);
	}

	val=PyList_GetItem(self->kw, 0);
	ASSERT_VAL(val);
	*result=PyInt_AsLong(val);

	Py_DECREF(self->kw);
	return 1;
}

int kw_GetLongArray(Getkw_t *self, const char *path, long **result)
{	
	PyObject *val;
	int i;

	if (!get_kw_obj(self, path, "INT_ARRAY")) 
		MSG_ERROR2_0("Invalid key:", path);

	*result=(long *) malloc(self->len*sizeof(long));
	ASSERT_MEM(*result);

	for (i=0; i < self->len; i++) {
		val=PyList_GetItem(self->kw, i);
		ASSERT_VAL(val);
		(*result)[i]=PyInt_AsLong(val);
	}

	Py_DECREF(self->kw);
	return self->len;
}

/** Get 'int' keyword 
 *
 * Note that Python stores ints as longs, the return value is thus cast to int
 *
 * @param self Getkw_t Object handle
 * @param path Path name of keyword to return
 * @param result Pointer to where to store result
 * @return 1 on succes, 0 on failure
 */
int kw_GetInt(Getkw_t *self, const char *path, int *result)
{	
	PyObject *val;

	if (!get_kw_obj(self, path, "INT")) 
		MSG_ERROR2_0("Invalid key:", path);
	if (self->len > 1) {
		Py_DECREF(self->kw);
		MSG_ERROR2_0("Keyword is not a scalar:", path);
	}

	val=PyList_GetItem(self->kw, 0);
	ASSERT_VAL(val);
	*result=(int) PyInt_AsLong(val);

	Py_DECREF(self->kw);
	return 1;
}

/** Get 'int array' keyword 
 *
 * Note: The returned array has been allocated and needs to be explicitly
 * deallocated by the user.
 *
 * @param self Getkw_t Object handle
 * @param path Path name of keyword to return
 * @param result Pointer to array where to store result
 * @return number of elements in array on succes, 0 on failure
 */
int kw_GetIntArray(Getkw_t *self, const char *path, int **result)
{	
	PyObject *val;
	int i;

	if (!get_kw_obj(self, path, "INT_ARRAY")) 
		MSG_ERROR2_0("Invalid key:", path);

	*result=(int *) malloc(self->len*sizeof(int));
	ASSERT_MEM(*result);

	for (i=0; i < self->len; i++) {
		val=PyList_GetItem(self->kw, i);
		ASSERT_VAL(val);
		(*result)[i]=(int) PyInt_AsLong(val);
	}

	Py_DECREF(self->kw);
	return self->len;
}

/** Get 'bool' keyword 
 *
 * @param self Getkw_t Object handle
 * @param path Path name of keyword to return
 * @param result Pointer to where to store result
 * @return 1 on succes, 0 on failure
 */
int kw_GetBool(Getkw_t *self, const char *path, int *result)
{	
	PyObject *val;

	if (!get_kw_obj(self, path, "BOOL")) 
		MSG_ERROR2_0("Invalid key:", path);
	if (self->len > 1) {
		Py_DECREF(self->kw);
		MSG_ERROR2_0("Keyword is not a scalar:", path);
	}

	val=PyList_GetItem(self->kw, 0);
	ASSERT_VAL(val);
	*result=(int) PyInt_AsLong(val);

	Py_DECREF(self->kw);
	return 1;
}

/** Get 'bool array' keyword 
 *
 * Note: The returned array has been allocated and needs to be explicitly
 * deallocated by the user.
 *
 * @param self Getkw_t Object handle
 * @param path Path name of keyword to return
 * @param result Pointer to array where to store result
 * @return number of elements in array on succes, 0 on failure
 */
int kw_GetBoolArray(Getkw_t *self, const char *path, int **result)
{	
	PyObject *val;
	int i;

	if (!get_kw_obj(self, path, "BOOL_ARRAY")) 
		MSG_ERROR2_0("Invalid key:", path);

	*result=(int *) malloc(self->len*sizeof(int));
	ASSERT_MEM(*result);

	for (i=0; i < self->len; i++) {
		val=PyList_GetItem(self->kw, i);
		ASSERT_VAL(val);
		(*result)[i]=(int) PyInt_AsLong(val);
	}

	Py_DECREF(self->kw);
	return self->len;
}

/** Get 'double' keyword 
 *
 * @param self Getkw_t Object handle
 * @param path Path name of keyword to return
 * @param result Pointer to where to store result
 * @return 1 on succes, 0 on failure
 */
int kw_GetDouble(Getkw_t *self, const char *path, double *result)
{	
	PyObject *val;

	if (!get_kw_obj(self, path, "DBL")) 
		MSG_ERROR2_0("Invalid key:", path);
	if (self->len > 1) {
		Py_DECREF(self->kw);
		MSG_ERROR2_0("Keyword is not a scalar:", path);
	}

	val=PyList_GetItem(self->kw, 0);
	ASSERT_VAL(val);
	*result=PyFloat_AsDouble(val);

	Py_DECREF(self->kw);
	return 1;
}

/** Get 'double array' keyword 
 *
 * Note: The returned array has been allocated and needs to be explicitly
 * deallocated by the user.
 *
 * @param self Getkw_t Object handle
 * @param path Path name of keyword to return
 * @param result Pointer to array where to store result
 * @return number of elements in array on succes, 0 on failure
 */
int kw_GetDoubleArray(Getkw_t *self, const char *path, double **result)
{	
	PyObject *val;
	int i;

	if (!get_kw_obj(self, path, "DBL_ARRAY")) 
		MSG_ERROR2_0("Invalid key:", path);

	*result=(double *) malloc(self->len*sizeof(double));
	ASSERT_MEM(*result);

	for (i=0; i < self->len; i++) {
		val=PyList_GetItem(self->kw, i);
		ASSERT_VAL(val);
		(*result)[i]=PyFloat_AsDouble(val);
	}

	Py_DECREF(self->kw);
	return self->len;
}

/** Get 'string' keyword 
 *
 * Note: The returned string has been allocated and needs to be explicitly
 * deallocated by the user.
 *
 * @param self Getkw_t Object handle
 * @param path Path name of keyword to return
 * @param result Pointer to where to store result
 * @return 1 on succes, 0 on failure
 */
int kw_GetString(Getkw_t *self, const char *path, char **result)
{	
	PyObject *val;
	int strsz;
	char *pystr;

	if (!get_kw_obj(self, path, "STR")) 
		MSG_ERROR2_0("Invalid key:", path);

	val=PyList_GetItem(self->kw, 0);
	if (except_error(val)) return 0;
	strsz=PyString_Size(val)+1;
	*result=(char *) malloc(strsz*sizeof(char));
	ASSERT_MEM(*result);
	pystr=PyString_AsString(val);
	strncpy(*result, pystr, strsz);

	Py_DECREF(self->kw);
	return self->len;
}

/** Get 'string array' keyword 
 *
 * Note: The returned array, and each element of the array have been 
 * allocated and needs to be explicitly deallocated by the user. 
 *
 * @param self Getkw_t Object handle
 * @param path Path name of keyword to return
 * @param result Pointer to array of pointers where to store result
 * @return number of elements in array on succes, 0 on failure
 */
int kw_GetStringArray(Getkw_t *self, const char *path, char ***result)
{	
	PyObject *val;
	int i, strsz;
	char *pystr;

	if (!get_kw_obj(self, path, "STR_ARRAY")) 
		MSG_ERROR2_0("Invalid key:", path);

	*result=(char **) malloc(self->len*sizeof(char *));
	ASSERT_MEM(*result);
	for (i=0; i < self->len; i++) {
		val=PyList_GetItem(self->kw, i);
		ASSERT_VAL(val);
		strsz=PyString_Size(val)+1;
		(*result)[i]=(char *) malloc(strsz*sizeof(char));
		ASSERT_MEM((*result)[i]);
		pystr=PyString_AsString(val);
		strncpy((*result)[i], pystr, strsz);
	}
	Py_DECREF(self->kw);
	return self->len;
}

/** Get 'data' keyword 
 *
 * This function is essentially identical to kw_GetStringArray, but checks for
 * keywords of type 'DATA'.
 * Note: The returned array, and each element of the array have been 
 * allocated and needs to be explicitly deallocated by the user. 
 *
 * @param self Getkw_t Object handle
 * @param path Path name of keyword to return
 * @param result Pointer to array of pointers where to store result
 * @return number of elements in array on succes, 0 on failure
 */
int kw_GetData(Getkw_t *self, const char *path, char ***result)
{
	PyObject *val;
	int i, strsz;
	char *pystr;

	if (!get_kw_obj(self, path, "DATA")) 
		MSG_ERROR2_0("Invalid key:", path);

	*result=(char **) malloc(self->len*sizeof(char *));
	ASSERT_MEM(*result);
	for (i=0; i < self->len; i++) {
		val=PyList_GetItem(self->kw, i);
		ASSERT_VAL(val);
		strsz=PyString_Size(val)+1;
		(*result)[i]=(char *) malloc(strsz*sizeof(char));
		ASSERT_MEM((*result)[i]);
		pystr=PyString_AsString(val);
		strncpy((*result)[i], pystr, strsz);
	}
	Py_DECREF(self->kw);
	return self->len;
}

/** Temporarily push a path on to the section stack
 *
 * Useful for implementing objects which can be multiply defined 
 * (e.g. grids, coords, basis sets, etc.), but in different sections.
 *
 * @param self Getkw_t Object handle
 * @param path Path name of section to push
 * @return 1 on succes, 0 on failure
 */
int kw_PushSection(Getkw_t *self, const char *path)
{
	PyObject *tmp;

	tmp=PyObject_CallMethod(self->getkw, "push_sect", "s", path);
	if (tmp == NULL || tmp == Py_None) {
		Py_DECREF(tmp);
		if (verbose) MSG_ERROR2_0("Push failed:", path);
		return 0;
	}
	
	return 1;
}

/** Pop a path from the section stack
 *
 * @param self Getkw_t Object handle
 * @return 1 on succes, 0 on failure
 */
int kw_PopSection(Getkw_t *self)
{
	PyObject *tmp;
	tmp=PyObject_CallMethod(self->getkw, "pop_sect", NULL);
	if (tmp == NULL || tmp == Py_None) {
		Py_DECREF(tmp);
		if (verbose) MSG_ERROR_0("Pop failed!");
		return 0;
	}

	return 1;
}

/** Test if a section has a particular keyword
 *
 * @param self Getkw_t Object handle
 * @param path Path name of keyword to find
 * @return 1 on succes, 0 on failure
 */
int kw_HasKeyword(Getkw_t *self, const char *path)
{
	PyObject *obj;

	obj=PyObject_CallMethod(self->getkw, "getkw", "s", path);
	if (obj == NULL) return 0;
	if (obj == Py_None) {
		Py_DECREF(obj);
		return 0;
	}
	Py_DECREF(obj);
	return 1;
}

/** Test if a section has a particular sub-section
 *
 * @param self Getkw_t Object handle
 * @param path Path name of section to find 
 * @return 1 on succes, 0 on failure
 */
int kw_HasSection(Getkw_t *self, const char *path)
{
	PyObject *obj;

	obj=PyObject_CallMethod(self->getkw, "find_sect", "s", path);
	if (obj == NULL) return 0;
	if (obj == Py_None) {
		Py_DECREF(obj);
		return 0;
	}
	Py_DECREF(obj);
	return 1;
}

/** Test if a section has been explicitly set in the input, or whether it will
 * return the default values.
 *
 * @param self Getkw_t Object handle
 * @param path Path name of section to find 
 * @return 1 on succes, 0 on failure
 */
int kw_SectionIsSet(Getkw_t *self, const char *path)
{
	PyObject *obj, *result;
	long i;

	obj=PyObject_CallMethod(self->getkw, "find_sect", "s", path);
	if (except_error(obj)) return -1;
	result=PyObject_CallMethod(obj, "is_set", NULL);
	if (except_error(result)) {
		Py_DECREF(obj);
		return -1;
	}
	i=PyInt_AsLong(result);

	Py_DECREF(obj);
	Py_DECREF(result);
	return i;
}

/** Test if a keyword has been explicitly set in the input, or whether it will
 * return the default value.
 *
 * @param self Getkw_t Object handle
 * @param path Path name of keyword to find 
 * @return 1 on succes, 0 on failure
 */
int kw_KeywordIsSet(Getkw_t *self, const char *path)
{
	PyObject *obj, *result;
	long i;

	obj=PyObject_CallMethod(self->getkw, "get_keyword", "s", path);
	if (except_error(obj)) return -1;
	result=PyObject_CallMethod(obj, "is_set", NULL);
	if (except_error(result)) {
		Py_DECREF(obj);
		return -1;
	}
	i=PyInt_AsLong(result);

	Py_DECREF(obj);
	Py_DECREF(result);
	return i;
}


