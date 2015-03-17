/** @file cgetkw.c
 * @brief C-interface to the Python getkw module
 *
 * Written by Jonas Juselius <jonas.juselius@chem.uit.no>
 * CTCC, University of Tromsø, April 2008
 *
 */

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <assert.h>
#include "getkw.h"

#define MAX_PATH 256
#define MAX_BUF 1001
#define MAX_KEYLEN 33

#define ASSERT_KEY(x,p) if (x == NULL ) { if (strict) { \
	printf("Error: Key not found: %s\n", p); abort(); } else { return 0; } }
#define ASSERT_TYPE(a,b,path) if (a->type != b) { \
	printf("Error: Invalid type for %s, expected %d got %d\n", path, b,\
	a->type); if (strict) abort(); return 0; }
#define ASSERT_MEM(X) {if (X == NULL) { printf("Critical: %s(), line %d: %s\n", __func__, __LINE__, "Memory allocation failed!"); abort();} }


static int verbose=1;
static int strict=1;

static Keyword_t *new_keyword(const char *name, int type, int len)
{
	Keyword_t *self;
	int strln;

	self=(Keyword_t *) calloc(1,sizeof(Keyword_t));
	ASSERT_MEM(self);
	strln=strlen(name);
	self->name=(char *) calloc(strln+1,sizeof(char));
	ASSERT_MEM(self->name);
	strncpy(self->name, name, strln);
	self->type=type;
	self->len=len;

	switch(type) {
		case INT:
		case DBL:
		case BOOL:
		case STR:
			if (len != 1) {
				MSG_WARN("not an array");
			}
			break;
		case INT_ARRAY:
			self->val.iarray=(int *) calloc(len, sizeof(int));
			ASSERT_MEM(self->val.iarray);
			break;
		case DBL_ARRAY:
			self->val.darray=(double *) calloc(len, sizeof(double));
			ASSERT_MEM(self->val.darray);
			break;
		case BOOL_ARRAY:
			self->val.larray=(int *) calloc(len, sizeof(int));
			ASSERT_MEM(self->val.larray);
			break;
		case STR_ARRAY:
		case DATA:
			self->val.data=(char **) calloc(len, sizeof(char *));
			ASSERT_MEM(self->val.data);
			break;
		default:
			MSG_ERROR_0("Invalid type");
			break;
	}
	return self;
}

static void del_keyword(Keyword_t *self)
{
	int i;

	switch(self->type) {
		case INT:
		case DBL:
		case BOOL:
			break;
		case STR:
			free(self->val.str);
			break;
		case INT_ARRAY:
			free(self->val.iarray);
			break;
		case DBL_ARRAY:
			free(self->val.darray);
			break;
		case BOOL_ARRAY:
			free(self->val.larray);
			break;
		case STR_ARRAY:
		case DATA:
			for (i=0; i < self->len; i++) {
				free(self->val.data[i]);
			}
			free(self->val.data);
			break;
		default:
			MSG_ERROR("Invalid type");
			break;
	}
	free(self->name);
	free(self);
}

/* Very simple implementation, since we know everything we need beforehand,
 * no need for elaborate linked lists or anything like that. */
static Section_t *new_section(const char *name,
		int nkeys, int nsect)
{
	Section_t *self;

	self=(Section_t *) calloc(1,sizeof(Section_t));
	ASSERT_MEM(self);

	int strln=strlen(name);
	self->name=(char *) calloc(strln+1,sizeof(char));
	ASSERT_MEM(self->name);
	strncpy(self->name, name, strln);

	if (nkeys > 0) {
		self->kw=(Keyword_t **) calloc(nkeys,sizeof(Keyword_t *));
		ASSERT_MEM(self->kw);
		self->nkeys=nkeys;
	}
	if (nsect > 0) {
		self->sect=(Section_t **) calloc(nsect,sizeof(Section_t *));
		ASSERT_MEM(self->sect);
		self->nsect=nsect;
	}
	return self;
}

/* recursively delete current section and all its subsections */
static void del_section(Section_t *self)
{
	int i;

	for (i=0; i < self->nkeys; i++) {
		del_keyword(self->kw[i]);
	}
	for (i=0; i < self->nsect; i++) {
		del_section(self->sect[i]);
	}
	free(self->name);
	free(self->kw);
	free(self->sect);
	free(self);
}

static Keyword_t *getkw(Section_t *self, const char *name)
{
	int i, len;

	for (i=0; i < self->nkeys; i++) {
		len=strlen(self->kw[i]->name);
		if (strncmp(self->kw[i]->name, name, len) == 0 ) {
			return self->kw[i];
		}
	}
	return NULL;
}

static Section_t *getsect(Section_t *self, const char *name)
{
	int i, len;

	for (i=0; i < self->nsect; i++) {
/*        printf("kuk %d %p %s\n", i, self->sect[i], self->sect[i]->name);*/
/*        printf("name=%s\n", name);*/
		len=strlen(self->sect[i]->name);
		if (strncmp(self->sect[i]->name, name, len) == 0 ) {
/*            printf("foundit!\n");*/
			return self->sect[i];
		}
	}
	return NULL;
}

static Section_t *findsect(Section_t *self, const char *path)
{
	const char *name;
	char tmp[MAX_PATH];
	Section_t *sect;
	int len, i;

	name=index(path, '.');
	if (name == NULL) {
		if (strcmp(path, self->name) == 0) {
			return self;
		}
		for (i=0; i < self->nsect; i++) {
			if (strcmp(path, self->sect[i]->name) == 0) {
				return self->sect[i];
			}
		}
		return self;
	}

	len=(name-path)/sizeof(char);
	if (len > MAX_PATH-1) {
		MSG_ERROR("Too long path name");
		return NULL;
	}
	strncpy(tmp, path, len);
	tmp[len]='\0';

	sect=getsect(self, tmp);
	if (sect == NULL) {
		MSG_ERROR("Invalid section");
		return NULL;
	}

	return findsect(sect, name+1);
}

static Keyword_t *findkw(Section_t *self, const char *path)
{
	int i;
	const char *name;
	Section_t *sect;

	name=rindex(path, '.');
	if (name == NULL) {
		name=path;
		sect=self;
	} else {
		name++; /* get rid of the dot */
		sect=findsect(self, path);
		if (sect == NULL) {
			return NULL;
		}
	}
	for (i=0; i < sect->nkeys; i++) {
		if (strcmp(sect->kw[i]->name, name) == 0 ) {
			return sect->kw[i];
		}
	}

	return NULL;
}

static int conv_type(const char *type)
{
	if (strncmp(type, "INT_ARRAY", 9) == 0) {
		return INT_ARRAY;
	}
	if (strncmp(type, "DBL_ARRAY", 9) == 0) {
		return DBL_ARRAY;
	}
	if (strncmp(type, "BOOL_ARRAY", 10) == 0) {
		return BOOL_ARRAY;
	}
	if (strncmp(type, "STR_ARRAY", 9) == 0) {
		return DATA;
	}
	if (strncmp(type, "DATA", 4) == 0) {
		return DATA;
	}
	if (strncmp(type, "INT", 3) == 0) {
		return INT;
	}
	if (strncmp(type, "DBL", 3) == 0) {
		return DBL;
	}
	if (strncmp(type, "BOOL", 4) == 0) {
		return BOOL;
	}
	if (strncmp(type, "STR", 3) == 0) {
		return STR;
	}
	return -1;
}

static int conv_bool(char p)
{
	if (p == 'F') return 0;
	return 1;
}

static char *scan_line(FILE *fd)
{
	char buf[MAX_BUF];
	char *line;
	int len;

	len=fscanf(fd, "%[^\n]\n", buf);
	len=strlen(buf);
	line=(char *) calloc(len+1,sizeof(char));
	strncpy(line, buf, len);
	return line;
}


static Keyword_t *read_key(FILE *fd)
{
	Keyword_t *kw;
	char name[MAX_KEYLEN];
	char type[11];
	int len, i, typ, n;
	char set[6];

	n=fscanf(fd, "%s %s %d %s\n", type, name, &len, set);
/*    printf("KEY: %s %s %d %c\n", type, name, len, set[0]);*/

	typ=conv_type(type);
	kw=new_keyword(name, typ, len);

	kw->set=conv_bool(set[0]);

	for (i=0; i < kw->len; i++) {
		switch(kw->type) {
			case INT:
				n=fscanf(fd, "%d\n", &kw->val.ival);
				break;
			case DBL:
				n=fscanf(fd, "%lf\n", &kw->val.dval);
				break;
			case BOOL:
				n=fscanf(fd, "%s\n", set);
				kw->val.lval=conv_bool(set[0]);
				break;
			case STR:
				kw->val.str=scan_line(fd);
				break;
			case INT_ARRAY:
				n=fscanf(fd, "%d\n", &kw->val.iarray[i]);
				break;
			case DBL_ARRAY:
				n=fscanf(fd, "%lf\n", &kw->val.darray[i]);
				break;
			case BOOL_ARRAY:
				n=fscanf(fd, "%d\n", &kw->val.larray[i]);
				break;
			case STR_ARRAY:
			case DATA:
				kw->val.data[i]=scan_line(fd);
				break;
			default:
				MSG_ERROR("Invalid type");
				return NULL;
				break;
		}
	}
	return kw;
}

static Section_t *read_sect(FILE *fd)
{
	Section_t *sect;
	char name[MAX_KEYLEN];
	char tagname[2*MAX_KEYLEN+2];
	char tag[MAX_KEYLEN];
	int nsect, nkeys, i;
	char set[6];

	i=fscanf(fd, "%*s %s %d %s\n", name, &nsect, set);
	i=fscanf(fd, "%*s %1c %*s %d\n", tag, &nkeys);
/*    printf("SECT: %s %c %d %d %c \n", name, tag[0], nsect, nkeys, set[0]);*/
	if (tag[0] == 'T') {
		i=fscanf(fd, "%s\n", tag);
		sprintf(tagname, "%s(%s)", name, tag);
/*        printf("TAG: %s -> %s\n", tag, tagname);*/
		sect=new_section(tagname, nkeys, nsect);
	} else {
		sect=new_section(name, nkeys, nsect);
	}

	sect->set=conv_bool(set[0]);

	for (i=0; i < nkeys; i++) {
		sect->kw[i]=read_key(fd);
	}
	for (i=0; i < nsect; i++) {
		sect->sect[i]=read_sect(fd);
	}

	return sect;
}

static Section_t *read_input(FILE *fd)
{
	Section_t *top;
	top=read_sect(fd);

	return top;
}

static void print_key(Keyword_t *kw)
{
	int i;

	printf("TYPE(%d) %s %d %d\n", kw->type, kw->name, kw->len, kw->set);
	for (i=0; i < kw->len; i++) {
		switch(kw->type) {
			case INT:
				printf("%d\n", kw->val.ival);
				break;
			case DBL:
				printf("%lf\n", kw->val.dval);
				break;
			case BOOL:
				printf("%d\n", kw->val.lval);
				break;
			case STR:
				printf("%s\n", kw->val.str);
				break;
			case INT_ARRAY:
				printf("%d\n", kw->val.iarray[i]);
				break;
			case DBL_ARRAY:
				printf("%lf\n", kw->val.darray[i]);
				break;
			case BOOL_ARRAY:
				printf("%d\n", kw->val.larray[i]);
				break;
			case STR_ARRAY:
			case DATA:
				printf("%s\n", kw->val.data[i]);
				break;
			default:
				MSG_ERROR("Invalid type");
				break;
		}
	}
}

void kw_PrintSection(Section_t *self)
{
	int i;

	printf("SECT: %s sects=%d keys=%d set=%d\n", self->name, self->nsect,
			self->nkeys, self->set);
	for (i=0; i < self->nkeys; i++)
	{
		print_key(self->kw[i]);
	}
	for (i=0; i < self->nsect; i++)
	{
		kw_PrintSection(self->sect[i]);
	}
}


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


/** Constructor for the Getkw "class"
 *
 * Allocates memory for the Getkw_t struct
 *
 * @param getkw Python Getkw object pointer
 * @return Getkw_t object
 */
Getkw_t *kw_InitGetkw(const char *file)
{
	Getkw_t *self;
	FILE *fd;

	self=(Getkw_t *) calloc(1,sizeof(Getkw_t));
	ASSERT_MEM(self);

	if ( file == NULL || strncasecmp(file, "stdin", 5) == 0) {
		fd=stdin;
	} else {
		fd=fopen(file, "r");
	}

	if ( fd == NULL ) {
		MSG_CRITICAL("fopen() failed", 256);
	}

	self->toplevel=read_input(fd);
	self->cur=self->toplevel;
	self->sp=0;

	if (fd != stdin) {
		fclose(fd);
	}

	return self;
}

/** Destructor for the Getkw "class"
 *
 * @param self Python Getkw object pointer to be deallocated
 */
void kw_DelGetkw(Getkw_t *self)
{
	del_section(self->toplevel);
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
	Keyword_t *key;

	key=findkw(self->cur, path);
	ASSERT_KEY(key,path);
	ASSERT_TYPE(key, INT,path);
	*result=(long) key->val.ival;
	return key->len;
}

int kw_GetLongArray(Getkw_t *self, const char *path, long **result)
{
	Keyword_t *key;

	key=findkw(self->cur, path);
	ASSERT_KEY(key,path);
	ASSERT_TYPE(key, INT_ARRAY,path);
	*result=(long *) key->val.iarray;
	return key->len;
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
	Keyword_t *key;

	key=findkw(self->cur, path);
	ASSERT_KEY(key,path);
	ASSERT_TYPE(key, INT,path);
	*result=key->val.ival;
	return key->len;
}

int kw_GetIntRef(Getkw_t *self, const char *path, int **result)
{
	Keyword_t *key;

	key=findkw(self->cur, path);
	ASSERT_KEY(key,path);
	ASSERT_TYPE(key, INT,path);
	*result=&key->val.ival;
	return key->len;
}

/** Get 'int array' keyword
 *
 * @param self Getkw_t Object handle
 * @param path Path name of keyword to return
 * @param result Pointer to array where to store result
 * @return number of elements in array on succes, 0 on failure
 */
int kw_GetIntArray(Getkw_t *self, const char *path, int **result)
{
	Keyword_t *key;

	key=findkw(self->cur, path);
	ASSERT_KEY(key,path);
	ASSERT_TYPE(key, INT_ARRAY,path);
	*result= key->val.iarray;
	return key->len;
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
	Keyword_t *key;

	key=findkw(self->cur, path);
	ASSERT_KEY(key,path);
	ASSERT_TYPE(key, BOOL,path);
	*result=key->val.ival;
	return key->len;
}

/** Get 'bool array' keyword
 *
 * @param self Getkw_t Object handle
 * @param path Path name of keyword to return
 * @param result Pointer to array where to store result
 * @return number of elements in array on succes, 0 on failure
 */
int kw_GetBoolArray(Getkw_t *self, const char *path, int **result)
{
	Keyword_t *key;

	key=findkw(self->cur, path);
	ASSERT_KEY(key,path);
	ASSERT_TYPE(key, BOOL_ARRAY,path);
	*result=key->val.iarray;
	return key->len;
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
	Keyword_t *key;

	key=findkw(self->cur, path);
	ASSERT_KEY(key,path);
	ASSERT_TYPE(key, DBL,path);
	*result=key->val.dval;
	return key->len;
}

int kw_GetDoubleRef(Getkw_t *self, const char *path, double **result)
{
	Keyword_t *key;

	key=findkw(self->cur, path);
	ASSERT_KEY(key,path);
	ASSERT_TYPE(key, DBL,path);
	*result=&key->val.dval;
	return key->len;
}

/** Get 'double array' keyword
 *
 * @param self Getkw_t Object handle
 * @param path Path name of keyword to return
 * @param result Pointer to array where to store result
 * @return number of elements in array on succes, 0 on failure
 */
int kw_GetDoubleArray(Getkw_t *self, const char *path, double **result)
{
	Keyword_t *key;

	key=findkw(self->cur, path);
	ASSERT_KEY(key,path);
	ASSERT_TYPE(key, DBL_ARRAY,path);
	*result=key->val.darray;
	return key->len;
}

/** Get 'string' keyword
 *
 * @param self Getkw_t Object handle
 * @param path Path name of keyword to return
 * @param result Pointer to where to store result
 * @return 1 on succes, 0 on failure
 */
int kw_GetString(Getkw_t *self, const char *path, char **result)
{
	Keyword_t *key;

	key=findkw(self->cur, path);
	ASSERT_KEY(key,path);
	ASSERT_TYPE(key, STR,path);
	*result=key->val.str;
	return key->len;
}

/** Get 'string array' keyword
 *
 * @param self Getkw_t Object handle
 * @param path Path name of keyword to return
 * @param result Pointer to array of pointers where to store result
 * @return number of elements in array on succes, 0 on failure
 */
int kw_GetStringArray(Getkw_t *self, const char *path, char ***result)
{
	Keyword_t *key;

	key=findkw(self->cur, path);
	ASSERT_KEY(key,path);
	ASSERT_TYPE(key, STR_ARRAY,path);
	*result=key->val.data;
	return key->len;
}

/** Get 'data' keyword
 *
 * This function is essentially identical to kw_GetStringArray, but checks for
 * keywords of type 'DATA'.
 *
 * @param self Getkw_t Object handle
 * @param path Path name of keyword to return
 * @param result Pointer to array of pointers where to store result
 * @return number of elements in array on succes, 0 on failure
 */
int kw_GetData(Getkw_t *self, const char *path, char ***result)
{
	Keyword_t *key;

	key=findkw(self->cur, path);
	ASSERT_KEY(key,path);
	ASSERT_TYPE(key, DATA,path);
	*result=key->val.data;
	return key->len;
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
	Section_t *sect;

	sect=findsect(self->cur, path);
	if (sect == NULL) {
		MSG_ERROR("Section lookup failed!");
		return 0;
	}

	if (self->sp < MAX_SECT_STACK) {
		self->stack[self->sp]=self->cur;
		self->cur=sect;
		self->sp++;
		return 1;
	}

	MSG_ERROR_0("Stack overflow!");
}

/** Pop a path from the section stack
 *
 * @param self Getkw_t Object handle
 * @return 1 on succes, 0 on failure
 */
int kw_PopSection(Getkw_t *self)
{
	if (self->sp < 0) {
		self->cur=self->toplevel;
		MSG_ERROR_0("Stack underflow!");
	}
	self->sp--;
	self->cur=self->stack[self->sp];
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
	Keyword_t *key;

	key=findkw(self->cur, path);
	if (key == NULL) {
		return 0;
	}
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
	Section_t *sect;

	sect=findsect(self->cur, path);
	if (sect == NULL) {
		return 0;
	}
	if (sect == self->cur) {
		if (strcmp(path, sect->name) == 0) {
			return 1;
		}
		return 0;
	}
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
	Section_t *sect;

	sect=findsect(self->cur, path);
	if (sect == NULL) {
		return 0;
	}
	return sect->set;
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
	Keyword_t *key;

	key=findkw(self->cur, path);
	if (key == NULL) {
		return 0;
	}
	return key->set;
}


#ifdef TEST

static void write_test_inp(void)
{
	char *foo[]={"SECT toplevel 2 False",
"TAG F KW 2",
"INT apa 1 True",
"1",
"DBL rel_prec 1 True",
"0.0001",
"SECT raboof 1 True",
"TAG F KW 1",
"DATA COORD 3 True",
"1 1.0",
"2 2.0",
"3 3.0",
"SECT foobar 0 True",
"TAG T KW 2",
"gnat",
"INT foo 1 True",
"42",
"INT bar 1 True",
"84",
"SECT defs 0 True",
"TAG T KW 2",
"apa",
"INT_ARRAY foo 4 True",
"-42",
"-42",
"-42",
"-42",
"INT bar 1 True",
"-84",
NULL };
	int i=0;
	FILE *fd=fopen("test.inp", "w");
	while (foo[i] != NULL) {
		fprintf(fd, "%s\n", foo[i]);
		i++;
	}
	fclose(fd);
	printf("Wrote 'test.inp'\n");
}

int main(void)
{
	Getkw_t *input;
	int i, n;
	int *a;
	char **data;

	write_test_inp();
	input=kw_InitGetkw("test.inp");
	printf("*******************************************************\n");
	printf("Input %p\n", input);
	printf("*******************************************************\n");
	kw_PrintSection(input->toplevel);
	printf("*******************************************************\n");

	n=kw_GetInt(input, "raboof.foobar(gnat).foo", &i);
	printf("raboof.foobar(gnat).foo %d, len=%d\n", i, n);

	kw_PushSection(input, "raboof");
	n=-1;
	i=-1;
	n=kw_GetInt(input, "foobar(gnat).bar", &i);
	printf("foobar(gnat).bar %d, len=%d\n", i, n);
	kw_PopSection(input);

	n=kw_GetIntArray(input, "defs(apa).foo", &a);
	printf("defs(apa): %d\n", n);
	for (i=0; i < n; i++) {
		printf("%d ", a[i]);
	}
	printf("\n");

	n=kw_GetData(input, "raboof.COORD", &data);
	printf("COORD: %d\n", n);
	for (i=0; i < n; i++) {
		printf("%s\n", data[i]);
	}
	printf("\n");
	printf("foot(): %d\n", kw_HasSection(input, "foot"));
	printf("raboof(): %d\n", kw_HasSection(input, "raboof"));
	printf("raboof() set: %d\n", kw_SectionIsSet(input, "raboof"));
	printf("has key: %d\n", kw_HasKeyword(input, "rabd"));

	kw_DelGetkw(input);
	unlink("test.inp");
	printf("done.\n");

	return 1;
}
#endif
