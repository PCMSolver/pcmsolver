#ifndef GETKW_H
#define GETKW_H

#define MEM_ERROR 256
#define MAX_SECT_STACK 10

#define MSG_DEBUG(X) {printf("Debug: %s(), line %d in %s: %s\n", __func__, __LINE__, __FILE__, X);}
#define MSG_INFO(X) {printf("Info: %s(): %s\n", __func__, X);}
#define MSG_WARN(X) {printf("Warning: %s(): %s\n", __func__, X);}
#define MSG_ERROR(X) {printf("Error: %s(): %s\n", __func__, X);}
#define MSG_ERROR_0(X) {printf("Error: %s(): %s\n", __func__, X); return 0;}
#define MSG_ERROR2_0(X,A) {printf("Error: %s(): %s: %s\n", __func__, X,A); return 0;}
#define MSG_CRITICAL(X,Y) {printf("Critical: %s(), line %d: %s\n", __func__, __LINE__, X); exit(Y);}
#define ASSERT_MEM(X) {if (X == NULL) { printf("Critical: %s(), line %d: %s\n", __func__, __LINE__, "Memory allocation failed!"); abort();} }

enum kw_type {
	INT,
	DBL,
	STR,
	BOOL,
	INT_ARRAY,
	DBL_ARRAY,
	STR_ARRAY,
	BOOL_ARRAY,
	DATA
};

typedef struct kw_keyword {
	char *name;
	int type;
	int len;
	int set;
	union {
		int ival;
		int lval;
		double dval;
		char *str;
		int *iarray;
		int *larray;
		double *darray;
		char **data;
	} val;
} Keyword_t;

typedef struct kw_section {
	char *name;
	int set;
	int nkeys;
	int nsect;
	struct kw_section **sect;
	struct kw_keyword **kw;
} Section_t;

typedef struct {
	Section_t *toplevel;
	Section_t *cur;
	Section_t *stack[MAX_SECT_STACK];
	int sp;
} Getkw_t;

Getkw_t *kw_InitGetkw(char *file);
void kw_DelGetkw(Getkw_t *self);

int kw_GetLong(Getkw_t *inp, const char *path, long *result);
int kw_GetLongArray(Getkw_t *inp, const char *path, long **result);
int kw_GetInt(Getkw_t *inp, const char *path, int *result);
int kw_GetIntRef(Getkw_t *inp, const char *path, int **result);
int kw_GetIntArray(Getkw_t *inp, const char *path, int **result);
int kw_GetDouble(Getkw_t *inp, const char *path, double *result);
int kw_GetDoubleRef(Getkw_t *inp, const char *path, double **result);
int kw_GetDoubleArray(Getkw_t *inp, const char *path, double **result);
int kw_GetBool(Getkw_t *inp, const char *path, int *result);
int kw_GetBoolArray(Getkw_t *inp, const char *path, int **result);
int kw_GetString(Getkw_t *self, const char *path, char **result);
int kw_GetStringArray(Getkw_t *self, const char *path, char ***result);
int kw_GetData(Getkw_t *self, const char *path, char ***result);

int kw_PushSection(Getkw_t *self, const char *path);
int kw_PopSection(Getkw_t *self);
int kw_HasKeyword(Getkw_t *self, const char *path);
int kw_HasSection(Getkw_t *self, const char *path);
int kw_KeywordIsSet(Getkw_t *self, const char *path);
int kw_SectionIsSet(Getkw_t *self, const char *path);

int kw_set_verbose(int v);
int kw_set_strict(int s);
void kw_PrintSection(Section_t *self);

#endif
