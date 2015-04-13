#include <Python.h>
#include <getkw.h>

extern Getkw_t *input;

int getkw_tester(void)
{
	int i=0, n=0;
	int *ia;
	double *da;
	double d;
	char *str;
	char **str_a;
	
	i=kw_GetInt(input, "K_order", &n);
	printf("K_order: %d (%d)\n", n,i);
	
	i=kw_GetInt(input, "I_or_L", &n);
	printf("I_or_L: %d (%d)\n", n,i);
	
	i=kw_GetInt(input, "uniform", &n);
	printf("uniform: %d (%d)\n", n,i);
	
	i=kw_GetInt(input, "finest_scale", &n);
	printf("finest_scale: %d (%d)\n", n,i);
	
	i=kw_GetInt(input, "f_thr_type", &n);
	printf("f_thr_type: %d (%d)\n", n,i);

	i=kw_GetDouble(input, "rel_prec", &d);
	printf("rel_prec: %f (%d)\n", d,i);

	i=kw_GetDouble(input, "a_f", &d);
	printf("a_f: %f (%d)\n", d,i);

	i=kw_GetDouble(input, "shift_f", &d);
	printf("shift_f: %f (%d)\n", d,i);

	i=kw_PushSection(input, "func(inp)");
	n=0;
	i=kw_GetInt(input, "dim", &n);
	printf("push: func(inp), dim=%d (%d)\n", n, i);
	i=kw_PopSection(input);

/*
	i=kw_GetIntArray(input, "int_array", &ia);
	printf("int: [%d %d %d] (%d)\n", ia[0], ia[1], ia[2], i);
	free(ia);

	i=kw_GetDouble(input, "dbl", &d);
	printf("dbl: %f (%d)\n", d,i);

	i=kw_GetDoubleArray(input, "dbl_array", &da);
	printf("dbl_array: [%f %f %f] (%d)\n", da[0], da[1], da[2], i);
	free(da);

	i=kw_GetString(input, "str", &str);
	printf("str: %s (%d)\n", str, i);
	free(str);

	i=kw_GetStringArray(input, "str_array", &str_a);
	printf("str: [%s, %s, %s] (%d)\n", str_a[0], str_a[1], str_a[2], i);
	free(str_a[0]);
	free(str_a[1]);
	free(str_a[2]);
	free(str_a);

	i=kw_GetBool(input, "bool", &n);
	printf("bool: %d (%d)\n", n,i);

	i=kw_KeywordIsSet(input, "basis.atoms");
	printf("kw is set: %d basis.atoms\n", i);

	i=kw_KeywordIsSet(input, "sect1.foo");
	printf("kw is set: %d\n", i);

	i=kw_GetInt(input, "sect1.foo", &n);
	printf("sect1.foo: %d (%d)\n", n, i);
	
	i=kw_PushSection(input, "sect1");
	n=0;
	i=kw_GetInt(input, "foo", &n);
	printf("sect1.foo (push): %d (%d)\n", n, i);
	i=kw_PopSection(input);

	i=kw_GetIntArray(input, "sect2", &ia);
	printf("sect2.arg: [%d %d %d] (%d)\n", ia[0], ia[1], ia[2], i);
	free(ia);

	i=kw_GetData(input, "data", &str_a);
	printf("$data: (%d)\n", i);
	for (n=0; n < i; n++) {
		printf("%s\n", str_a[n]);
		free(str_a[n]);
	}
	free(str_a);
	printf("$end data\n");
*/	
	return 42;
}
