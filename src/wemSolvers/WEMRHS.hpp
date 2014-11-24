class GenericAnsatzFunction;

void WEMRHS1(double **rhs, GenericAnsatzFunction *af);
void WEMRHS2(double **rhs, GenericAnsatzFunction *af);
void WEMRHS2M(double **rhs, double *potential, GenericAnsatzFunction *af);
void WEMRHS2M_test(double **rhs, double *potential, GenericAnsatzFunction *af);
