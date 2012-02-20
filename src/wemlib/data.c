/************
 *  data.c  *
 ************/


#include <math.h>
#include "vector3.h"
#include "data.h"


/* number of charges */
const unsigned int k = 12;

/* position of charges: vectors (x,y,z) */
const vector3 x[12] = {
    {9.96641558322, 3.77756252386, -16.1911734389},
    {12.5232150303, 3.81346732024, -15.5127617601},
    {13.9197226367, 1.56658295762, -15.4995336772},
    {12.759430796, -0.716206201373,-16.1647172732},
    {10.2026313488,-0.752110997747,-16.8412392259},
    {8.80612374247, 1.49477336487, -16.8544673088},
    {8.88927169197, 5.51044138048, -16.2006220696},
    {13.4189452136, 5.57469206873, -15.0006459802},
    {15.8925967112, 1.59492884949, -14.9779692667},
    {13.8365746872,-2.44908505799, -16.1533789164},
    {9.30690116559,-2.51333574624, -17.3533550058},
    {6.83324966798, 1.466427473,   -17.3779214454}
};

/* charges: double values */
const double alpha[12] = { 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1 };


double f(a)
vector3 a;
{
    double c = 0;
    unsigned int i;
    for (i = 0; i < k; i++)
        c += alpha[i] / vector3_norm(vector3_make(a.x - x[i].x, a.y - x[i].y, a.z - x[i].z));
    return (c);
}


vector3 df(a)
vector3 a;
{
    unsigned int i;
    vector3 c, r, v;
    c.x = c.y = c.z = 0;
    for (i = 0; i < k; i++) {
        r = vector3_make(a.x - x[i].x, a.y - x[i].y, a.z - x[i].z);
        v = vector3_Smul(alpha[i] / pow(vector3_norm(r), 3), r);
        c.x -= v.x;
        c.y -= v.y;
        c.z -= v.z;
    }
    return (c);
}
