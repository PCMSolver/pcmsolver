#include "Vector3.hpp"
#include "Vector2.hpp"
#include <stdlib.h>
#include <string.h>
#ifdef DEBUG
  #include <cstdio>
#endif
#include "Interpolation.hpp"
Interpolation::Interpolation(Vector3*** U, int gradeIn, const int /* type */, unsigned int nLevelsIn, unsigned int noPatchIn){
  noPatch = noPatchIn;
  n = 1<<(nLevelsIn-gradeIn);
  grade = (1<<gradeIn);
  coeff = (Vector3*) malloc(sizeof(Vector3)*(grade+1));
  if (gradeIn < 0) grade = 0;
  h = 1./grade;
  if (grade == 0) n = 1<<(nLevelsIn);
  unsigned int i1, i2, i3, iSize, el;
  unsigned int i;

  //initialization -- maybe allocate in one go
  pSurfaceInterpolation = (Vector3****) malloc(noPatch*sizeof
      (Vector3***));//+noPatch*n*(Vector3**) + noPatch*n*n*(Vector3*) + noPatch*n*n*(grade+*(grade+1)*(Vector3));
  for (i1=0; i1<noPatch; i1++){
    pSurfaceInterpolation[i1] = (Vector3***) malloc(n*sizeof(Vector3**) + n*n*sizeof
        (Vector3*) + n*n*(grade+1)*(grade+1)*sizeof(Vector3));
    for (i2=0; i2 < n; ++i2) { // rowwise counting of the patch zi = (i1,i2,i3)
      pSurfaceInterpolation[i1][i2] = (Vector3**) (pSurfaceInterpolation[i1]+n)+i2*n;
      for (i3=0; i3<n; ++i3) 	{
        pSurfaceInterpolation[i1][i2][i3] = (Vector3*) (pSurfaceInterpolation[i1]+n+n*n)
          +i2*n*(grade+1)*(grade+1) + i3*(grade+1)*(grade+1);
        for (iSize = 0; iSize <= grade; iSize++) {
          for( el = 0; el <= grade; el++) {
            //calculate 1D interpolation polinomials
            if(grade == 0) {
							pSurfaceInterpolation[i1][i2][i3][(grade+1)*iSize+el] = U[i1][i2][i3];
						} else {
							pSurfaceInterpolation[i1][i2][i3][(grade+1)*iSize+el] =
                U[i1][grade*i2+iSize][grade*i3+el];
						}
					}
        }
        if(grade!=0) {
					for (iSize = 0; iSize <= grade; iSize++) {
						for( el = 1; el <= grade; el++) {
							for(i = grade; i >= el; --i) {
								pSurfaceInterpolation[i1][i2][i3][(grade+1)*iSize+i] = vector3SMul(1./(h*el),
                    vector3Sub(pSurfaceInterpolation[i1][i2][i3][(grade+1)*iSize+i],
                      pSurfaceInterpolation[i1][i2][i3][(grade+1)*iSize+i-1]));
							}
						}
					}
					for(iSize = 1; iSize <= grade; ++iSize) {
						// column iSize - Newton Scheme on each column
						for(el = 0; el <= grade; ++el) {
							// Newton Scheme for calculating interpolation polinomial
							for( i = grade; i >= iSize; --i) {
								pSurfaceInterpolation[i1][i2][i3][(grade+1)*i+el]= vector3SMul(1./(h*iSize),
                    vector3Sub(pSurfaceInterpolation[i1][i2][i3][(grade+1)*i+el],
                      pSurfaceInterpolation[i1][i2][i3][(grade+1)*(i-1)+el]));
							}
						}
					}
				}
      }
    }
  }
#ifdef DEBUG2
  FILE * debugFile = fopen("debug.out", "a");
  fprintf(debugFile, ">>> INTERPOLATION_POLINOMIALS\n");
  fprintf(debugFile, "%d %lf\n", grade, h);
  for (i1=0; i1<noPatch; ++i1) {
    for (i2=0; i2 < n; ++i2) { // rowwise counting of the patch zi = (i1,i2,i3)
      for (i3=0; i3<n; ++i3) 	{
        for (iSize = 0; iSize <= grade; iSize++) {
          for( el = 0; el <= grade; el++) {
            fprintf(debugFile, "%d %d %d  %d %lf %lf %lf\n",i1, i2, i3, (grade+1)*iSize+el,
                pSurfaceInterpolation[i1][i2][i3][(grade+1)*iSize+el].x,
                pSurfaceInterpolation[i1][i2][i3][(grade+1)*iSize+el].y,
                pSurfaceInterpolation[i1][i2][i3][(grade+1)*iSize+el].z);
					}
				}
				fprintf(debugFile, "\n");
			}
		}
	}
  fprintf(debugFile, "<<< INTERPOLATION_POLINOMIALS\n");
  fclose(debugFile);
#endif
  return;
}

// calculate the value in one point
Vector3 Interpolation::Chi(Vector2 a, int patch){
	unsigned int x, y;
  Vector3 c(0.0, 0.0, 0.0), d(0.0, 0.0, 0.0);

  x = (unsigned int) floor(a.x*n);
  y = (unsigned int) floor(a.y*n);

  if(x == n) --x;
  if(y == n) --y;
  a.x = n*a.x - floor(a.x*n);
  a.y = n*a.y - floor(a.y*n);
	//if(x > 0) if((grade*x/h - a.x)*(grade*x/h - a.x) < 1e-10) --x;
  //if(y > 0) if((grade*y/h - a.y)*(grade*y/h - a.y) < 1e-10) --y;

	c = pSurfaceInterpolation[patch][y][x][grade*(grade+1)+grade];
	for(int j = grade-1; j >=0; --j) {
		c = vector3SMul(a.x-j*h, c);
		c = vector3Add(c, pSurfaceInterpolation[patch][y][x][grade*(grade+1)+j]);
	}

	for(int i = grade-1; i >=0; --i) {
		d = pSurfaceInterpolation[patch][y][x][i*(grade+1)+grade];
		for(int j = grade-1; j >=0; --j) {
			d = vector3SMul(a.x-j*h, d);
			d = vector3Add(d, pSurfaceInterpolation[patch][y][x][i*(grade+1)+j]);
		}
		c = vector3SMul(a.y-i*h,c);
		c = vector3Add(c,d);
	}
  return c;
}

// calculate the x derivative in one point
Vector3 Interpolation::dChi_dx(Vector2 a, int patch){
  unsigned int x, y;

  Vector3 c(0.0,0.0,0.0), dc_dx(0.0,0.0,0.0);

  x = (unsigned int)floor(a.x*n);
  y = (unsigned int)floor(a.y*n);

  //if(x > 0) if((grade*x/h - a.x)*(grade*x/h - a.x) < 1e-10) --x;
  //if(y > 0) if((grade*y/h - a.y)*(grade*y/h - a.y) < 1e-10) --y;
  /** @note that here there are two possibilities to adjust the "corner"
   * points, either they belong to the old interpolation or the new one. In my
   * case I chose to make it belong to the new one, except of course for the
   * last points which are here assigned to the old interpolation polinomial
   */
  if(x == n) --x;
  if(y == n) --y;
  a.x = n*a.x - floor(a.x*n);
  a.y = n*a.y - floor(a.y*n);

  coeff[grade] = pSurfaceInterpolation[patch][y][x][grade*(grade+1)+grade];
  dc_dx = coeff[grade];
	for(int j = grade-1; j >=1; --j) {
    coeff[j] = vector3Add( pSurfaceInterpolation[patch][y][x][grade*(grade+1)+j],vector3SMul((a.x-j*h), coeff[j+1]));
		dc_dx = vector3SMul(a.x-(j-1)*h, dc_dx);
		dc_dx = vector3Add(dc_dx, coeff[j]);
	}

	for(int i = grade-1; i >=0; --i) {
    coeff[grade] = pSurfaceInterpolation[patch][y][x][i*(grade+1)+grade];
		c = coeff[grade];
		for(int j = grade-1; j >=1; --j) {
      coeff[j] = vector3Add(pSurfaceInterpolation[patch][y][x][i*(grade+1)+j],vector3SMul((a.x-j*h), coeff[j+1]));
			c = vector3SMul(a.x-(j-1)*h, c);
			c = vector3Add(c, coeff[j]);
		}
		dc_dx = vector3SMul(a.y-i*h,dc_dx);
		dc_dx = vector3Add(dc_dx,c);
	}
  dc_dx = vector3SMul(n,dc_dx);
  return dc_dx;
}

// calculate the y derivative in one point
Vector3 Interpolation::dChi_dy(Vector2 a, int patch){
  unsigned int x, y;

  Vector3 c(0.0,0.0,0.0), dc_dx(0.0,0.0,0.0), dc_dy(0.0,0.0,0.0), res(0.0,0.0,0.0);

  x = (unsigned int)floor(a.x*n);
  y = (unsigned int)floor(a.y*n);

  //if(x > 0) if((grade*x/h - a.x)*(grade*x/h - a.x) < 1e-10) --x;
  //if(y > 0) if((grade*y/h - a.y)*(grade*y/h - a.y) < 1e-10) --y;
  /** @note that here there are two possibilities to adjust the "corner"
   * points, either they belong to the old interpolation or the new one. In my
   * case I chose to make it belong to the new one, except of course for the
   * last points which are here assigned to the old interpolation polinomial
   */
  if(x == n) --x;
  if(y == n) --y;
  a.x = n*a.x - floor(a.x*n);
  a.y = n*a.y - floor(a.y*n);

  coeff[grade] = pSurfaceInterpolation[patch][y][x][grade*(grade+1)+grade];
  dc_dy = coeff[grade];
	for(int j = grade-1; j >=0; --j) {
    coeff[j] = pSurfaceInterpolation[patch][y][x][grade*(grade+1)+j];
		dc_dy = vector3SMul(a.x-((j)*h), dc_dy);
		dc_dy = vector3Add(dc_dy, coeff[j]);
	}

	for(int i = grade-1; i >=1; --i) {
    coeff[grade] = vector3Add(pSurfaceInterpolation[patch][y][x][i*(grade+1)+grade],vector3SMul((a.y-(i)*h),coeff[grade]));
		c = coeff[grade];
		for(int j = grade-1; j >=0; --j) {
      coeff[j] = vector3Add(pSurfaceInterpolation[patch][y][x][i*(grade+1)+j],vector3SMul((a.y-(i)*h),coeff[j]));
			c = vector3SMul(a.x-(j)*h, c);
			c = vector3Add(c, coeff[j]);
		}
		dc_dy = vector3SMul(a.y-(i-1)*h,dc_dy);
		dc_dy = vector3Add(dc_dy,c);
	}
  dc_dy = vector3SMul(n,dc_dy);
  return dc_dy;
}

// calculate the y second derivative in one point
Vector3 Interpolation::d2Chi_dy2(Vector2 a, int patch){
  unsigned int x, y;
  double s, s1, p;

  Vector3 dc_dx(0.0,0.0,0.0), res(0.0,0.0,0.0);

  x = (unsigned int)floor(a.x*n);
  y = (unsigned int)floor(a.y*n);

  //if(x > 0) if((grade*x/h - a.x)*(grade*x/h - a.x) < 1e-10) --x;
  //if(y > 0) if((grade*y/h - a.y)*(grade*y/h - a.y) < 1e-10) --y;
  /** @note that here there are two possibilities to adjust the "corner"
   * points, either they belong to the old interpolation or the new one. In my
   * case I chose to make it belong to the new one, except of course for the
   * last points which are here assigned to the old interpolation polinomial
   */
  if(x == n) --x;
  if(y == n) --y;
  a.x = n*a.x - floor(a.x*n);
  a.y = n*a.y - floor(a.y*n);

  for(unsigned int dy = 2; dy <= grade; dy++){
		res = pSurfaceInterpolation[patch][y][x][dy*(grade+1)+grade];
		for(int j = grade-1; j >=0; --j){
			res = vector3SMul(a.x-j*h, res);
			res = vector3Add(res, pSurfaceInterpolation[patch][y][x][dy*(grade+1)+j]);
		}
		s = 0;
		for(unsigned int j = 0; j <= dy-1; ++j){
			s1 = 0;
      for(unsigned int k = 0; k <= dy-1; ++k){
        if ( k != j) {
          p = 1;
			    for(unsigned int i = 0; i <=dy-1; ++i){
				    if((j!=i) && (i != k)){
					    p *= (a.y-i*h);
            }
          }
          s1+=p;
				}
			}
			s+=s1;
		}
		dc_dx =vector3Add(dc_dx, vector3SMul(s,res));
	}

  return dc_dx;
}

// calculate the second derivative in x in one point
Vector3 Interpolation::d2Chi_dx2(Vector2 a, int patch){
  unsigned int x, y;
  double s, s1, p;

  Vector3 c(0.0,0.0,0.0), dc_dx(0.0,0.0,0.0), dc_dy(0.0,0.0,0.0), res(0.0,0.0,0.0);

  x = (unsigned int)floor(a.x*n);
  y = (unsigned int)floor(a.y*n);
  memcpy(coeff, pSurfaceInterpolation[patch][y][x], (grade+1)*(grade+1)*sizeof(Vector3));

  //if(x > 0) if((grade*x/h - a.x)*(grade*x/h - a.x) < 1e-10) --x;
  //if(y > 0) if((grade*y/h - a.y)*(grade*y/h - a.y) < 1e-10) --y;
  /** @note that here there are two possibilities to adjust the "corner"
   * points, either they belong to the old interpolation or the new one. In my
   * case I chose to make it belong to the new one, except of course for the
   * last points which are here assigned to the old interpolation polinomial
   */
  if(x == n) --x;
  if(y == n) --y;
  a.x = n*a.x - floor(a.x*n);
  a.y = n*a.y - floor(a.y*n);

  for(unsigned int dy = 2; dy <= grade; dy++){
		res = pSurfaceInterpolation[patch][y][x][grade*(grade+1)+dy];
		for(int j = grade-1; j >=0; --j){
			res = vector3SMul(a.y-j*h, res);
			res = vector3Add(res, pSurfaceInterpolation[patch][y][x][j*(grade+1)+dy]);
		}
		s = 0;
		for(unsigned int j = 0; j <= dy-1; ++j){
			s1 = 0;
      for(unsigned int k = 0; k <= dy-1; ++k){
        if ( k != j) {
          p = 1;
			    for(unsigned int i = 0; i <=dy-1; ++i){
				    if((j!=i) && (i != k)){
					    p *= (a.x-i*h);
            }
          }
          s1+=p;
				}
			}
			s+=s1;
		}
		dc_dx =vector3Add(dc_dx, vector3SMul(s,res));
	}

  return dc_dx;
}

// calculate the second derivative in one point
Vector3 Interpolation::d2Chi_dxy(Vector2 a, int patch){
  unsigned int x, y;
  double s,p, s2, p2;

  Vector3 c(0.0,0.0,0.0), dc_dx(0.0,0.0,0.0), dc_dy(0.0,0.0,0.0), res(0.0,0.0,0.0);

  x = (unsigned int)floor(a.x*n);
  y = (unsigned int)floor(a.y*n);

  //if(x > 0) if((grade*x/h - a.x)*(grade*x/h - a.x) < 1e-10) --x;
  //if(y > 0) if((grade*y/h - a.y)*(grade*y/h - a.y) < 1e-10) --y;
  /** @note that here there are two possibilities to adjust the "corner"
   * points, either they belong to the old interpolation or the new one. In my
   * case I chose to make it belong to the new one, except of course for the
   * last points which are here assigned to the old interpolation polinomial
   */
  if(x == n) --x;
  if(y == n) --y;
  a.x = n*a.x - floor(a.x*n);
  a.y = n*a.y - floor(a.y*n);

	for(int dy = grade; dy >= 1; --dy){
		res = pSurfaceInterpolation[patch][y][x][dy*(grade+1)+1];
		// calculate the other elements
		for(unsigned int der = 2; der <= grade; ++der){
			//coefficient of pSurfaceInterpolation[y][x][grade*(grade+1)+grade];
			s = 0;
			for(int i = der-1; i >=0; --i){
				// all combinations of i elements
				p = 1;
				for(unsigned int j = 0; j <=der-1; ++j){
					if((int)j!=i){
						p *= (a.x-j*h);
					}
				}
				s+=p;
			}
			res = vector3Add(res, vector3SMul(s, pSurfaceInterpolation[patch][y][x][dy*(grade+1)+der]));
		}
    s2 = 0;
    for(int i = dy-1; i >=0; --i){
      // all combinations of i elements
      p2 = 1;
      for(int j = 0; j <=dy-1;++j){
        if(j!=i){
          p2*=(a.y-j*h);
        }
      }
      s2+=p2;
    }
		//dc_dy = vector3SMul((a.y-dy*h),dc_dy);
		res = vector3SMul(s2,res);
		dc_dy = vector3Add(res,dc_dy);
	}

  return dc_dy;
}
// calculate the normal in one point
Vector3 Interpolation::n_Chi(Vector2 a, int patch){
  unsigned int x, y;

  Vector3 c(0.0,0.0,0.0), dc_dx(0.0,0.0,0.0), dc_dy(0.0,0.0,0.0), res(0.0,0.0,0.0);

  x = floor(a.x*n);
  y = floor(a.y*n);

  //if(x > 0) if((grade*x/h - a.x)*(grade*x/h - a.x) < 1e-10) --x;
  //if(y > 0) if((grade*y/h - a.y)*(grade*y/h - a.y) < 1e-10) --y;
  /** @note that here there are two possibilities to adjust the "corner"
   * points, either they belong to the old interpolation or the new one. In my
   * case I chose to make it belong to the new one, except of course for the
   * last points which are here assigned to the old interpolation polinomial
   */
  if(x == n) --x;
  if(y == n) --y;
  a.x = n*a.x - floor(a.x*n);
  a.y = n*a.y - floor(a.y*n);

  coeff[grade] = pSurfaceInterpolation[patch][y][x][grade*(grade+1)+grade];
  dc_dy = coeff[grade];
	for(int j = grade-1; j >=0; --j) {
    coeff[j] = pSurfaceInterpolation[patch][y][x][grade*(grade+1)+j];
		dc_dy = vector3SMul(a.x-((j)*h), dc_dy);
		dc_dy = vector3Add(dc_dy, coeff[j]);
	}

	for(int i = grade-1; i >=1; --i) {
    coeff[grade] = vector3Add(pSurfaceInterpolation[patch][y][x][i*(grade+1)+grade],vector3SMul((a.y-(i)*h),coeff[grade]));
		c = coeff[grade];
		for(int j = grade-1; j >=0; --j) {
      coeff[j] = vector3Add(pSurfaceInterpolation[patch][y][x][i*(grade+1)+j],vector3SMul((a.y-(i)*h),coeff[j]));
			c = vector3SMul(a.x-(j)*h, c);
			c = vector3Add(c, coeff[j]);
		}
		dc_dy = vector3SMul(a.y-(i-1)*h,dc_dy);
		dc_dy = vector3Add(dc_dy,c);
	}
  dc_dy = vector3SMul(n,dc_dy);

  coeff[grade] = pSurfaceInterpolation[patch][y][x][grade*(grade+1)+grade];
  dc_dx = coeff[grade];
	for(int j = grade-1; j >=1; --j) {
    coeff[j] = vector3Add( pSurfaceInterpolation[patch][y][x][grade*(grade+1)+j],vector3SMul((a.x-(j)*h), coeff[j+1]));
		dc_dx = vector3SMul(a.x-((j-1)*h), dc_dx);
		dc_dx = vector3Add(dc_dx, coeff[j]);
	}

	for(int i = grade-1; i >=0; --i) {
    coeff[grade] = pSurfaceInterpolation[patch][y][x][i*(grade+1)+grade];
		c = coeff[grade];
		for(int j = grade-1; j >=1; --j) {
      coeff[j] = vector3Add(pSurfaceInterpolation[patch][y][x][i*(grade+1)+j],vector3SMul((a.x-(j)*h), coeff[j+1]));
			c = vector3SMul(a.x-(j-1)*h, c);
			c = vector3Add(c, coeff[j]);
		}
		dc_dx = vector3SMul(a.y-(i)*h,dc_dx);
		dc_dx = vector3Add(dc_dx,c);
	}
  dc_dx = vector3SMul(n,dc_dx);

	//c.x = (dc_dy.y*dc_dx.z - dc_dy.z*dc_dx.y);
	//c.y = (dc_dy.z*dc_dx.x - dc_dy.x*dc_dx.z);
	//c.z = (dc_dy.x*dc_dx.y - dc_dy.y*dc_dx.x);
	c.x = (dc_dy.z*dc_dx.y - dc_dy.y*dc_dx.z);
	c.y = (dc_dy.x*dc_dx.z - dc_dy.z*dc_dx.x);
	c.z = (dc_dy.y*dc_dx.x - dc_dy.x*dc_dx.y);

  return c;
}

Interpolation::~Interpolation(){
  for(unsigned int i = 0; i < noPatch; ++i) free(pSurfaceInterpolation[i]);
  free(pSurfaceInterpolation);
  free(coeff);
}


