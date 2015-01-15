#include "Interpolation.hpp"

#include <cmath>
#include <cstdio>
#include <cstdlib>

#include "Vector3.hpp"
#include "Vector2.hpp"

using namespace std;

Interpolation::Interpolation(Vector3*** U, int gradeIn, const int type,
                             unsigned int noHlpElementsIn, unsigned int noPatchIn)
{
    noPatch = noPatchIn;
    h = (1<<noHlpElementsIn);
    noHlpElements = noHlpElementsIn;
    grade = (1<<gradeIn);
    if (gradeIn < 0) grade = 0;
    n = 1<<(noHlpElements-gradeIn);
    if (grade == 0) n = 1<<(noHlpElements);
    unsigned int i1, i2, i3, iSize, el;
    unsigned int i;

    //initialization -- maybe allocate in one go
    pSurfaceInterpolation = (Vector3****) malloc(noPatch*sizeof(
                                Vector3***));//+noPatch*n*(Vector3**) + noPatch*n*n*(Vector3*) + noPatch*n*n*(grade+*(grade+1)*(Vector3));
    for (i1=0; i1<noPatch; i1++) {
        pSurfaceInterpolation[i1] = (Vector3***) malloc(n*sizeof(Vector3**) + n*n*sizeof(
                                        Vector3*) + n*n*(grade+1)*(grade+1)*sizeof(Vector3));
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
                                pSurfaceInterpolation[i1][i2][i3][(grade+1)*iSize+i] = vector3SMul(h/el,
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
                                pSurfaceInterpolation[i1][i2][i3][(grade+1)*i+el]= vector3SMul(h/iSize,
                                        vector3Sub(pSurfaceInterpolation[i1][i2][i3][(grade+1)*i+el],
                                                   pSurfaceInterpolation[i1][i2][i3][(grade+1)*(i-1)+el]));
                            }
                        }
                    }
                }
            }
        }
    }
#ifdef DEBUG
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
}

// calculate the value in one point
Vector3 Interpolation::Chi(Vector2 a, int patch)
{
    unsigned int x, y;
    Vector3 c(0.0, 0.0, 0.0), d(0.0, 0.0, 0.0);
    Vector2 aTmp(a.y, a.x);
    a = aTmp;
  
    x = (unsigned int) floor((a.x*h)/grade);
    y = (unsigned int) floor((a.y*h)/grade);
  
    if(x > 0)if(x == (unsigned int)(n)) --x;
    if(y > 0)if(y == (unsigned int)(n)) --y;
	
    c = pSurfaceInterpolation[patch][x][y][grade*(grade+1)+grade];
    for(int j = grade-1; j >=0; --j) {
        c = vector3SMul(a.y-((y*grade+j)/h), c);
        c = vector3Add(c, pSurfaceInterpolation[patch][x][y][grade*(grade+1)+j]);
    }

    for(int i = grade-1; i >=0; --i) {
        d = pSurfaceInterpolation[patch][x][y][i*(grade+1)+grade];
        for(int j = grade-1; j >=0; --j) {
            d = vector3SMul(a.y-(y*grade+j)/h, d);
            d = vector3Add(d, pSurfaceInterpolation[patch][x][y][i*(grade+1)+j]);
        }
        c = vector3SMul(a.x-(x*grade+i)/h,c);
        c = vector3Add(c,d);
    }
    return c;
}

// calculate the x derivative in one point
Vector3 Interpolation::dChi_dx(Vector2 a, int patch){
  unsigned int x, y;
  double s,p;
  Vector2 aTmp(a.y, a.x);
  a = aTmp;
#ifdef DEBUG
  //printf("%lf %lf ", a.x, a.y);
#endif

  Vector3 c(0.0,0.0,0.0), dc_dx(0.0,0.0,0.0), dc_dy(0.0,0.0,0.0), res(0.0,0.0,0.0);

  x = floor((a.x*h)/grade);
  y = floor((a.y*h)/grade);
    
  //if(x > 0) if((grade*x/h - a.x)*(grade*x/h - a.x) < 1e-10) --x;
  //if(y > 0) if((grade*y/h - a.y)*(grade*y/h - a.y) < 1e-10) --y;
  /** @note that here there are two possibilities to adjust the "corner"
   * points, either they belong to the old interpolation or the new one. In my
   * case I chose to make it belong to the new one, except of course for the
   * last points which are here assigned to the old interpolation polinomial
   */
  if(x == (unsigned int)(n)) --x;
  if(y == (unsigned int)(n)) --y;

  dc_dy = pSurfaceInterpolation[patch][x][y][grade*(grade+1)+1];
	// calculate the other elements
	for(unsigned int der = 2; der <= grade; ++der){
		//coefficient of pSurfaceInterpolation[x][y][grade*(grade+1)+grade];
		s = 0;
		for(int i = der-1; i >=0; --i){
			p = 1;
			for(unsigned int j = 0; j <=der-1; ++j){
				if((int)j!=i){
					p *= (a.y-(y*grade+j)/h);
				}
			}
			s+=p;
		}
		dc_dy =vector3Add(dc_dy, vector3SMul(s,pSurfaceInterpolation[patch][x][y][grade*(grade+1)+der])) ;
	}
	for(int dy = grade-1; dy >= 0; --dy){
		res = pSurfaceInterpolation[patch][x][y][dy*(grade+1)+1];
		// calculate the other elements
		for(unsigned int der = 2; der <= grade; ++der){
			//coefficient of pSurfaceInterpolation[x][y][grade*(grade+1)+grade];
			s = 0;
			for(int i = der-1; i >=0; --i){
				// all combinations of i elements
				p = 1;
				for(unsigned int j = 0; j <=der-1; ++j){
					if((int)j!=i){
						p *= (a.y-(y*grade+j)/h);
					}
				}
				s+=p;
			}
			res = vector3Add(res, vector3SMul(s, pSurfaceInterpolation[patch][x][y][dy*(grade+1)+der]));
		}
		dc_dy = vector3SMul((a.x-(x*grade+dy)/h),dc_dy);
		dc_dy = vector3Add(res,dc_dy);
	}
		
    
#ifdef DEBUG
  //printf("%lf %lf %lf\n", c.x, c.y, c.z);
#endif
	//dc_dy = vector3SMul(0.5,dc_dy);
  return dc_dy;
}

// calculate the y derivative in one point
Vector3 Interpolation::dChi_dy(Vector2 a, int patch){
  unsigned int x, y;
  double s,p;
  Vector2 aTmp(a.y, a.x);
  a = aTmp;
#ifdef DEBUG
  //printf("%lf %lf ", a.x, a.y);
#endif

  Vector3 dc_dx(0.0,0.0,0.0), res(0.0,0.0,0.0);

  x = floor((a.x*h)/grade);
  y = floor((a.y*h)/grade);
    
  //if(x > 0) if((grade*x/h - a.x)*(grade*x/h - a.x) < 1e-10) --x;
  //if(y > 0) if((grade*y/h - a.y)*(grade*y/h - a.y) < 1e-10) --y;
  /** @note that here there are two possibilities to adjust the "corner"
   * points, either they belong to the old interpolation or the new one. In my
   * case I chose to make it belong to the new one, except of course for the
   * last points which are here assigned to the old interpolation polinomial
   */
  if(x == (unsigned int)(n)) --x;
  if(y == (unsigned int)(n)) --y;

  for(unsigned int dy = grade; dy >= 1; --dy){
		res = pSurfaceInterpolation[patch][x][y][dy*(grade+1)+grade];
		for(int j = grade-1; j >=0; --j){
			res = vector3SMul(a.y-((y*grade+j)/h), res);
			res = vector3Add(res, pSurfaceInterpolation[patch][x][y][dy*(grade+1)+j]);
		}
		s = 0;
		for(int i = dy-1; i >=0; --i){
			p = 1;
			for(unsigned int j = 0; j <=dy-1; ++j){
				if((int)j!=i){
					p *= (a.x-(x*grade+j)/h);
				}
			}
			s+=p;
		}
		dc_dx =vector3Add(dc_dx, vector3SMul(s,res)) ;
	}
		
#ifdef DEBUG
  //printf("%lf %lf %lf\n", c.x, c.y, c.z);
#endif
	//dc_dx = vector3SMul(0.5,dc_dx);
  return dc_dx;
}

// calculate the y second derivative in one point
Vector3 Interpolation::d2Chi_dy2(Vector2 a, int patch){
  unsigned int x, y;
  double s, s1, p;
  Vector2 aTmp(a.y, a.x);
  a = aTmp;
#ifdef DEBUG
  //printf("%lf %lf ", a.x, a.y);
#endif

  Vector3 dc_dx(0.0,0.0,0.0), res(0.0,0.0,0.0);

  x = floor((a.x*h)/grade);
  y = floor((a.y*h)/grade);
    
  //if(x > 0) if((grade*x/h - a.x)*(grade*x/h - a.x) < 1e-10) --x;
  //if(y > 0) if((grade*y/h - a.y)*(grade*y/h - a.y) < 1e-10) --y;
  /** @note that here there are two possibilities to adjust the "corner"
   * points, either they belong to the old interpolation or the new one. In my
   * case I chose to make it belong to the new one, except of course for the
   * last points which are here assigned to the old interpolation polinomial
   */
  if(x == (unsigned int)(1<<noHlpElements)) --x;
  if(y == (unsigned int)(1<<noHlpElements)) --y;

  for(unsigned int dy = 2; dy <= grade; dy++){
		res = pSurfaceInterpolation[patch][x][y][dy*(grade+1)+grade];
		for(int j = grade-1; j >=0; --j){
			res = vector3SMul(a.y-((y*grade+j)/h), res);
			res = vector3Add(res, pSurfaceInterpolation[patch][x][y][dy*(grade+1)+j]);
		}
		s = 0;
		for(int j = 0; j <= dy-1; ++j){
			s1 = 0;
      for(int k = 0; k <= dy-1; ++k){
        if ( k != j) {
          p = 1;
			    for(unsigned int i = 0; i <=dy-1; ++i){
				    if((j!=i) && (i != k)){
					    p *= (a.x-(x*grade+i)/h);
            }
          }
          s1+=p;
				}
			}
			s+=s1;
		}
		dc_dx =vector3Add(dc_dx, vector3SMul(s,res));
	}
		
#ifdef DEBUG
  //printf("%lf %lf %lf\n", c.x, c.y, c.z);
#endif
  //dc_dx = vector3SMul(0.25, dc_dx);
  return dc_dx;
}

// calculate the second derivative in x in one point
Vector3 Interpolation::d2Chi_dx2(Vector2 a, int patch){
  unsigned int x, y;
  double s, s1, p;
  Vector2 aTmp(a.y, a.x);
  a = aTmp;
#ifdef DEBUG
  //printf("%lf %lf ", a.x, a.y);
#endif

  Vector3 c(0.0,0.0,0.0), dc_dx(0.0,0.0,0.0), dc_dy(0.0,0.0,0.0), res(0.0,0.0,0.0);

  x = floor((a.x*h)/grade);
  y = floor((a.y*h)/grade);
    
  //if(x > 0) if((grade*x/h - a.x)*(grade*x/h - a.x) < 1e-10) --x;
  //if(y > 0) if((grade*y/h - a.y)*(grade*y/h - a.y) < 1e-10) --y;
  /** @note that here there are two possibilities to adjust the "corner"
   * points, either they belong to the old interpolation or the new one. In my
   * case I chose to make it belong to the new one, except of course for the
   * last points which are here assigned to the old interpolation polinomial
   */
  if(x == (unsigned int)(1<<noHlpElements)) --x;
  if(y == (unsigned int)(1<<noHlpElements)) --y;
/*
  dc_dx = pSurfaceInterpolation[patch][x][y][grade*(grade+1)+2];
	// calculate the other elements
	for(unsigned int der = 3; der <= grade; ++der){
		//coefficient of pSurfaceInterpolation[x][y][grade*(grade+1)+grade];
		s = 0;
		for(unsigned int i = 0; i <= der-1; ++i){
      s1 = 0;
      for(unsigned int k = 0; k <= der-1; ++k){
        if(i!=k){
			    // all combinations of i elements
			    p = 1;
			    for(unsigned int j = 0; j <=der-1; ++j){
				    if((j!=i) && (k!=j)){
					    p *= (a.y-(y*grade+j)/h);
            }
          }
          s1 +=p;
				}
			}
			s+=s1;
		}
		dc_dx =vector3Add(dc_dx, vector3SMul(s,pSurfaceInterpolation[patch][x][y][grade*(grade+1)+der])) ;
	}
	for(int dy = grade-1; dy >= 0; --dy){
		res = pSurfaceInterpolation[patch][x][y][dy*(grade+1)+2];
		// calculate the other elements
		for(unsigned int der = 3; der <= grade; ++der){
			//coefficient of pSurfaceInterpolation[x][y][grade*(grade+1)+grade];
			s = 0;
			for(unsigned int i = 0; i <= der-1; ++i){
        s1 = 0;
        for(unsigned int k = 0; k <= der-1; ++k){
          if(i!=k){
				    // all combinations of i elements
				    p = 1;
				    for(unsigned int j = 0; j <=der-1; ++j){
					    if((j!=i) && (k!=j)){
						    p *= (a.y-(y*grade+j)/h);
              }
            }
            s1 +=p;
					}
				}
				s+=s1;
			}
			res = vector3Add(res, vector3SMul(s, pSurfaceInterpolation[patch][x][y][dy*(grade+1)+der]));
		}
		dc_dx = vector3SMul((a.x-(x*grade+dy)/h),dc_dx);
		dc_dx = vector3Add(res,dc_dx);
	}
*/

  for(unsigned int dy = 2; dy <= grade; dy++){
		res = pSurfaceInterpolation[patch][x][y][grade*(grade+1)+dy];
		for(int j = grade-1; j >=0; --j){
			res = vector3SMul(a.x-((x*grade+j)/h), res);
			res = vector3Add(res, pSurfaceInterpolation[patch][x][y][j*(grade+1)+dy]);
		}
		s = 0;
		for(int j = 0; j <= dy-1; ++j){
			s1 = 0;
      for(int k = 0; k <= dy-1; ++k){
        if ( k != j) {
          p = 1;
			    for(unsigned int i = 0; i <=dy-1; ++i){
				    if((j!=i) && (i != k)){
					    p *= (a.y-(y*grade+i)/h);
            }
          }
          s1+=p;
				}
			}
			s+=s1;
		}
		dc_dx =vector3Add(dc_dx, vector3SMul(s,res));
	}

  //dc_dx = vector3SMul(0.25, dc_dx);
  return dc_dx;
}

// calculate the second derivative in one point
Vector3 Interpolation::d2Chi_dxy(Vector2 a, int patch){
  unsigned int x, y;
  double s,p, s2, p2;
  Vector2 aTmp(a.y, a.x);
  a = aTmp;
#ifdef DEBUG
  //printf("%lf %lf ", a.x, a.y);
#endif

  Vector3 c(0.0,0.0,0.0), dc_dx(0.0,0.0,0.0), dc_dy(0.0,0.0,0.0), res(0.0,0.0,0.0);

  x = floor((a.x*h)/grade);
  y = floor((a.y*h)/grade);
    
  //if(x > 0) if((grade*x/h - a.x)*(grade*x/h - a.x) < 1e-10) --x;
  //if(y > 0) if((grade*y/h - a.y)*(grade*y/h - a.y) < 1e-10) --y;
  /** @note that here there are two possibilities to adjust the "corner"
   * points, either they belong to the old interpolation or the new one. In my
   * case I chose to make it belong to the new one, except of course for the
   * last points which are here assigned to the old interpolation polinomial
   */
  if(x == (unsigned int)(n)) --x;
  if(y == (unsigned int)(n)) --y;

/*
  dc_dy = pSurfaceInterpolation[patch][x][y][grade*(grade+1)+1];

  // calculate the other elements
	for(unsigned int der = 2; der <= grade; ++der){
		//coefficient of pSurfaceInterpolation[x][y][grade*(grade+1)+grade];
		s = 0;
		for(int i = der-1; i >=0; --i){
			p = 1;
			for(unsigned int j = 0; j <=der-1; ++j){
				if((int)j!=i){
					p *= (a.y-(y*grade+j)/h);
				}
			}
			s+=p;
		}
		dc_dy =vector3Add(dc_dy, vector3SMul(s,pSurfaceInterpolation[patch][x][y][grade*(grade+1)+der])) ;
	}
*/
	for(int dy = grade; dy >= 1; --dy){
		res = pSurfaceInterpolation[patch][x][y][dy*(grade+1)+1];
		// calculate the other elements
		for(unsigned int der = 2; der <= grade; ++der){
			//coefficient of pSurfaceInterpolation[x][y][grade*(grade+1)+grade];
			s = 0;
			for(int i = der-1; i >=0; --i){
				// all combinations of i elements
				p = 1;
				for(unsigned int j = 0; j <=der-1; ++j){
					if((int)j!=i){
						p *= (a.y-(y*grade+j)/h);
					}
				}
				s+=p;
			}
			res = vector3Add(res, vector3SMul(s, pSurfaceInterpolation[patch][x][y][dy*(grade+1)+der]));
		}
    s2 = 0;
    for(int i = dy-1; i >=0; --i){
      // all combinations of i elements
      p2 = 1;
      for(unsigned int j = 0; j <=dy-1;++j){
        if((int)j!=i){
          p2*=(a.x-(x*grade+j)/h);
        }
      }
      s2+=p2;
    }
		//dc_dy = vector3SMul((a.x-(x*grade+dy)/h),dc_dy);
		res = vector3SMul(s2,res);
		dc_dy = vector3Add(res,dc_dy);
	}
  //dc_dy = vector3SMul(0.25, dc_dy);
	
  return dc_dy;
}
// calculate the normal in one point
Vector3 Interpolation::n_Chi(Vector2 a, int patch)
{
    unsigned int x, y;
    double s,p;
    Vector2 aTmp(a.y, a.x);
    a = aTmp;

    Vector3 c(0.0,0.0,0.0), dc_dx(0.0,0.0,0.0), dc_dy(0.0,0.0,0.0), res(0.0,0.0,0.0);

    x = floor((a.x*h)/grade);
    y = floor((a.y*h)/grade);

    //if(x > 0) if((grade*x/h - a.x)*(grade*x/h - a.x) < 1e-10) --x;
    //if(y > 0) if((grade*y/h - a.y)*(grade*y/h - a.y) < 1e-10) --y;
    /** @note that here there are two possibilities to adjust the "corner"
     * points, either they belong to the old interpolation or the new one. In my
     * case I chose to make it belong to the new one, except of course for the
     * last points which are here assigned to the old interpolation polinomial
     */
    if(x == ((unsigned int)n)) --x;
    if(y == ((unsigned int)n)) --y;

    dc_dy = pSurfaceInterpolation[patch][x][y][grade*(grade+1)+1];
    // calculate the other elements
    for(unsigned int der = 2; der <= grade; ++der) {
        //coefficient of pSurfaceInterpolation[x][y][grade*(grade+1)+grade];
        s = 0;
        for(int i = der-1; i >=0; --i) {
            p = 1;
            for(unsigned int j = 0; j <=der-1; ++j) {
                if((int)j!=i) {
                    p *= (a.y-(y*grade+j)/h);
                }
            }
            s+=p;
        }
        dc_dy =vector3Add(dc_dy, vector3SMul(s,
                                             pSurfaceInterpolation[patch][x][y][grade*(grade+1)+der])) ;
    }
    for(int dy = grade-1; dy >= 0; --dy) {
        res = pSurfaceInterpolation[patch][x][y][dy*(grade+1)+1];
        // calculate the other elements
        for(unsigned int der = 2; der <= grade; ++der) {
            //coefficient of pSurfaceInterpolation[x][y][grade*(grade+1)+grade];
            s = 0;
            for(int i = der-1; i >=0; --i) {
                // all combinations of i elements
                p = 1;
                for(unsigned int j = 0; j <=der-1; ++j) {
                    if((int)j!=i) {
                        p *= (a.y-(y*grade+j)/h);
                    }
                }
                s+=p;
            }
            res = vector3Add(res, vector3SMul(s,
                                              pSurfaceInterpolation[patch][x][y][dy*(grade+1)+der]));
        }
        dc_dy = vector3SMul((a.x-(x*grade+dy)/h),dc_dy);
        dc_dy = vector3Add(res,dc_dy);
    }

    for(unsigned int dy = 1; dy <= grade; dy++) {
        res = pSurfaceInterpolation[patch][x][y][dy*(grade+1)+grade];
        for(int j = grade-1; j >=0; --j) {
            res = vector3SMul(a.y-((y*grade+j)/h), res);
            res = vector3Add(res, pSurfaceInterpolation[patch][x][y][dy*(grade+1)+j]);
        }
        s = 0;
        for(int i = dy-1; i >=0; --i) {
            p = 1;
            for(unsigned int j = 0; j <=dy-1; ++j) {
                if((int)j!=i) {
                    p *= (a.x-(x*grade+j)/h);
                }
            }
            s+=p;
        }
        dc_dx =vector3Add(dc_dx, vector3SMul(s,res)) ;
    }

    c.x = (-dc_dy.y*dc_dx.z + dc_dy.z*dc_dx.y);
    c.y = (-dc_dy.z*dc_dx.x + dc_dy.x*dc_dx.z);
    c.z = (-dc_dy.x*dc_dx.y + dc_dy.y*dc_dx.x);

    return c;
}

Interpolation::~Interpolation()
{
    for(unsigned int i = 0; i < noPatch; ++i) free(pSurfaceInterpolation[i]);
    free(pSurfaceInterpolation);
}


