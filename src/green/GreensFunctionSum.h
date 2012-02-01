#ifndef GREENSFUNCTIONSUM
#define GREENSFUNCTIONSUM

// I am for now limiting myself to one template parameter
// In principle three templates are necessary here
// The question is whether it is worth implementng the full version.

template<class T>
class GreensFunctionSum : public GreensFunction<T>
{
 public:
    GreensFunctionSum(GreensFunction<T> &first, GreensFunction<T> &second);
    GreensFunctionSum(Section green);
    ~GreensFunctionSum(){delete greenFirst; delete greenSecond;}
    double evald(Vector3d &direction, Vector3d &p1, Vector3d &p2);
    void compDiagonalElementS(double area);
    void compDiagonalElementD(double area, double radius);
 protected:
    T evalGreensFunction(T * source, T * probe);
    GreensFunctionInterface* greenFirst;
    GreensFunctionInterface* greenSecond;
};
#endif
