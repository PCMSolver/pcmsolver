#ifndef VACUUM
#define VACUUM

template<typename T>
class Vacuum : public GreensFunction<T>
{
 public:
    Vacuum(){GreensFunction<T>::uniformFlag = true;};
    ~Vacuum(){};
    double evald(Vector3d & direction, Vector3d & p1, Vector3d & p2);
    double compDiagonalElementS(double area);
    double compDiagonalElementD(double area, double radius);
 private:
    T evalGreensFunction(T * source, T * probe);
};

#endif
