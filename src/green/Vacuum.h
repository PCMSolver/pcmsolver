#ifndef VACUUM
#define VACUUM

class Vacuum : public GreensFunction
{
 public:
    Vacuum(){uniformFlag = true;};
    ~Vacuum(){};
    double evalf(Vector3d & p1, Vector3d & p2);
    double evald(Vector3d & direction, Vector3d & p1, Vector3d & p2, double delta = 0.001);
    double derivative(Vector3d & direction, Vector3d & p1, Vector3d & p2, double delta = 0.001);
    void gradient(Vector3d & gradient, Vector3d & p1, Vector3d & p2, double delta = 0.001);
    // private:
    //    double def;
};

#endif
