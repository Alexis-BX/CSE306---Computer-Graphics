#include "vector.cpp"

enum Materials{
    opaque,
    miror,
    transparent
};

struct TriangleIndices {
    int vtxindices[3] ;
    int normalindices[3] ;
    int uvindices[3] ;
};

class Mesh {
public:
    vector<Vector> vertices;
    vector<Vector> normals;
    vector<Vector> uvs;
    vector<TriangleIndices> triangles;
};

class Sphere{
public:
    Sphere(Vector& center, double r, Vector& color, Materials material=opaque){
        R = r;
        p = center;
        c = color;
        m = material;
    }
    Vector p, c;
    double R;
    Materials m;
};

class Ray {
public:
    Ray(Vector& point, Vector& direction){
        p = point;
        d = direction;
    }
    Vector p, d;
};

class Camera {
public:
    Camera(){}
    Camera(Vector& point, double fieldOfView, double aperatureRad=.1){
        p = point;
        fov = fieldOfView;
        f = w/(2.*tan(fov/2.*PI/180.));
        R = aperatureRad;
    }
    Vector p;
    double fov, f, R;
};

class Light {
public:
    Light(){}
    Light(Vector& point, int intensity, double radius){
        p = point;
        I = intensity;
        R = radius;
    }
    Vector p;
    double R;
    int I;
};