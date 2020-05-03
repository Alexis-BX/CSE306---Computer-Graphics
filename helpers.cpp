#ifndef HELP
#define HELP

#include "simple_obj_file_reader.cpp"

const double gamma {2.2}; // usual: 2.2
const double rgbCorrection = pow(255, (gamma-1)/gamma);
const Vector NULLVEC(0,0,0);

Camera cam;
vector<TriangleMesh> scene;
Vector getColor(Vector p, int si, int sj, Light light, int depth=10, Vector previous=cam.p);

double random(){
    return double(rand()) / double(RAND_MAX);
}

void writePPM(int w, int h, int* image){
    ofstream output;
    output.open("pic.ppm");
    output << "P3" << "\n";
    output << w << " " << h << "\n";
    output << "255\n";

    for(int i=0; i<h; i++){
        for(int j=0; j<w; j++){
            int tmp = i*w*3 + j*3;
            output << image[tmp] << " " << image[tmp+1] << " " << image[tmp+2] << "\n";
        }
    }

    output.close();
}

void buildScene(){
    TriangleMesh mesh;
    mesh.readOBJ("model_cat/cat.obj");
    scene.push_back(mesh);
}

Vector getPixCoord(double i, double j){
    double x = i;
    double y = h-j-1;
    return Vector(cam.p[0]+x+0.5-w/2., cam.p[1]+y+0.5-h/2., cam.p[2]-cam.f);
}

Vector intersect(Ray ray, Vector A, Vector B, Vector C){
    Vector N = normalize(cross(A-B, A-C));

    if(dot(N, ray.d)==0) return NULLVEC;
    double t = dot(N,A-ray.p) / dot(N, ray.d);
    if (t < 0) return NULLVEC;
    Vector ret = ray.p + t * ray.d;

    Vector tmp = cross(B - A, ret - A);
    if (dot(N,tmp) < 0) return NULLVEC;

    tmp = cross(C - B, ret - B);
    if (dot(N,tmp) < 0)  return NULLVEC;

    tmp = cross(A - C, ret - C);
    if (dot(N,tmp) < 0) return NULLVEC;
    return ret;
}

struct sphereIpointP{
    int i = -1;
    int j = -1;
    Vector inter = NULLVEC;
};

sphereIpointP intersectMesh(Ray ray, int i){
    std::vector<TriangleIndices> mesh = scene[i].indices;
    sphereIpointP res;
    double closest = 0;

    for(int j=0; j<int(mesh.size()); j++){
        TriangleIndices tmp = mesh[j];
        Vector A = scene[i].vertices[tmp.vtxi];
        Vector B = scene[i].vertices[tmp.vtxj];
        Vector C = scene[i].vertices[tmp.vtxk];
        Vector inter = intersect(ray, A, B, C);
        if (inter!=NULLVEC){
            double d = norm(cam.p-inter);
            if(d < closest || res.j==-1){
                res.j = j;
                res.inter = inter;
                closest = d;
            }

        }
    }
    if (res.j!=-1) res.i=i;
    return res;
}

sphereIpointP intersectScene(Ray ray,int skip=-1){
    sphereIpointP res;
    double closest = 0;

    for(int i=0; i<int(scene.size()); i++){
        if (i==skip) continue;
        sphereIpointP tmp = intersectMesh(ray, i);
        if (tmp.i!=-1){
            double d = norm(cam.p-tmp.inter);
            if(d < closest || res.i==-1){
                res.i = tmp.i;
                res.j = tmp.j;
                res.inter = tmp.inter;
                closest = d;
            }
        }
        
    }
    return res;
}

double visibility(Vector p, Vector light, int si){
    Vector tmp = light-p;
    double d = norm(tmp);
    tmp = tmp/d;
    Ray beam(p, tmp); 
    sphereIpointP res = intersectScene(beam, si);

    if (res.i==-1) return 1;
    if(norm(light-p)>norm(res.inter-p)) return 0;
    return 1;
}

Vector normalTriangle(int i, int j){
    TriangleIndices tmp = scene[i].indices[j];
    Vector A = scene[i].vertices[tmp.vtxi];
    Vector B = scene[i].vertices[tmp.vtxj];
    Vector C = scene[i].vertices[tmp.vtxk];
    return normalize(cross(A-B, A-C));
}

Vector colorTriangle(int i, int j){
    TriangleIndices tmp = scene[i].indices[j];
    return Vector(255,255,255);
}

Vector lambertian(Vector p, int si, int sj, Light light){
    Vector n = normalTriangle(si, sj);
    Vector tmp = light.p-p;
    double d = norm(tmp);
    tmp = (light.I/(4.*PI*PI*d*d) * visibility(p, light.p, si) * max(dot(n, tmp/d), 0.)) * colorTriangle(si,sj);
    return min(tmp, Vector(255,255,255));
}
/*
Vector mirorSurface(Vector p, int si, Light light, int depth, Vector previous){
    Vector omegaI = normalize(p-previous);
    Vector n = normalize(p - scene[si].p);
    Vector omegaR = omegaI - 2 * dot(omegaI, n) * n;
    omegaR = normalize(omegaR);

    Ray ray(p, omegaR);
    sphereIpointP best = intersectAll(ray);

    if (best.i == -1) return NULLVEC;
    return getColor(best.inter, best.i, light, depth, p);
}
*/
Vector gammaCor(Vector color){
    double x = pow(color[0],1/gamma) * rgbCorrection;
    double y = pow(color[1],1/gamma) * rgbCorrection;
    double z = pow(color[2],1/gamma) * rgbCorrection;
    return Vector(x, y, z);
}
/*
Vector intersectSelf(Sphere s, Ray r){
    Vector tmp = r.p-s.p;
    double t = dot(r.d, tmp);
    double delta = t*t;
    delta -= (norm(tmp)*norm(tmp) - s.R*s.R);
    t = -t;

    if (delta<0 || t<0) return NULLVEC;

    if (delta != 0.)
        delta = sqrt(delta);
    t += delta; 
    if (norm(t*r.d) < s.R/1000.) return NULLVEC;
    return r.p + t*r.d;
}

Vector refract(Vector p, int si, Light light, int depth, Vector previous){
    double n1 = scene[si].c[1];
    double n2 = scene[si].c[2];

    Vector omegaI = normalize(p-previous);
    Vector n = normalize(p - scene[si].p);
    double tmpDot = dot(omegaI, n);

    double R = 1. - 4.*n1*n2/(n1*n1+2.*n1*n2+n2*n2);
    R = R + (1-R)*pow(1.-abs(tmpDot), 5.);

    if (random() < R){
        return mirorSurface(p, si, light, depth+1, previous);
    }
    
    if (tmpDot>0){
        n = -1*n;
        double a = n1;
        n1 = n2;
        n2 = a;
    }
    tmpDot = dot(omegaI, n);
    double n12 = n1/n2;

    double rad = 1-n12*n12*(1-tmpDot*tmpDot);
    if(rad<0){
        return mirorSurface(p, si, light, depth+1, previous);
    }

    Vector omegaT = n12*(omegaI-tmpDot*n);
    omegaT = omegaT - n*sqrt(rad);
    omegaT = normalize(omegaT);
    
    Ray ray(p, omegaT);
    sphereIpointP best = intersectAll(ray, si);
    Vector inter = intersectSelf(scene[si], ray);

    if (inter!=NULLVEC){
        if (best.i == -1) return getColor(inter, si, light, depth, p);
        double d1 = norm(p-best.inter);
        double d2 = norm(p-inter);
        if(d2 < d1){
            return getColor(inter, si, light, depth, p);
        }
    }
    
    if (best.i == -1)
        return NULLVEC;
    
    return getColor(best.inter, best.i, light, depth, p);
}
*/
Vector randomVec(const Vector& n){
    const double r1 = random();
    const double r2 = random();
    const double x = cos(2.*PI*r1)*sqrt(1.-r2);
    const double y = sin(2.*PI*r1)*sqrt(1.-r2);
    const double z = sqrt(r2);
    Vector T1;
    if(abs(n[0])<abs(n[1])){
        if(abs(n[0])<abs(n[2])){
            T1 = Vector(0, -n[2], n[1]);
        } else{
            T1 = Vector(-n[1], n[0], 0);
        }
    } else{
        if(abs(n[1])<abs(n[2])){
            T1 = Vector(-n[2], 0, n[0]);
        } else{
            T1 = Vector(-n[1], n[0], 0);
        }
    }
    T1 = normalize(T1);
    Vector T2 = normalize(cross(n,T1));
    return normalize(x*T1 + y*T2 + z*n);
}

Vector indirectLight(Vector p, int si, int sj, Light light, int depth){
    if (depth<=0) return NULLVEC;
    depth -= 1;
    Vector diffusion(0,0,0);
    Vector n = normalTriangle(si, sj);

    Vector ranedVec = randomVec(n);
    Ray ray(p, ranedVec);
    sphereIpointP tmpbest = intersectScene(ray);
    if (tmpbest.i != -1)                                
        diffusion += getColor(tmpbest.inter, tmpbest.i, tmpbest.j, light, depth, p);
    return diffusion;
}

void boxMuller(double& x, double& y, double stdev=0.3){
    const double r1 = random();
    const double r2 = random();
    x = x+sqrt(-2*log(r1))*cos(2*PI*r2)*stdev;
    y = y+sqrt(-2*log(r1))*sin(2*PI*r2)*stdev;
}

Vector sphericalLight(Vector p, int si, int sj, Light light){
    Vector n = normalTriangle(si, sj);
    Vector tmp = p-light.p;
    double d = norm(tmp);
    Vector nPrime = randomVec(tmp/d);
    Vector xPrime = light.p + light.R*nPrime;
    Vector omegaI = normalize(xPrime - p);
    double vis = visibility(p, xPrime, si);
    double R = light.R;
    double pdf = dot(nPrime, tmp/d)/PI/R/R;

    if (pdf==0) return NULLVEC; 
    tmp = light.I/(4.*PI*PI*R*R) * colorTriangle(si,sj)/PI * vis * max(dot(n, omegaI), 0.) * max(dot(nPrime, -1*omegaI), 0.) / (pow(norm(xPrime-p), 2)*pdf);
    return tmp;
}

Vector getColor(Vector p, int si, int sj, Light light, int depth, Vector previous){
    if (p==previous || depth<=0 || p!=p) return NULLVEC;
    depth -= 1;
    return (sphericalLight(p, si, sj, light) + colorTriangle(si,sj) * indirectLight(p, si, sj, light, depth-1) / 255) / 2;
/*    switch (scene[si].m){
    case opaque: // lambertian or sphericalLight
        return (sphericalLight(p, si, light) + scene[si].c * indirectLight(p, si, light, depth-1) / 255) / 2;
    case miror:
        return mirorSurface(p, si, light, depth, previous);
    case transparent:
        return refract(p, si, light, depth, previous);
    default:
        return NULLVEC;
    }*/
}

#endif