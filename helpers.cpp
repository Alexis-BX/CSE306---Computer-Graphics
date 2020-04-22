#include "objects.cpp"

double PI{3.14159265};
const int w {300};
const int h {400};
const double gamma {2.2}; // usual: 2.2
const double rgbCorrection = pow(255, (gamma-1)/gamma);

Camera cam;
vector<Sphere> scene;

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

void addSphere(double x, double y, double z, double r, double g, double b, double R, Materials m=opaque){
    Vector tmpPoint(x,y,z);
    Vector tmpColor(r,g,b);
    scene.push_back(Sphere(tmpPoint, R, tmpColor, m));
}

void addHollowSphere(double x, double y, double z, double t, double n1, double n2, double R){
    addSphere(x, y, z, t, n1, n2, R-0.1, transparent);
    addSphere(x, y, z, t, n1, n2, R, transparent);
}

void buildScene(){
    addSphere(0, 0, 1000, 255, 0, 255, 940);
    addSphere(0, -1000, 0, 0, 0, 255, 990);
    addSphere(0, 1000, 0, 255, 0, 0, 940);
    addSphere(0, 0, -1000, 0, 255, 0, 940);
    addSphere(0, 0, 0, 100, 100, 100, 10);
    addSphere(-20, 0, 0, 100, 100, 100, 10, miror);

    //Hollow sphere
    //for transparents color: (transparency (0 opaque, 1 transparent), n1, n2)
    addSphere(20, 0, 0, 1, 1, 2, 10, transparent);
}

Vector getPixCoord(int i, int j){
    int x = i;
    int y = h-j-1;
    return Vector(cam.p[0]+x+0.5-w/2., cam.p[1]+y+0.5-h/2., cam.p[2]-w/(2.*tan(cam.fov/2.*PI/180.)));
}

Vector intersect(Sphere s, Ray r){
    Vector tmp = r.p-s.p;
    double t = dot(r.d, tmp);
    double delta = t*t;
    delta -= (norm(tmp)*norm(tmp) - s.R*s.R);
    t = -t;

    Vector nothing(0,0,0);
    if (delta<0 || t<0) return nothing;

    if (delta == 0.){
        return r.p + t*r.d;
    }

    delta = sqrt(delta);
    if (t<delta){
        t += delta;
    }
    else{
        t -= delta;
    }

    return r.p + t*r.d;
}

struct sphereIpointP{
    int i = -1;
    Vector inter;
};

sphereIpointP intersectScene(Ray ray, int skip=-1){
    sphereIpointP res;

    double closest = 0;

    for(int k=0; k<int(scene.size()); k++){
        if (k==skip) continue;
        Sphere& tmpSphere = scene[k];
        Vector inter = intersect(tmpSphere, ray);
        if (inter[0]!=0. || inter[1]!=0. || inter[2]!=0.){
            Vector tmpInter = cam.p-inter;
            double d = norm(tmpInter);
            if(d < closest || res.i==-1){
                res.i = k;
                res.inter = inter;
                closest = d;
            }

        }
    }
    return res;
}

Vector normalSatP(Vector p, Sphere s){
    Vector v = p-s.p;
    return v/norm(v);
}

double visibility(Vector p, Vector light, int si){
    Vector tmp = light-p;
    double d = norm(tmp);
    tmp = tmp/d;
    Ray beam(p, tmp); 
    double passThrough = 1;
    for(int i=0; i<int(scene.size()); i++){
        if (i==si) {
            continue;
        }
        Vector inter = intersect(scene[i], beam);
        if (inter[0]!=0. || inter[1]!=0. || inter[2]!=0.){
            inter = p-inter;
            if(norm(inter) < d){
                switch (scene[i].m) {
                case transparent:
                    passThrough = scene[i].c[0];
                case opaque:    
                default:
                    passThrough = 0.;
                    break;
                }
            }
        }
    }
    return passThrough;
}

Vector lambertian(Vector p, int si, Vector light, int I){
    //point, sphere index, light, intensity, visibility of light
    // p is not ofset as we ignore the sphere p belongs to when looking for shadows
    Vector n = normalSatP(p, scene[si]);
    Vector tmp = light-p;
    double d = norm(tmp);
    return (I/(4.*PI*PI*d*d) * visibility(p, light, si) * max(dot(n, tmp/d), 0.)) * scene[si].c;
}

Vector getColor(Vector p, int si, Vector light, int I, int depth=10, Vector previous=cam.p);

Vector mirorSurface(Vector p, int si, Vector light, int I, int depth, Vector previous){
    if (depth<=0) return Vector(0,0,0);
    depth -= 1;

    Vector omegaI = p-previous;
    omegaI = omegaI/norm(omegaI);
    Vector n = normalSatP(p, scene[si]);
    Vector omegaR = omegaI - 2 * dot(omegaI, n) * n;
    omegaR = omegaR/norm(omegaR);

    Ray ray(p, omegaR);
    sphereIpointP best = intersectScene(ray);

    if (best.i == -1)
        return Vector(0,0,0);
    else
        return getColor(best.inter, best.i, light, I, depth, p);
}

Vector gammaCor(Vector color){
    double x = pow(color[0],1/gamma) * rgbCorrection;
    double y = pow(color[1],1/gamma) * rgbCorrection;
    double z = pow(color[2],1/gamma) * rgbCorrection;
    return Vector(x, y, z);
}

Vector intersectSelf(Sphere s, Ray r){
    Vector tmp = r.p-s.p;
    double t = dot(r.d, tmp);
    double delta = t*t;
    delta -= (norm(tmp)*norm(tmp) - s.R*s.R);
    t = -t;

    Vector nothing(0,0,0);
    if (delta<0 || t<0) return nothing;

    if (delta != 0.) delta = sqrt(delta);

    t += delta;
    return r.p + t*r.d;
}

Vector refract(Vector p, int si, Vector light, int I, int depth, Vector previous){
    if (depth<=0) return Vector(0,0,0);
    depth -= 1;

    double n1 = scene[si].c[1];
    double n2 = scene[si].c[2];

    Vector omegaI = p-previous;
    omegaI = omegaI/norm(omegaI);
    Vector n = normalSatP(p, scene[si]);
    double tmpDot = dot(omegaI, n);

    if (n2<n1) return mirorSurface(p, si, light, I, depth+1, previous);

    if(tmpDot>0.) n = -1*n;
    tmpDot = dot(omegaI, n);
    double n12 = n1/n2;

    double R = 1. - 4*n1*n2/(n1*n1+2*n1*n2+n2*n2);
    R = R + (1-R)*pow(1.-abs(tmpDot), 5.);
    
    if (double(rand()) / double(RAND_MAX) < R){
        return mirorSurface(p, si, light, I, depth+1, previous);
    }
    

    Vector omegaT = n12*(omegaI-tmpDot*n);
    omegaT = omegaT - n*sqrt(1-n12*n12*(1-tmpDot*tmpDot));
    omegaT = omegaT/norm(omegaT);

    Ray ray(p, omegaT);
    sphereIpointP best = intersectScene(ray, si);

    Vector inter = intersectSelf(scene[si], ray);
    if (inter[0]!=0. || inter[1]!=0. || inter[2]!=0.){
        Vector tmp = p-best.inter;
        double d1 = norm(tmp);
        tmp = p-inter;
        double d2 = norm(tmp);
        if(d2 < d1){
            return getColor(inter, si, light, I, depth, p);
        }
    }
    
    if (best.i == -1)
        return Vector(0,0,0);
    
    return getColor(best.inter, best.i, light, I, depth, p);
}

Vector randomVec(const Vector& n){
    const double r1 = double(rand()) / double(RAND_MAX);
    const double r2 = double(rand()) / double(RAND_MAX);
    const double x = cos(2.*PI*r1)*sqrt(1-r2);
    const double y = sin(2.*PI*r1)*sqrt(1-r2);
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
    T1 = T1 / norm(T1);
    Vector T2 = cross(n,T1);
    //T2 = T2/norm(T2);
    //Vector tmp = x*T1 + y*T2 + z*n;
    //tmp = tmp/norm(tmp);
    return x*T1 + y*T2 + z*n;
}

Vector indirectLight(Vector p, int si, Vector light, int I, int depth){
    if (depth<=0) return Vector(0,0,0);
    depth -= 1;
    Vector diffusion(0,0,0);
    int amount = 10;
    Vector n = normalSatP(p, scene[si]);

    for (int k=0; k<amount; k++){
        Vector ranedVec = randomVec(n);
        Ray ray(p, ranedVec);
        sphereIpointP tmpbest = intersectScene(ray);
        if (tmpbest.i != -1)                                
            diffusion += getColor(tmpbest.inter, tmpbest.i, light, I, depth);
    }
    return diffusion/amount;
}

Vector getColor(Vector p, int si, Vector light, int I, int depth, Vector previous){
    switch (scene[si].m){
    case opaque:
        return lambertian(p, si, light, I);
    case miror:
        return mirorSurface(p, si, light, I, depth, previous);
    case transparent:
        return refract(p, si, light, I, depth, previous);
    default:
        return Vector(0,0,0);
    }
}

















