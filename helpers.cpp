#include "objects.cpp"

double PI{3.1415};
const int w {300};
const int h {400};
const double gamma {2.2}; // usual: 2.2

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

void buildScene(){
    Vector tmpPoint;
    Vector tmpColor;

    tmpPoint = Vector(0,0,1000);
    tmpColor = Vector(255,0,255);
    scene.push_back(Sphere(tmpPoint, 940, tmpColor));

    tmpPoint = Vector(0,-1000,0);
    tmpColor = Vector(0,0,255);
    scene.push_back(Sphere(tmpPoint, 990, tmpColor));

    tmpPoint = Vector(0,1000,0);
    tmpColor = Vector(255,0,0);
    scene.push_back(Sphere(tmpPoint, 940, tmpColor));

    tmpPoint = Vector(0,0,-1000);
    tmpColor = Vector(0,255,0);
    scene.push_back(Sphere(tmpPoint, 940, tmpColor));

    tmpPoint = Vector(0,0,0);
    tmpColor = Vector(100,100,100);
    scene.push_back(Sphere(tmpPoint, 10, tmpColor));

    tmpPoint = Vector(-20,0,0);
    tmpColor = Vector(100,100,100);
    scene.push_back(Sphere(tmpPoint, 10, tmpColor, miror));

    tmpPoint = Vector(20,0,0);
    tmpColor = Vector(128,128,128);
    scene.push_back(Sphere(tmpPoint, 10, tmpColor, transparent));
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

sphereIpointP intersectScene(Ray ray){
    sphereIpointP res;

    double closest = 0;

    for(int k=0; k<int(scene.size()); k++){
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
                if(scene[i].m==opaque){
                    passThrough = 0;
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

Vector intersectSelf(Sphere s, Ray r){
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
    t += delta;

    return r.p + t*r.d;
}

Vector refract(Vector p, int si, Vector light, int I, int depth, Vector previous){
    if (depth<=0) return Vector(0,0,0);
    depth -= 1;
    
    Vector omegaI = p-previous;
    omegaI = omegaI/norm(omegaI);
    Vector n = normalSatP(p, scene[si]);
    double tmpDot = dot(omegaI, n);
    double n1, n2, nSphere{10};
    if(tmpDot>0.){
        n1 = 1;
        n2 = nSphere;
        n = -1*n;
    }
    else{
        n1 = 1;
        n2 = nSphere;
    }

    if (n2<n1){
        return mirorSurface(p, si, light, I, depth+1, previous);
    }
    double n12 = n1/n2;
    tmpDot = dot(omegaI, n);

    Vector omegaT = n12*(omegaI-tmpDot*n);
    omegaT = omegaT - n*sqrt(1-n12*n12*(1-tmpDot*tmpDot));
    omegaT = omegaT/norm(omegaT);

    Ray ray(p, omegaT);
    sphereIpointP best = intersectScene(ray);

    Vector inter = intersectSelf(scene[si], ray);
    if (inter[0]!=0. || inter[1]!=0. || inter[2]!=0.){
        return getColor(inter, si, light, I, depth, p);
        /*
        Vector tmp = p-best.inter;
        double d1 = norm(tmp);
        tmp = p-inter;
        double d2 = norm(tmp);
        if(d2 < d1){
            return getColor(inter, si, light, I, depth, p);
        }*/
    }
    
    if (best.i == -1)
        return Vector(0,0,0);
    
    return getColor(best.inter, best.i, light, I, depth, p);
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

















