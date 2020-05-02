#include "helpers.cpp"

int main(int argc, char *argv[]){
    int image[w*h*3];
        
    // camera
    Vector camPoint = Vector(0, 0, 55);
    cam = Camera(camPoint, 60.);

    // spheres
    buildScene();

    // lights
    
    Vector lightSource(-10, 20, 40);
    int lightI = 100000; //max 1000000000
    Light light(lightSource, lightI, 5);
    //when working with more lights, make dedicated class (point+I) and list of sources

    for (int i=0; i<w*h*3; i++){
        image[i]=0;
    }

    #pragma omp parallel for
    for(int i=0; i<w; i++){
        #pragma omp parallel for
        for(int j=0; j<h; j++){
            Vector color(0,0,0);
            int amount = 10;
            
            #pragma omp parallel for
            for (int k=0; k<amount; k++){
                double tmpI = double(i);
                double tmpJ = double(j);
                boxMuller(tmpI, tmpJ);
                Vector pixel = getPixCoord(tmpI,tmpJ);

                double D = norm(cam.p - pixel);
                D = cam.f*D/(D-cam.f);
                Vector u = normalize(pixel-cam.p);
                pixel = cam.p + D / abs(u[2]) * u;
                Camera localCam = cam;
                double r = random()*localCam.R;
                double omega = random()*2*PI;
                localCam.p = Vector(cam.p[0]+cos(omega)*r, cam.p[1]+sin(omega)*r, cam.p[2]);
                
                Vector direction = normalize(pixel-localCam.p);
                Ray ray(localCam.p, direction);
                sphereIpointP best = intersectAll(ray);

                if (best.i != -1){
                    color += getColor(best.inter, best.i, light, 10, localCam.p);
                }
            }
            
            color = color/amount;
            color = gammaCor(color);
            
            image[(j*w+i)*3+0] = min(max(int(color[0]),0), 255);
            image[(j*w+i)*3+1] = min(max(int(color[1]),0), 255);
            image[(j*w+i)*3+2] = min(max(int(color[2]),0), 255);                
        }
    }

    writePPM(w, h, image);
    print("Done");
    return 0;
}