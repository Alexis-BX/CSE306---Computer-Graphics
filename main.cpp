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
    int lightI = 80000; //max 1000000000
    //when working with more lights, make dedicated class (point+I) and list of sources

    double rgbCorrection = pow(255, (gamma-1)/gamma);

    for (int i=0; i<w*h*3; i++){
        image[i]=0;
    }

    for(int i=0; i<w; i++){
        for(int j=0; j<h; j++){
            Vector pixel = getPixCoord(i,j);
            Vector direction = pixel-cam.p;
            direction = direction/norm(direction);
            Ray ray(cam.p, direction);

            sphereIpointP best = intersectScene(ray);

            if (best.i == -1){
                image[(j*w+i)*3+0] = 0;
                image[(j*w+i)*3+1] = 0;
                image[(j*w+i)*3+2] = 0;
            } else {
                Vector color = getColor(best.inter, best.i, lightSource, lightI, 2);
                image[(j*w+i)*3+0] = int(min(pow(color[0],1/gamma) * rgbCorrection, 255.));
                image[(j*w+i)*3+1] = int(min(pow(color[1],1/gamma) * rgbCorrection, 255.));
                image[(j*w+i)*3+2] = int(min(pow(color[2],1/gamma) * rgbCorrection, 255.));                
            }

            
        }
    }

    writePPM(w, h, image);
    std::cout << "Done" << std::endl;
    return 0;
}








































