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

    for (int i=0; i<w*h*3; i++){
        image[i]=0;
    }
    #pragma omp parallel for
    for(int i=0; i<w; i++){
        #pragma omp parallel for
        for(int j=0; j<h; j++){
            Vector pixel = getPixCoord(i,j);
            Vector direction = normalize(pixel-cam.p);
            Ray ray(cam.p, direction);

            sphereIpointP best = intersectScene(ray);

            Vector color(0,0,0);
            if (best.i != -1){    
                int amount = 20;
                #pragma omp parallel for
                for (int k=0; k<amount; k++){
                    color += getColor(best.inter, best.i, lightSource, lightI, 3);
                }
                color = color/amount;
                
                color = gammaCor(color);
            }
            
            image[(j*w+i)*3+0] = min(max(int(color[0]),0), 255);
            image[(j*w+i)*3+1] = min(max(int(color[1]),0), 255);
            image[(j*w+i)*3+2] = min(max(int(color[2]),0), 255);                
        }
    }

    writePPM(w, h, image);
    std::cout << "Done" << std::endl;
    return 0;
}








































