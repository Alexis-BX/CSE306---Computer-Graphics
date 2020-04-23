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
            Vector direction = pixel-cam.p;
            direction = direction/norm(direction);
            Ray ray(cam.p, direction);

            sphereIpointP best = intersectScene(ray);

            if (best.i == -1){
                image[(j*w+i)*3+0] = 0;
                image[(j*w+i)*3+1] = 0;
                image[(j*w+i)*3+2] = 0;
            } else {
                Vector color(0,0,0);
                switch (scene[best.i].m){
                case opaque:{
                    color = getColor(best.inter, best.i, lightSource, lightI, 5);
                    color = color + scene[best.i].c * indirectLight(best.inter, best.i, lightSource, lightI, 4);
                    break;
                    }
                case miror:{
                    color = getColor(best.inter, best.i, lightSource, lightI, 5);
                    break;
                    }
                case transparent:{
                    int amount = 10;
                    for (int k=0; k<amount; k++){
                        color += getColor(best.inter, best.i, lightSource, lightI, 10);
                    }
                    color = color/amount;
                    break;
                    }
                default:
                    break;
                }
                
                color = gammaCor(color);

                image[(j*w+i)*3+0] = minMaxInt(color[0]);
                image[(j*w+i)*3+1] = minMaxInt(color[1]);
                image[(j*w+i)*3+2] = minMaxInt(color[2]);                
            }

            
        }
    }

    writePPM(w, h, image);
    std::cout << "Done" << std::endl;
    return 0;
}








































