#include "svg.cpp"

double rand01(){
    return double(rand()) / double(RAND_MAX);
}

int main(int argc, char *argv[]){
    
    int N = 3;
    std::vector<Vector> points(N);
    for(int i = 0; i < N; i ++) {
        points[i] = Vector(rand01(), rand01(), 0);
    }
    
    std::vector<Polygon> voronois = vple(points);
    print("Done");
    save_svg(voronois, "res.svg");
}