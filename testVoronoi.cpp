#include "svg.cpp"

double rand01(){
    return double(rand()) / double(RAND_MAX);
}

std::vector<Vector> tutte(std::vector<Vector> points, std::vector<std::vector<int>> adjcancy, std::vector<int> boundaryIndex) {
    int n = boundaryIndex.size();
    
    double s = 0;
    for (int i=0; i<n; i++)
        s += norm(points[boundaryIndex[i]] - points[boundaryIndex[(i+1)<n?i+1:0]]);
    
    double cs = 0;
    std::vector<Vector> v;

    for (int i=0; i<n; i++){
        double teta = 2*PI*cs/s;
        points[boundaryIndex[i]] = Vector(cos(teta), sin(teta), 0);
        cs += norm(points[boundaryIndex[i]] - points[boundaryIndex[(i+1)<n?i+1:0]]);
    }

    for (int j=0; j<100; j++){
        std::vector<Vector> tmp = points;
        for(int i = 0; i < n; i ++) {
            tmp[i] = NULLVEC;
            for(int k=0; k<adjcancy[i].size(); k++) 
                tmp[i] += points[adjcancy[i][k]];
            tmp[i] = tmp[i] / double(adjcancy[i].size());
        } 
        for(int i=0; i<n; i++) tmp[boundaryIndex[i]] = points[boundaryIndex[i]];
        points = tmp;
    }
    return points;
}

int main(int argc, char *argv[]){
    
    int N = 10;
    std::vector<Vector> points(N);
    for(int i = 0; i < N; i ++) {
        points[i] = Vector(rand01(), rand01(), 0);
    }
    
    std::vector<Polygon> voronois = vple(points);
    print("Done");
    save_svg(voronois, "res.svg");
}