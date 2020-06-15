#include "vector.cpp"

class Polygon{
public:
    std::vector<Vector> p;
    Polygon(){}
    Polygon(std::vector<Vector> points){
        p = points;
    }
    Vector& operator[](int i) {return (i>0)?p[i%p.size()]:p[p.size()+i%p.size()];}  
    const Vector& operator[](int i) const {return (i>0)?p[i%p.size()]:p[p.size()+i%p.size()];}
    const unsigned int size() {return p.size();}
};

std::vector<Polygon> vple(std::vector<Vector> points){
    std::vector<Polygon> set(points.size());
    if (points.size()==1) {
        set[0] = Polygon(points);
        return set;
    }

    Vector boxCornerMin = points[0];
	Vector boxCornerMax = points[0];
    for(int i=1; i<int(points.size()); i++){
        boxCornerMin = min(boxCornerMin, points[i]);
	    boxCornerMax = max(boxCornerMax, points[i]);
    }
    std::vector<Vector> tmp = {boxCornerMin, 
                               Vector(boxCornerMin[0], boxCornerMax[1], 0), 
                               boxCornerMax,
                               Vector(boxCornerMax[0], boxCornerMin[1], 0)};
    tmp = {Vector(0, 0, 0), Vector(1, 0, 0), Vector(1, 1, 0), Vector(0, 1, 0)};

    Polygon boundingBox(tmp);
    double e = 10e-6;
    
    #pragma omp parallel for
    for(int i = 0; i < int(points.size()); i++) {
        Vector center = points[i];
        Polygon cell = boundingBox.p;
        
        for(int j = 0; j < int(points.size()); j++) {
            if(i == j) continue;

            Vector N = normalize(points[j] - center);
            double t = dot(center + ((points[j] - center) / 2), N);

            int n = cell.size();
            std::vector<double> ts(n);
            for(int k = 0; k < n; k ++) ts[k] = dot(N, cell[k]) - t;
            Polygon ncell;
            for(int k = 0; k < n; k ++) {
                int l = (k+1)%n;
                int inter = 0;
                if(ts[k] < e) {
                    ncell.p.push_back(cell[k]);
                    if(ts[l] > -e) inter = 1;
                } else if(ts[l] < e) { inter = -1;}
                    
                if(inter!=0) {
                    Vector v = cell[l] - cell[k];
                    double abN = dot(v, N);
                    Vector x;
                    if(std::abs(abN) < e) {
                        if(inter == 1) x = cell[l];
                        else x = cell[k];
                    } else {
                        x = cell[k] + v * (t - dot(cell[k], N)) / abN;
                    }
                    ncell.p.push_back(x);
                }
            }
            cell = ncell;
            for(int k = 0; k < cell.size(); k ++) cell[k][2] = 0;
        }
        set[i] = cell;
    }
    return set;
    ///////////////////////////
    /*std::vector<Vector> tmp = {Vector(0, 0, 0), Vector(1, 0, 0), Vector(1, 1, 0), Vector(0, 1, 0)};

    Polygon boundingBox(tmp);
    double M = 0;
    double EPS = 10e-6;
    for(auto point : points) M = std::max(point[2], M);
    for(auto point : points) point[2] = sqrt(M - point[2]);
    std::vector<Polygon> voronois(points.size());

    for(int iC = 0; iC < int(points.size()); iC++) {
        Vector C = points[iC];
        std::vector<Vector> cell = boundingBox.p;

        for(int iP = 0; iP < int(points.size()); iP++) {
            if(iP == iC) continue;

            Vector P = points[iP];
            Vector CP = P - C;

            Vector N = normalize(CP);
            double T = dot(C + (CP / 2), N);

            int n = cell.size();
            std::vector<double> ts(n);
            for(int i = 0; i < n; i ++) ts[i] = dot(N, cell[i]) - T;

            std::vector<Vector> ncell;
            for(int i = 0; i < n; i ++) {
                int j = (i+1) % n;
                int inter = 0;
                if(ts[i] < EPS) {
                    ncell.push_back(cell[i]);
                    if(ts[j] > -EPS) inter = 1;
                } else if(ts[j] < EPS) inter = 2;
                if(inter) {
                    Vector a = cell[i];
                    Vector b = cell[j];
                    Vector ab = b - a;
                    double abN = dot(ab, N);
                    Vector x;
                    if(std::abs(abN) < EPS) {
                        if(inter == 1) x = b;
                        else x = a;
                    } else {
                        double t = (T - dot(a, N)) / abN;
                        x = a + ab * t;
                    }
                    ncell.push_back(x);
                }
            }
            cell = ncell;
            for(int i = 0; i < n; i ++) cell[i][2] = 0;
        }
        voronois[iC].p = cell;
    }

    return voronois;
*/
}