#include "../vector.cpp"

class Polygon{
public:
    std::vector<Vector> p;
    Polygon(){}
    Polygon(std::vector<Vector> points){
        p = points;
    }
    Vector& operator[](int i) {return p[i];}  
    const Vector& operator[](int i) const {return p[i];}
    const unsigned int size() {return p.size();}
    void add(Vector tmp) {p.push_back(tmp);}
};

std::vector<Polygon> vple(std::vector<Vector> points){
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
    Polygon boundingBox(tmp);

    const double e = 10e-6;
    std::vector<Polygon> set(points.size());
    
    #pragma omp parallel for
    for(int ic = 0; ic < points.size(); ic++) {
        Vector center = points[ic];
        Polygon cell = boundingBox;

        for(int ip = 0; ip < points.size(); ip++) {
            if(ic == ip) continue;

            Vector point = points[ip];

            Vector N = normalize(point - center);
            double t = dot(center + ((point - center) / 2), N);

            int n = cell.size();
            std::vector<double> ts(n);
            for(int i = 0; i < n; i ++) ts[i] = dot(N, cell[i]) - t;

            Polygon ncell;
            for(int i = 0; i < n; i ++) {
                int j = (i+1)%n;
                int inter = 0;
                if(ts[i] < e) {
                    ncell.add(cell[i]);
                    if(ts[j] > -e) 
                        inter = 1;
                } else if(ts[j] < e) 
                    inter = -1;

                if(inter) {
                    double ijN = dot((cell[j] - cell[i]), N);
                    ncell.add((std::abs(ijN) < e)? cell[(inter==1)?j:i] : cell[i]+(cell[j]-cell[i])*(t-dot(cell[i],N))/ijN);
                }
            }
            for(int i = 0; i < ncell.size(); i ++) 
                ncell[i][2] = 0;
            cell = ncell;
            
        }
        set[ic] = cell;
    }
    return set;

}