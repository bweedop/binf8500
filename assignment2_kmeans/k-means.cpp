#include <iostream>
#include <vector>
#include <fstream>
#include <chrono>
#include <random>

using namespace std;

const double XMAX = 25.0;
const double YMAX = 50.0;
const int NCLUST = 3;

struct Point {
    int label;
    double x, y, dist;
};

struct Centroid {
    int label;
    double x, y;
};

void averages (vector<Point> data, double &x_avg, double &y_avg) {
    double x_sum, y_sum = 0.0;
    double size = data.size();

    if (size != 0) {
        for (int i = 0; i < data.size(); i++) {
            x_sum = x_sum + data[i].x;
            y_sum = y_sum + data[i].y;
        }
        x_avg = x_sum / size;
        y_avg = y_sum / size;
    } else {
        throw domain_error("One of your clusters has no points");
    }
}

double distance (const Point pnt1, const Centroid pnt2) {
    double distance = sqrt(pow((pnt2.x-pnt1.x), 2.0) + pow((pnt2.y-pnt1.y), 2.0));
    return distance;
}

void all_distances (vector<Point> &pnts, vector<Centroid> &clst) {
    for (int i = 0; i < NCLUST; i++) {
        for (int j = 0; j < pnts.size(); j++) {
            double tmp_dist;
            tmp_dist = distance(pnts[j], clst[i]);
            if (tmp_dist < pnts[j].dist) {
                pnts[j].dist = tmp_dist;
                pnts[j].label = clst[i].label;
            }
        }
    }
}

double random_number (double upper, double lower) {
    /* initialize random seed: */
    unsigned seed = chrono::system_clock::now().time_since_epoch().count();
    default_random_engine generator (seed);
    uniform_real_distribution<double> unif(lower, upper);
    double a_random_double = unif(generator);
    return a_random_double;
}

int main (int argc, char* argv[]) {
    
    vector<Point> points;
    vector<Point> pnt_clst[NCLUST];
    vector<Centroid> clusters;
    double avg_x, avg_y;

    ifstream data_file (argv[1]);
    if (data_file.is_open()) {
        while (!data_file.eof()) {
            Point tmp;
            data_file >> tmp.x >> tmp.y;
            tmp.dist = 10000.0; // Set it to something initially 
            points.push_back(tmp);
        }
        data_file.close();
    } else {
        throw domain_error("Could not read data file...");
    }

    for (int i = 0; i < NCLUST; i++) {
        Centroid tmp;
        tmp.label = i;
        tmp.x = random_number(XMAX, 0.0);
        tmp.y = random_number(YMAX, 0.0);
        clusters.push_back(tmp);
    }

    for (int i = 0; i < NCLUST; i++)
        cout << clusters[i].x << " " << clusters[i].y << endl;

    bool change = true;
    while (change) {
        vector<Centroid> tmp = clusters;
        int check = 0;
        all_distances(points, clusters);

        for (int j = 0; j < points.size(); j++) {
            int tmp = points[j].label;
            pnt_clst[tmp].push_back(points[j]);
        }

        for (int i = 0; i < NCLUST; i++) {
            averages(pnt_clst[i], clusters[i].x, clusters[i].y);
        }

        for (int i = 0; i < NCLUST; i++) {
            cout << clusters[i].x << " " << clusters[i].y << endl;
            if (round(clusters[i].x)/100 == round(tmp[i].x)/100 && round(clusters[i].y)/100 == round(tmp[i].y)/100) {
                check = check + 1;
            }
        }
        cout << check << endl;
        if (check == NCLUST) {
            change = false;
        }        
    }
    
    ofstream output ("k-output.txt");
    if (output.is_open()) {
        for (int j = 0; j < points.size(); j++) {
            output << points[j].label << " " << points[j].x << " " << points[j].y << endl;
    }
    }
    output.close();
    return 0;
}