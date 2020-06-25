#ifndef _FILEIO_HPP_
#define _FILEIO_HPP_

#include "typedefs.hpp"
#include "Vector3.hpp"
#include "Matrix3x3.hpp"
#include <vector>
#include <iostream>
#include <fstream>
#include <string>

using namespace std;

string fname = "positions_";

int N = 2000;
int fn;

void writeToFile(int n_points, bool allowedToBreak) {
    ofstream location_out;
    for (int i=0; i<N; i++) {
        string file = fname + to_string(n_points) + "_" + to_string(i) + ".txt";
        ifstream fin(file);
        if(!fin) {
            fn = i;
            break;
        }
    }
    location_out.open(fname + to_string(n_points) + "_" + to_string(fn) + ".txt", std::ios::out | std::ios::app);
    if (!location_out.is_open()) {
        cout << "Cannot open the file." << endl;
        return;
    }
    location_out << to_string(n_points) << endl;
    location_out << to_string(allowedToBreak) << endl;
    location_out.close();
}

void writeToFile(Real t, Real dt, int n_points) {
    ofstream location_out;
    location_out.open(fname + to_string(n_points) + "_" + to_string(fn) + ".txt", std::ios::out | std::ios::app);
    if (!location_out.is_open()) {
        cout << "Cannot open the file." << endl;
        return;
    }
    location_out << to_string(t) << " " << to_string(dt) << endl;
    location_out.close();
}

void writeToFile(vector<Vec3f> pos, int n_points) {
    ofstream location_out;
    location_out.open(fname + to_string(n_points) + "_" + to_string(fn) + ".txt", std::ios::out | std::ios::app);
    if (!location_out.is_open()) {
        cout << "Cannot open the file." << endl;
        return;
    }
    for (int i=0; i<n_points; i++) {
        location_out << to_string(pos[i][0]) << " " << to_string(pos[i][1]) << " " << to_string(pos[i][2]) << endl;
    }
    location_out.close();
}




#endif  /* _CALCULATIONS_HPP_ */
