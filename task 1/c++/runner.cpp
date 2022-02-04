#include <iostream>
#include <stdlib.h>
#include <list>
using namespace std;

int main() {
    int N = 20;
    int R = 7;
    float sigma = 0.7;
    float mu = 1;
    float u_th = 0.98;
    float dt = 0.01;

    list<float> u;
    list<float> circles;

    for (int i = 1; i <= N; i++) {
        float r = rand() % 10;
        u.push_back(r / 10.0);
        circles.push_back(0.0);
    }

    float t = 0;
    float time_max = 10;

    return 0;
}