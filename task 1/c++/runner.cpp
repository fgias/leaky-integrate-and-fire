#include <iostream>
#include <fstream>
using namespace std;
// using namespace std::chrono;

int main() {
    int N = 2000;
    int R = 700;
    double sigma = 0.7;
    double mu = 1;
    double u_th = 0.98;
    double dt = 0.01;

    double u[N];
    double circles[N];

    for (int i = 1; i <= N; i++) {
        double r = rand() % 10;
        u[i] = r / 10.0;
        circles[i] = 0;
    }

    double t = 0;
    double time_max = 300;

    // auto start = high_resolution_clock::now();

    while (t < time_max) {
        for (int i=0; i<N; i++) {
            u[i] += dt*(mu - u[i]);
            double coupling = 0;
            if (i-R<0){
                for (int j=((N + ((i-R)%N)) % N); j<N; j++) {
                    coupling += sigma*(u[i]-u[j])/(2*R);                  
                }
                for (int j=0; j<i+R+1; j++) {
                    coupling += sigma*(u[i]-u[j])/(2*R); 
                }
                u[i] += dt*coupling;
            }
            else if (i+R > N-1) {
                for (int j=i-R; j<N; j++) {
                    coupling += sigma*(u[i]-u[j])/(2*R);  
                }
                for (int j=0; j<((N + ((i+R)%N)) % N + 1); j++) {
                    coupling += sigma*(u[i]-u[j])/(2*R);  
                }
                u[i] += dt*coupling;
            }
            else {
                for (int j=i-R; j<i+R+1; j++){
                    coupling += sigma*(u[i]-u[j])/(2*R); 
                }
                u[i] += dt*coupling;
            }

            if (u[i] > u_th) {
                u[i] = 0;
                circles[i] += 1;
            }
        }
        cout << t << endl;
        t += dt;
    }

    // auto stop = high_resolution_clock::now();

    // auto duration = duration_cast<microseconds>(stop - start);

    // cout << duration.count() << endl;


    // int a = -1;
    // int b = 5;

    // cout << (b + (a%b)) % b << endl;

    ofstream arrayData("u.txt"); // create file

    for(int k=0; k<N; k++)
    {
        arrayData << u[k] << endl; // output to txt
    }

    return 0;
}