#include<iostream>
#include<cmath>
#include <random>
#define simn 10000000
#define A 6.0
#define ldb long double
#define pi 3.14159265358979
#define fc 10
#define simp 100
#define simt 1.0/simp
#define N 2.0
using namespace std;
ldb sinwv[simp],coswv[simp];
void genwv(){
    for(int i =0;i<simp;++i){
        sinwv[i] = sin(2*pi*simt*fc*i);
        coswv[i] = cos(2*pi*simt*fc*i);
    }
}
int main(){
    genwv();
    default_random_engine rgen;
    normal_distribution<double> dist(0,N);
    int n,crct=0;
    ldb mod[simp],sum,nois;
    for(int k=0;k<simn;++k){
        sum = 0;
        n = rand()%2?1:-1;
        nois = dist(rgen);
        for(int i=0;i<simp;++i)
            mod[i] = n*(A+nois)*coswv[i];
        for(int i=0;i<simp;++i)
            sum += mod[i]*coswv[i];
        crct+=sum*n<0;
    }
    cout<<crct<<" "<<simn<<endl;
    return 0;
}
