#define simn 10000000
#define A 20.0
#include<iostream>
#include<cmath>
#include <complex>
#include <random>
#define cpx complex <long double>
#define ldb long double
#define pi 3.14159265358979
#define fc 10
#define simp 100
#define simt 1.0/simp
#define N 2.0
using namespace std;
ldb sinwv[simp],coswv[simp];
cpx cir[8];
int mp[8]={0,1,3,2,7,6,4,5},iv[8]={0,1,3,2,6,7,5,4};
void genwv(){
    for(int i =0;i<simp;++i){
        sinwv[i] = sin(2*pi*simt*fc*i);
        coswv[i] = cos(2*pi*simt*fc*i);
    }
    for(int i=0;i<8;++i){
        cir[i] = cpx(cos(pi*i/4),sin(pi*i/4));
    }
}
int main(){
    srand(time(0));
    genwv();
    default_random_engine rgen;
    normal_distribution<double> dist(0,N);
    int n,scr,b0,b1,b2,be,out,osy;
    ldb mod[simp],sumr,sumi,nis1,nis2,rp,ip,pha;
    scr = b0 = b1 = b2 = be = 0;
    for(int k=0;k<simn;++k){
        sumr = sumi = 0;
        n = rand()%8;
        rp = A*real(cir[mp[n]]);
        ip = A*imag(cir[mp[n]]);
        nis1 = dist(rgen);
        nis2 = dist(rgen);
        for(int i=0;i<simp;++i){
            mod[i] = (rp+nis1)*coswv[i]+(ip+nis2)*sinwv[i];
        }
        for(int i=0;i<simp;++i){
            sumr += mod[i]*coswv[i];
            sumi += mod[i]*sinwv[i];
        }
        pha = arg(cpx(sumr,sumi));
        pha += (pha<-pi/8)?2*pi:0;
        osy = int((pha+pi/8)/(pi/4));
        out = iv[osy];
        scr += n != out;
        b0 += (((n^out)&4)==0)?0:1;
        b1 += (((n^out)&2)==0)?0:1;
        b2 += (((n^out)&1)==0)?0:1;
    }
    cout<<scr<<" "<<b0<<" "<<b1<<" "<<b2<<" "<<simn<<endl;
    return 0;
}
