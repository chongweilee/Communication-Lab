//#define A 10.0
#define simn 100000
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
#define flr(x) (int(x)<x?int(x):int(x)-1)
using namespace std;
ldb sinwv[simp],coswv[simp];
cpx cir[16],z;
int mp[16]={0,4,12,8,1,5,13,9,3,7,15,11,2,6,14,10};
int iv[16]={0,4,12,8,1,5,13,9,3,7,15,11,2,6,14,10};
void genwv(){
    for(int i =0;i<simp;++i){
        sinwv[i] = sin(2*pi*simt*fc*i);
        coswv[i] = cos(2*pi*simt*fc*i);
    }
    for(int i=0;i<16;++i)
        cir[i] = cpx(-3,-3)+cpx(2*(i%4),2*int(i/4));
}
int main(){
    srand(time(0));
    genwv();
    default_random_engine rgen;
    normal_distribution<double> dist(0,N);
    double A;
    for(int m =1;m<20;++m){
        A = 0.5*m;
        int n,scr,b0,b1,b2,b3,be,our,oui,out,osy;
        ldb mod[simp],sumr,sumi,nis1,nis2,rp,ip;
        scr = b0 = b1 = b2 = b3 = be = 0;
        for(int k=0;k<simn;++k){
            sumr = sumi = 0;
            n = rand()%16;
            rp = A*real(cir[mp[n]]);
            ip = A*imag(cir[mp[n]]);
            nis1 = dist(rgen);
            nis2 = dist(rgen);
            for(int i=0;i<simp;++i)
                mod[i] = (rp+nis1)*coswv[i]+(ip+nis2)*sinwv[i];
            for(int i=0;i<simp;++i){
                sumr += mod[i]*coswv[i];
                sumi += mod[i]*sinwv[i];
            }
            z = (cpx(sumr*2/simp/A,sumi*2/simp/A));
            our = flr(real(z));
            oui = flr(imag(z));
            if(our<-3 || our>3) our = (our>0)?3:-3;
            if(oui<-3 || oui>3) oui = (oui>0)?3:-3;
            osy = int((oui+4)/2)*4+(our+4)/2;
            out = iv[osy];
            scr += out != n;
            b0 += (((n^out)&8)==0)?0:1;
            b1 += (((n^out)&4)==0)?0:1;
            b2 += (((n^out)&2)==0)?0:1;
            b3 += (((n^out)&1)==0)?0:1;
        }
        cout<<A<<"\t"<<scr<<"\t";
        cout<<b0<<"\t"<<b1<<"\t"<<b2<<"\t"<<b3<<"\t"<<simn<<endl;
    }
    return 0;
}
