#include<iostream>
#include<cmath>
#include <random>
#define simn 10
#define A 4.0
#define ldb long double
#define pi 3.14159265358979
#define fc 10
#define simp 100
#define simt 1.0/simp
#define N 2.0
using namespace std;
ldb sinwv[simp],coswv[simp];
int tk[10] = {1,2,4,8,16,32,64,128,256,512};
void genwv(){
    for(int i =0;i<simp;++i){
        sinwv[i] = sin(2*pi*simt*fc*i);
        coswv[i] = cos(2*pi*simt*fc*i);
    }
}
struct state{
    static bool g1[9],g2[9];
    static int rgm,cwl,gel;
    int reg,bit,reo;
    bool c1,c2;
    state(int r=0){
        reo = reg = r;
    };
    void cnv(int b){
        bit = b;
        //cout<<b<<" ";
        reg = (reg<<1)+b;
        c1=0,c2=0;
        for(int i=0;i<3;++i){
            c1 = c1 ^ (((reg&tk[i])>0)&g1[i]);
            c2 = c2 ^ (((reg&tk[i])>0)&g2[i]);
        }
        //cout<<c1<<c2<<endl;
    };
    int show(bool& g1,bool& g2){
        g1 = c1, g2 = c2;
        //cout<<bit<<" "<<c1<<c2<<endl;
        return reg&3;
    }
};
int state::gel=3;
int state::rgm=2;
int state::cwl=2;
bool state::g1[9]={1,0,1};
bool state::g2[9]={1,1,1};

int main(){
    genwv();
    default_random_engine rgen;
    normal_distribution<double> dist(0,N);
    int n,rb,rg=0,cw[240];
    bool r1=0,r2=0,g1,g2;
    ldb mod[simp],sum,nois,biterr=0;
    int arr[10]={1,0,1,1,0,1,1,1,0,1};
    state* map[4][2];
    state s1(0);
    for(int i=0;i<4;++i){
        for(int j=0;j<2;++j){
            map[i][j] = new state(i);
            map[i][j]->cnv(j);
        }
    }
    for(int k=0;k<simn;++k){
        rb = arr[k]; //rand();
        //s1.cnv(rb);
        //n = rb%2?1:-1;
        //g1 = rb^r2;
        //g2 = rb^r1^r2;
        //cout<<rb<<" "<<g1<<g2<<endl;
        rg = map[rg][rb]->show(g1,g2);
        //r2 = r1;
        //r1 = rb;
        /*{
            sum = 0;
            nois = dist(rgen);
            for(int i=0;i<simp;++i)
                mod[i] = n*(A+nois)*coswv[i];
            for(int i=0;i<simp;++i)
                sum += mod[i]*coswv[i];
            biterr+=sum*n<0;
        }*/
    }
    //cout<<"SNR: "<<A*A/2/N<<", BER: "<<biterr/simn;
    //cout<<" with "<<simn<<" simulations"<<endl;
    return 0;
}
