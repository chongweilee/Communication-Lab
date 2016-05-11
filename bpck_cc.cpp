#include<iostream>
#include<cmath>
#include<vector>
#include<random>
#define simn 120
#define cycl 50000
#define A 7.5
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
    int reg,reo;
    bool c1,c2;
    state(int r, int b){
        reo = reg = r;
        reg = (reg<<1)+b;
        c1=0,c2=0;
        for(int i=0;i<3;++i){
            c1 = c1 ^ (((reg&tk[i])>0)&g1[i]);
            c2 = c2 ^ (((reg&tk[i])>0)&g2[i]);
        }
    };
    int show(bool& g1,bool& g2){
        g1 = c1, g2 = c2;
        return reg&3;
    }
};
bool state::g1[9]={1,0,1};
bool state::g2[9]={1,1,1};

int main(){
    srand(time(0));
    genwv();
    default_random_engine rgen;
    normal_distribution<double> dist(0,N);
    int n,rb,rg=0,mn;
    bool cw[240];
    bool r1=0,r2=0,g1,g2;
    ldb mod[simp],sum,nois,biterr=0;
    int arr[120],erc;
    state* map[4][2];
    int vitc[4],rgt,vitcn[4];
    vector<bool> od[4],nw[4];
    
    for(int k=0;k<4;++k){
        map[k][0] = new state(k,0);
        map[k][1] = new state(k,1);
    }

    for(int k=0;k<cycl;++k){
        //Encoding
        for(int i=0;i<simn;++i){
            arr[i] = rand()%2;
            rg = map[rg][arr[i]]->show(cw[i*2],cw[i*2+1]);
        }
        //Transmitting
        for(int b=0;b<2*simn;++b){
            n = (cw[b]>0)?1:-1;
            sum = 0;
            nois = dist(rgen);
            for(int i=0;i<simp;++i)
                mod[i] = n*(A+nois)*coswv[i];
            for(int i=0;i<simp;++i)
                sum += mod[i]*coswv[i];
            cw[b] = sum>0;
        }
        //Decoding
        for(int i=0;i<4;++i)od[i].clear();
        vitc[0]=0, vitc[1]=vitc[2]=vitc[3]=-1;
        for(int i=0;i<simn;++i){
            for(int j=0;j<4;++j)vitcn[j] = -1;
            for(int j=0;j<4;++j){
                if(vitc[j]<0)continue;
                for(int b=0;b<2;++b){
                    rgt = map[j][b]->show(g1,g2);
                    erc = (g1!=cw[i*2]) + (g2!=cw[i*2+1]);
                    if(vitcn[rgt]<0 || vitcn[rgt] > vitc[j]+erc){
                        vitcn[rgt] = vitc[j] + erc;
                        nw[rgt] = od[j];
                        nw[rgt].push_back(b);
                    }
                }
            }
            for(int j=0;j<4;++j){
                od[j] = nw[j];
                vitc[j] = vitcn[j];
            }
        }
        //Bit Errors Counting
        mn = 0;
        for(int i=0;i<4;++i){
            if(vitc[i]<vitc[mn])
                mn = i;
        }
        for(int i=0;i<simn;++i){
            biterr += (arr[i]!=od[mn][i]);
        }
    }
    cout<<"SNR: "<<A*A/2/N<<", Errors: "<<biterr<<", BER: "<<biterr/(simn*cycl);
    cout<<" with "<<simn*cycl<<" simulations"<<endl;
    return 0;
}
