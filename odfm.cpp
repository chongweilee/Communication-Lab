#include<iostream>
#include<random>
#include<fstream>
#include<complex>
#include<cmath>
#include<vector>
#define pi 3.14159265358979
#define cpx complex <long double>
#define ldb long double
#define add push_back
#define vec vector
#define simn 125
#define Nc 128
#define N 1.0
#define L 30
#define pskb 1
#define psk8 3
#define qam16 4
#define md qam16
#define modnc(x) (((x%Nc)+Nc)%Nc)
#define modncl(x) (((x%(Nc+L-1))+(Nc+L-1))%(Nc+L-1))
#define flr(x) (int(x)<x?int(x):int(x)-1)
using namespace std;
vec<cpx> hl,hlt,ekn;
ldb A = 100;
int pw[10] = {1,2,4,8,16,32,64,128,256,512};
int mpe[8]={0,1,3,2,7,6,4,5},ive[8]={0,1,3,2,6,7,5,4};
int mpq[16]={0,4,12,8,1,5,13,9,3,7,15,11,2,6,14,10};
cpx cire[8],cirq[16];
default_random_engine rgen;
normal_distribution<double> dist(0,N);

void mdcal(int m){
    for(int i=0;i<8;++i)
        cire[i] = cpx(cos(pi*i/4),sin(pi*i/4));
    for(int i=0;i<16;++i)
        cirq[i] = cpx(-3,-3)+cpx(2*(i%4),2*int(i/4));
}
void gethl(){
    ldb x,y;
    ifstream ifl ("lti");
    for(int i=0;i<30;++i){
        ifl>>x>>y;
        hl.add(cpx(x,y));
    }
    for(int n=0;n<Nc;++n){
        cpx ps(0,0);
        for(int k=0;k<30;++k){
            ps += hl[k] * ekn[modnc(-1*k*n)];
        }
        hlt.add(ps);
    }
}
void getekn(){
    for(int i=0;i<Nc;++i)
        ekn.add(cpx(cos(2*pi*i/Nc),sin(2*pi*i/Nc)));
}

void dft(vec<cpx>& u,vec<cpx> ut,bool inverse){
    int j = ((inverse)?1:-1);
    for(int n=0;n<Nc;++n){
        cpx ps(0,0);
        for(int k=0;k<Nc;++k){
            ps += ut[k] * ekn[modnc(j*k*n)];
        }
        ps = ps / cpx(pow(Nc,0.5),0);
        u.add(ps);
    }
}

void cycins(vec<cpx>&x,vec<cpx> u){
    for(int i=1;i<L;++i)
        x.add(u[Nc-L+i]);
    for(int i=0;i<Nc;++i)
        x.add(u[i]);
}

void ltifc(vec<cpx>& y,vec<cpx> x){
    cpx zn,ps;
    for(int n=0;n<Nc+L-1;++n){
        zn = cpx(dist(rgen),dist(rgen));
        ps = cpx(0,0);
        for(int l=0;l<=n && l<30;++l){
            ps += hl[l]*x[n-l] + zn;
        }
        y.add(ps);
    }
}

void cycrem(vec<cpx>&u,vec<cpx> y){
    for(int i=0;i<Nc;++i)
        u.add(y[L+i-1]);
}

void detc(vec<cpx>& u,vec<cpx> ut){
    for(int i=0;i<Nc;++i)
        u.add(ut[i]/hlt[i]);
}
cpx symget(int num){
    cpx sym;
    if(md == pskb){
        sym = cpx(2*num-1,0);
    }
    else if(md == psk8){
        sym = cire[mpe[num]];
    }
    else if(md == qam16){
        sym = cirq[mpq[num]];
    }
    return cpx(real(sym)*A,imag(sym)*A);
}

void syminv(int& num,cpx sym){
    int our,oui;
    ldb pha;
    if(md == pskb){
        num = real(sym)>0;
    }
    else if(md == psk8){
        pha = arg(sym);
        pha += (pha<-pi/8)?2*pi:0;
        num = ive[int((pha+pi/8)/(pi/4))];
    }
    else if(md == qam16){
        our = flr(real(sym)/A);
        oui = flr(imag(sym)/A);
        if(our<-3 || our>3) our = (our>0)?3:-3;
        if(oui<-3 || oui>3) oui = (oui>0)?3:-3;
        num = mpq[((oui+4)/2)*4+(our+4)/2];
    }
}

void modu(vec<cpx>& m,bool psi[][8]){
    int num;
    for(int i=0;i<Nc;++i){
        num = 0;
        for(int j=0;j<md;++j){
            num += psi[i][j]*pw[j];
        }
        m.add(symget(num));
    }
}

void demodu(vec<cpx> d,bool pso[][8]){
    int num;
    for(int i=0;i<Nc;++i){
        syminv(num,d[i]);
        for(int j=0;j<md;++j){
            pso[i][j] = num&pw[j];
        }
    }
}
int main(){
    srand(time(0));
    bool psi[Nc][8],pso[Nc][8];
    getekn();
    gethl();
    mdcal(md);
    ldb bn=1,sp=2,r=0.5,ro;
    int cnt = 0;
    while(r>1*10e-6){
        cnt = 0;
        for(int s=0;s<simn*bn;++s){
            vec<cpx> ut,ifts,xn,yn,iffts,rts,yt;
            for(int i=0;i<Nc;++i){
                for(int j=0;j<md;++j){
                    psi[i][j] = rand()%2;
                }
            }
            modu(ut,psi);
            dft(ifts,ut,true);
            cycins(xn,ifts);
            ltifc(yn,xn);
            cycrem(iffts,yn);
            dft(yt,iffts,false);
            detc(rts,yt); 
            demodu(rts,pso);
            for(int i=0;i<Nc;++i){
                for(int j=0;j<md;++j){
                    cnt += psi[i][j] != pso[i][j];
                }
            }
        }
        ro = r;
        r = cnt*1.0/(Nc*simn*md*bn);
        cout<<A<<" "<<cnt<<" "<<Nc*simn*md*bn<<" "<<r<<" L:"<<L<<endl;
        break;
    }
}
