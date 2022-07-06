#include <stdlib.h>
#include <assert.h>
#include <math.h>
#include <complex.h>
#include <time.h>
#include <NTL/ZZ.h>
#include <NTL/ZZX.h>
#include <NTL/mat_ZZ.h>
#include <gmp.h>
#include <fstream>
#include <string>

#include "Sampling.h"
#include "params.h"
#include "FFT.h"
#include "Random.h"
#include "Algebra.h"

using namespace std;
using namespace NTL;

const ZZX phi = Cyclo();


//==============================================================================
//Generates from parameters N and q :
// - a public key : polynomial h
// - a private key : polynomials f,g,F,G
//==============================================================================
void Keygen(ZZ_pX& PublicKey, ZZX* PrivateKey)
{
    ZZ SqNorm;
    ZZX f,g,F,G;

    SqNorm = conv<ZZ>(1.36*q0/2);

    GenerateBasis(f, g, F, G, SqNorm);
    PrivateKey[0] = f;
    PrivateKey[1] = g;
    PrivateKey[2] = F;
    PrivateKey[3] = G;

    for(unsigned int i=0; i<4; i++)
    {
            PrivateKey[i].SetLength(N0);
    }

    PublicKey = Quotient(f, g);
    //cout<<"PublicKey"<<PublicKey<<endl;
}

//==============================================================================
//Computes the private basis B from private key PrivateKey and parameter N
//==============================================================================
void CompletePrivateKey(mat_ZZ& B, const ZZX * const PrivateKey)
{
    ZZX f,g,F,G;
    f = PrivateKey[0];
    g = PrivateKey[1];
    F = PrivateKey[2];
    G = PrivateKey[3];

    f = -f;
    F = -F;

    B = BasisFromPolynomials(g, f, G, F);
}

vec_ZZ RandomBinaryVector()
{
    vec_ZZ w;
    unsigned int i;
    w.SetLength(N0);
    for(i=0; i<N0; i++)
    {
        w[i] = conv<ZZ>(rand())%2;
    }
    return w;
}

void RandomSmallRingVector(ZZ_pX& w, long int div)
{
    unsigned int i;
    //cout<<"w: ";
    for(i=0; i<N0; i++)
    {
        w[i] = random() % div;
        //cout<<w[i]<<" ";
    }
    //cout<<endl;
    
}

void Hash1ID(ZZ_pX& h1id, const MPK_Data * const MPKD, vec_ZZ id)  
{
    h1id = MPKD->H1[0];
    for (unsigned int i=0;i<length;i++)
    {
        if (id[i] == 1)
            h1id = h1id + MPKD->H1[i];
        else
            h1id = h1id - MPKD->H1[i]; 
    }
}

void Hash2ID(ZZ_pX& h2id, const MPK_Data * const MPKD, vec_ZZ id)  
{
    h2id = MPKD->H2[0];
    for (unsigned int i=0;i<length;i++)
    {
        if (id[i] == 1)
            h2id = h2id + MPKD->H2[i];
        else
            h2id = h2id - MPKD->H2[i]; 
    }
}


void rng_uint64(uint64_t &t)
{
    ZZ temp = RandomLen_ZZ(64);
   //t = conv<uint64_t>(temp);
    t = conv<long int>(temp);
}


RR_t singleBitGaussian (const uint32_t stdev)
{
    static double const Pi=3.141592653589793238462643383279502884L;
    static long const bignum = 0xfffffff;
    double r1, r2, theta, rr;
    uint64_t t;
    double s;
    rng_uint64(t);
        r1 = (1+(t&bignum))/((double)bignum+1);
        r2 = (1+((t>>32)&bignum))/((double)bignum+1);
        theta = 2*Pi*r1;
        rr = sqrt(-2.0*log(r2))*stdev;
        s = rr*sin(theta) + 0.5;
    return s;
}


void randomShortVec_ZZ(vec_ZZ &s, const uint32_t stdev)
{
    //const uint32_t   stdev = MSKD->sigma; 
    const uint16_t  dim = N0;
    uint16_t d2 = N0/2;
    uint16_t i;
    uint64_t t;

    static double const Pi=3.141592653589793238462643383279502884L;
    static long const bignum = 0xfffffff;
    double r1, r2, theta, rr;

    for (i=0;i<d2;i++)
    {
        rng_uint64(t);
        r1 = (1+(t&bignum))/((double)bignum+1);
        r2 = (1+((t>>32)&bignum))/((double)bignum+1);
        theta = 2*Pi*r1;
        rr = sqrt(-2.0*log(r2))*stdev;
        s[2*i] = (int64_t) floor(rr*sin(theta) + 0.5);
        s[2*i+1] = (int64_t) floor(rr*cos(theta) + 0.5);
    }

    if (dim%2 == 1)
    {
        rng_uint64(t);
        r1 = (1+(t&bignum))/((double)bignum+1);
        r2 = (1+((t>>32)&bignum))/((double)bignum+1);
        theta = 2*Pi*r1;
        rr = sqrt(-2.0*log(r2))*stdev;
        s[dim-1] = (int64_t) floor(rr*sin(theta) + 0.5);
    }  
}

void Gram_Schmidt(RR_t Bst[2*N0][2*N0], RR_t B[2*N0][2*N0])
{   
    
    RR_t v[2*N0], v1[2*N0], C_k, D_k, C_ko, D_ko, aux;
    //RR_t C[2*N0], D[2*N0];
    unsigned int j, k;

    //Reducing first vector (obvious)
    for(j=0; j<2*N0; j++)
    {    Bst[0][j] = B[0][j];    }

    //Initialising the vector v = b_N - Proj(b_N, (b_1...b_k-2) )
    for(j=0; j<N0-1; j++)
    {    v[j] = Bst[0][j+1];
         v[j+N0] = Bst[0][j+1+N0];    }
    v[N0-1] = -Bst[0][0];
    v[2*N0-1] = -Bst[0][N0];

    for(j=0; j<2*N0; j++)
    {    v1[j] = v[j];    }


    //Initialising recurring variables
    C_k = DotProduct(Bst[0], v);
    D_k = DotProduct(v, v);

    //C[0] = C_k;
    //D[0] = D_k;
    //CD[0] = C[0]/D[0];


    //Reducing b_2 to b_N and updating v at the same time
    for(k=1; k<N0; k++)
    {
        //b~k <-- r(b~_{k-1}) - <b~_{k-1},b_N>/<v_{k-1},b_N> r(v)
        aux = C_k/D_k;
        Bst[k][0] = -Bst[k-1][N0-1] + aux*v[N0-1];
        Bst[k][N0] = -Bst[k-1][2*N0-1] + aux*v[2*N0-1];
        for(j=1; j<N0; j++)
        {
            Bst[k][j] = Bst[k-1][j-1] - aux*v[j-1];
            Bst[k][j+N0] = Bst[k-1][j+N0-1] - aux*v[j+N0-1];
        }

        //v <-- v - Proj(v, b~_{k-1} )
        for(j=0; j<2*N0; j++)
        {
            v[j] -= aux*Bst[k-1][j];
        }
        //sqnorm_v -= aux*aux*SquareNorm[k-1];

        C_ko = C_k;
        D_ko = D_k;

        C_k = DotProduct(Bst[k], v1);
        D_k = D_ko - C_ko*C_ko/D_ko;

        //C[k] = C_k;
        //D[k] = D_k;
        //CD[k] = C[k]/D[k];
        //printf ("C[%d]= %Lf		", k, C_k);
        //printf ("D[%d]= %Lf\n", k, D_k);
    }



    //Reducing second half!
    //cout << "aux = " << (1<<10)/D[N0-1] << endl;
    for(j=0; j<N0; j++)
    {    Bst[N0][N0+j] = Bst[N0-1][N0-1-j]*q0/D_k;
         Bst[N0][j] = -Bst[N0-1][2*N0-1-j]*q0/D_k;    }

    //Initialising the vector v = b_N - Proj(b_N, (b_1...b_k-2) )
    for(j=0; j<N0-1; j++)
    {    v[j] = Bst[N0][j+1];
         v[j+N0] = Bst[N0][j+1+N0];    }
    v[N0-1] = -Bst[N0][0];
    v[2*N0-1] = -Bst[N0][N0];

    for(j=0; j<2*N0; j++)
    {    v1[j] = v[j];    }


    //Initialising recursive variables
    C_k = DotProduct(Bst[N0], v1);
    D_k = DotProduct(Bst[N0], Bst[N0]);

    //C[N0] = C_k;
    //D[N0] = D_k;
    //CD[N0] = C[N0]/D[N0];


    //Reducing b_2 to b_N and updating v at the same time
    for(k=N0+1; k<2*N0; k++)
    {
        //b~k <-- r(b~_{k-1}) - <b~_{k-1},b_N>/<v_{k-1},b_N> r(v)
        aux = C_k/D_k;
        Bst[k][0] = -Bst[k-1][N0-1] + aux*v[N0-1];
        Bst[k][N0] = -Bst[k-1][2*N0-1] + aux*v[2*N0-1];
        for(j=1; j<N0; j++)
        {
            Bst[k][j] = Bst[k-1][j-1] - aux*v[j-1];
            Bst[k][j+N0] = Bst[k-1][j+N0-1] - aux*v[j+N0-1];
        }
        //SquareNorm[k] = SquareNorm[k-1] - aux*aux*sqnorm_v;


        //v <-- v - Proj(v, b~_{k-1} )
        for(j=0; j<2*N0; j++)
        {
            v[j] -= aux*Bst[k-1][j];
        }
        //sqnorm_v -= aux*aux*SquareNorm[k-1];

        C_ko = C_k;
        D_ko = D_k;

        C_k = DotProduct(Bst[k], v1);
        D_k = D_ko - C_ko*C_ko/D_ko;

        //C[k] = C_k;
        //D[k] = D_k;
        //CD[k] = C[k]/D[k];
    }
}


void SampleD(RR_t * v, const RR_t * const c, const RR_t dev, Basis * Base)
{

    int i;
    unsigned j;
    RR_t ci[2*N0], zi, cip, sip, aux;

    for(j=0; j<2*N0;j++)
    {
        ci[j] = c[j];
    } 
    Gram_Schmidt(Base->Bstar, Base->B);

    for(i=2*N0-1; i>=0; i--)
    {
        aux = DotProduct((Base->Bstar)[i],(Base->Bstar)[i]);
        aux = sqrt(aux);
        cip = DotProduct(ci, Base->Bstar[i])/(aux*aux);
        sip = dev/aux;
        //zi = Sample4(cip, sip*PiPrime);// To be replaced
        zi = ( (signed int) floor(cip) ) + (int64_t) floor(singleBitGaussian(sip));

        for(j=0; j<2*N0; j++)
        {
            ci[j] -= zi*(Base->B)[i][j];
        }
    }

    for(j=0; j<2*N0; j++)
    {
        v[j] = c[j] - ci[j];
    }

}


void GPV(RR_t * v, const RR_t * const c, const RR_t s, const MSK_Data * const MSKD)
{

    int i;
    unsigned j;
    RR_t ci[2*N0], zi, cip, sip, aux;

    for(j=0; j<2*N0;j++)
    {
        ci[j] = c[j];
    } 

    for(i=2*N0-1; i>=0; i--)
    {
        aux = (MSKD->GS_Norms)[i];
        cip = DotProduct(ci, MSKD->Bstar[i])/(aux*aux);
        sip = s/aux;
        zi = Sample4(cip, sip*PiPrime);// To be replaced
        //zi = cip + singleBitGaussian(sip);

        for(j=0; j<2*N0; j++)
        {
            ci[j] -= zi*(MSKD->B)[i][j];
        }
    }

    for(j=0; j<2*N0; j++)
    {
        v[j] = c[j] - ci[j];
    }

}







//==============================================================================
//==============================================================================
//                            MAIN PROGRAMS
//==============================================================================
//==============================================================================


void CompleteMSK(MSK_Data * MSKD, ZZX * MSK)
{
    unsigned int i, j;
    mat_ZZ B0;

    for(i=0; i<4; i++)
    {
        MSKD->PrK[i] = MSK[i];
        ZZXToFFT(MSKD->PrK_fft[i], MSK[i]);
    }

    CompletePrivateKey(B0, MSK);

    for(i=0; i<2*N0; i++)
    {
        for(j=0; j<2*N0; j++)
        {
            MSKD->B[i][j] = ( (RR_t) conv<double>(B0[i][j]) );
        }
    }

    for(i=0; i<1; i++)
    {
        FastMGS(MSKD->Bstar, MSKD->B);
    }

    for(i=0; i<2*N0; i++)
    {
        MSKD->GS_Norms[i] = sqrt( DotProduct(MSKD->Bstar[i], MSKD->Bstar[i]) );
    }

    MSKD->sigma = 2*MSKD->GS_Norms[0];


}



void CompleteMPK(MPK_Data * MPKD, ZZ_pX MPK)
{
    MPKD->h = MPK;    
    for (int i=0; i<length+1; i++)
    {
        ZZ_pX temp1, temp2;
        random(temp1, N0);
        random(temp2, N0);        
        MPKD->H1[i] = temp1;
        MPKD->H2[i] = temp2; 
        //cout<<"H1|H2"<<i<<": "<<MPKD->H1[i]<<"|"<<MPKD->H2[i]<<endl;  
    }
    ZZ_pX temp3;
    random(temp3, N0);
    MPKD->u = temp3;   
    //cout<<"u: "<<MPKD->u<<endl; 
}


void IBE_Extract(SKID_Data * SKID, vec_ZZ id, const MSK_Data * const MSKD, const MPK_Data *const MPKD)
{
    unsigned int i;
    RR_t c[2*N0], sk[2*N0], sigma;
    ZZX f,g,aux;
    vec_ZZ temp;
    ZZ_pX mid1,mid2,u,h1,h2,res,zeroTest;
    RR_t B[2*N0][2*N0];
    unsigned long stddev = 1<<30;
    Basis *Base = new Basis;


    f = MSKD -> PrK[0];
    g = MSKD -> PrK[1];
    sigma = MSKD->sigma;
    u = MPKD->u;
    for(i=0; i<5; i++)
    {
        SKID->SKID[i].SetLength(N0);
    }

    temp.SetLength(N0);
    mid1.SetLength(N0);
    mid2.SetLength(N0);
    zeroTest.SetLength(N0);
    for (i=0; i<N0;i++)
    {
        zeroTest[i] = 0;
    }

    Hash1ID(h1, MPKD, id);
    Hash2ID(h2, MPKD, id);
    SKID->h1id = h1;
    SKID->h2id = h2;

    randomShortVec_ZZ(temp, stddev);
    //cout<<"sk3: ";
    for (i=0; i<N0; i++)
    {
        SKID->SKID[2][i] = conv<long>(temp[i]);
        //cout<<temp[i]<<" ";
    }
    //cout<<endl;
    if( SKID->SKID[2] == zeroTest)
    {
        cout<<"Sampling error3"<<endl;
        return;
    }

    randomShortVec_ZZ(temp, stddev);
    //cout<<"sk4: ";
    for (i=0; i<N0; i++)
    {
        SKID->SKID[3][i] = conv<long>(temp[i]);
        //cout<<temp[i]<<" ";
    }
   //cout<<endl;

    if( SKID->SKID[3] == zeroTest)
    {
        cout<<"Sampling error4"<<endl;
        return;
    }
    cout<<sigma<<endl;
    MulMod(mid1,SKID->SKID[2],h1, MPKD->Phi);
    MulMod(mid2,SKID->SKID[3],h2, MPKD->Phi);
    res = u-mid1-mid2;
    for(i=0;i<N0;i++)
    {
        ZZ temp = conv<ZZ>(res[i]);
        c[i] = ((RR_t) conv<double>(temp)) ;
        c[i+N0] = 0;
    }
    
    //GPV(sk, c, sigma*2, MSKD);
    //cout<<"SampleD begins"<<endl;

    for (i = 0; i<2*N0; i++){
        for (int j=0; j<2*N0; j++){
            //B[i][j] = 1;
            Base->B[i][j] = (MSKD->B)[i][j];
            Base->Bstar[i][j] = 0;
        }
    } 

    SampleD(sk, c, sigma, Base);
    //cout<<"SampleD ends"<<endl;

    for(i=0; i<N0; i++)
    {
        sk[i] = c[i] - sk[i];
        sk[i+N0] = - sk[i+N0];
    }
    //cout<<"s1|s2: ";
    for(i=0; i<N0; i++)
    {
        SKID->SKID[0][i] = sk[i];
        //cout<<SKID->SKID[0][i]<<"|";
        SKID->SKID[1][i] = sk[i+N0];
        //cout<<SKID->SKID[1][i]<<endl;
    }
    //cout<<endl;
    if( SKID->SKID[0] == zeroTest)
    {
        cout<<"Sampling error1"<<endl;
        return;
    }
    if( SKID->SKID[1] == zeroTest)
    {
        cout<<"Sampling error2"<<endl;
        return;
    }
    
}


unsigned long IBE_Verify_Key(const SKID_Data * SKID, const vec_ZZ id, const MSK_Data * const MSKD, const MPK_Data * const MPKD)
{
    //want to check: s1 + s2* h + s3* h1id + s4*h2id = u
    ZZ_pX h1,h2,u,h;
    ZZ_pX sk[4], temp, res;
    //string temp1;
    
    Hash1ID(h1, MPKD, id);
    Hash2ID(h2, MPKD, id);
    for (int i=0;i<4;i++)
        sk[i] = SKID->SKID[i];

    /*    
    ifstream hin("h.txt");
    ifstream uin("u.txt");
    ifstream h1in("h1.txt");
    ifstream h2in("h2.txt");
    ifstream sk1("sk1.txt");
    ifstream sk2("sk2.txt");
    ifstream sk3("sk3.txt");
    ifstream sk4("sk4.txt");

    h.SetLength(N0);
    u.SetLength(N0);
    h1.SetLength(N0);
    h2.SetLength(N0);
    sk[0].SetLength(N0);
    sk[1].SetLength(N0);
    sk[2].SetLength(N0);
    sk[3].SetLength(N0);

    for (int i=0; i<N0; i++){
        hin>>temp1;
        h[i]=stoll(temp1);
        uin>>temp1;
        u[i] = stoll(temp1);
        h1in>>temp1;
        h1[i] = stoll(temp1);
        h2in>>temp1;
        h2[i] = stoll(temp1);
        sk1>>temp1;
        sk[0][i] = stoll(temp1);
        sk2>>temp1;
        sk[1][i] = stoll(temp1);
        sk3>>temp1;
        sk[2][i] = stoll(temp1);
        sk4>>temp1;
        sk[3][i] = stoll(temp1);

    }*/
    temp = sk[0];
    MulMod(res,sk[1],MPKD->h, MPKD->Phi);
    //MulMod(res,sk[1],h, MPKD->Phi);
    temp = temp + res;
    MulMod(res,sk[2],h1, MPKD->Phi);
    temp = temp + res;
    MulMod(res,sk[3],h2,MPKD->Phi);
    temp = temp + res;
    
    if (MPKD->u==temp)
        return 0;
    else
        return 1;

}


void IBE_Encrypt(long C[5][N0], const long m[N0], const vec_ZZ id, const MPK_Data * const MPKD)
{
    // compute h1id, h2id, t, e_1, e_2, e_3, e_4, e_5
    // random gen v in R_q  
    // c_1 = vt + e_1
    // c_2 = (vh)t + e_2
    // c_3 = (vh1id)t + e_3
    // c_4 = (vh2id)t + e_4
    // c_5 = (vu)t + e_5 + q/2 m

    unsigned long i;
    ZZ_pX h,h1id, h2id, v, vt,t,mq;
    ZZ_pX c[5], e[5];
    for (i = 0; i<5; i++)
    {
        e[i].SetLength(N0);
        c[i].SetLength(N0);
    }
    v.SetLength(N0);
    t.SetLength(N0);
    vt.SetLength(N0);
    mq.SetLength(N0);
    random(v, N0);
    Hash1ID(h1id, MPKD, id);
    Hash2ID(h2id, MPKD, id);
    
    /*
    for(i=0; i<N0; i++)
    {
        cout<<v[i]<<endl;
    }*/
    RandomSmallRingVector(t,3);
    for (i = 0; i<5; i++)
    {
        RandomSmallRingVector(e[i], 3);
    }
    //cout<<"mq: ";
    for(i = 0; i<N0; i++)
    {
        if (m[i] ==0)
            mq[i] = q0/2;
        else
            mq[i] = 0;
        //cout<<mq[i]<<" ";
    }
    //cout<<endl;
    MulMod(vt,v,t,MPKD->Phi);
    c[0] = vt+e[0];
    //cout<<c[0]<<endl;
    MulMod(c[1],vt,MPKD->h,MPKD->Phi);
    c[1] = c[1]+e[1];
    //cout<<c[1]<<endl;
    MulMod(c[2],vt,h1id,MPKD->Phi);
    c[2] = c[2]+e[2];
    //cout<<c[2]<<endl;
    MulMod(c[3],vt,h2id,MPKD->Phi);
    c[3] = c[3]+e[3];
    //cout<<c[3]<<endl;
    MulMod(c[4],vt,MPKD->u,MPKD->Phi);
    c[4] = c[4]+e[4]+mq;
    //cout<<c[4]<<endl;
    //cout<<"C: ";
    for(i=0; i<5; i++)
    {
        for(int j=0; j<N0; j++)
        {
            C[i][j] = conv<long int>(c[i][j]);
            //cout<<C[i][j]<<" ";
        }
        //cout<<endl;
    }

}

void IBE_Decrypt(long message[N0], const long C[5][N0], SKID_Data * SKID, const MPK_Data * const MPKD)
{
    unsigned int i;
    //unsigned long message[N0];
    ZZ_pX c[5],aux;
    for(i=0; i<5;i++)
    {
        c[i].SetLength(N0);
    }
    aux.SetLength(N0);
   

    for(i=0; i<5; i++)
    {
        for(int j=0; j<N0; j++)
        {
            c[i][j] = C[i][j];
            //cout<<C[i][j]<<" ";
        }
        //cout<<endl;
    }
    MulMod(aux,c[0],SKID->SKID[0],MPKD->Phi);
    c[4] = c[4]-aux;
    MulMod(aux,c[1],SKID->SKID[1],MPKD->Phi);
    c[4] = c[4]-aux;
    MulMod(aux,c[2],SKID->SKID[2],MPKD->Phi);
    c[4] = c[4]-aux;
    MulMod(aux,c[3],SKID->SKID[3],MPKD->Phi);
    c[4] = c[4]-aux;
    long temp;
    for(i=0; i<N0; i++)
    {
        //message[i] = C[1][i] - message[i];
        temp = conv<long int>(c[4][i]);
        temp = temp%q0;
        temp = temp-q0/2;
        if(temp>q0/4||temp<-q0/4)
            message[i] = 1;
        else
            message[i] = 0;
    }

}



//==============================================================================
//==============================================================================
//                             BENCHES AND TESTS
//                   FOR EXTRACTION AND ENCRYPTION/DECRYPTION
//==============================================================================
//==============================================================================


void Extract_Bench(const unsigned int nb_extr, MSK_Data * MSKD, MPK_Data *MPKD)
{
    clock_t t1, t2;
    float diff;
    unsigned int i;
    vec_ZZ id;
    id = RandomBinaryVector();
    SKID_Data *SKID = new SKID_Data;

    t1 = clock();

    cout << "0%" << flush;
    for(i=0; i<nb_extr; i++)
    {
        id = RandomBinaryVector();

        IBE_Extract(SKID, id, MSKD,MPKD);
        if((i+1)%(nb_extr/10)==0)
        {
            cout << "..." << (i+1)/(nb_extr/10) << "0%" << flush;
        }
    }

    t2 = clock();
    diff = ((float)t2 - (float)t1)/1000000.0F;
    cout << "\n\nIt took " << diff << " seconds to extract " << nb_extr << " keys." << endl;
    cout << "That's " << (diff/nb_extr)*1000 << " milliseconds per key." << endl << endl;
}


void Encrypt_Bench(const unsigned int nb_cryp, MPK_Data * MPKD, MSK_Data * MSKD)
{
    clock_t t1, t2;
    double diff;
    unsigned int i,j;
    vec_ZZ id;
    //ZZX SK_id[2], w;
    long int message[N0], decrypted[N0];
    long int Ciphertext[5][N0];
    SKID_Data *SKID = new SKID_Data;

    id = RandomBinaryVector();
    
    IBE_Extract(SKID, id, MSKD,MPKD);
    IBE_Verify_Key(SKID,id, MSKD, MPKD);

    t1 = clock();

    cout << "0%" << flush;
    for(i=0; i<nb_cryp; i++)
    {

        for(j=0; j<N0; j++)
        {
            message[j] = (rand()%2);
        }

        IBE_Encrypt(Ciphertext, message, id, MPKD);
        IBE_Decrypt(decrypted, Ciphertext, SKID, MPKD);

        if((i+1)%(nb_cryp/10)==0)
        {
            cout << "..." << (i+1)/(nb_cryp/10) << "0%" << flush;
        }
    }

    t2 = clock();
    diff = ((double)t2 - (double)t1)/1000000.0l;
    cout << "\n\nIt took " << diff << " seconds to do " << nb_cryp << " encryptions and decryptions." << endl;
    cout << "That's " << (diff/nb_cryp)*1000 << " milliseconds per encryption+decryption." << endl;
    cout << "That's " << (diff/nb_cryp)*1000*1024/N0 << " milliseconds per encryption+decryption per Kilobit." << endl << endl;
}


void Extract_Test(const unsigned int nb_extr, MSK_Data * MSKD, MPK_Data *MPKD)
{
    unsigned int i, rep;
    vec_ZZ id;

    rep = 0;

    cout << "0%" << flush;
    for(i=0; i<nb_extr; i++)
    {
        id = RandomBinaryVector();
        SKID_Data *SKID = new SKID_Data;

        IBE_Extract(SKID, id, MSKD,MPKD);

        rep += IBE_Verify_Key(SKID,id, MSKD, MPKD);
        if((i+1)%(nb_extr/10)==0)
        {
            cout << "..." << (i+1)/(nb_extr/10) << "0%" << flush;
        }
    }

    cout << endl;
    if(rep == 0)
    {    cout << endl << nb_extr << " extractions successfully performed!" << endl << endl;    }
    else
    {    cout << endl << rep << " out of " << nb_extr << " extractions failed miserabily!" << endl << endl;    }
}


void Encrypt_Test(const unsigned int nb_cryp, MPK_Data * MPKD, MSK_Data * MSKD)
{
    unsigned int i, j, rep;
    vec_ZZ id;
    ZZX m, w;
    long int Ciphertext[5][N0];
    long int message[N0], decrypted[N0];


    id = RandomBinaryVector();
    SKID_Data *SKID = new SKID_Data;
    IBE_Extract(SKID, id, MSKD,MPKD);
    IBE_Verify_Key(SKID,id, MSKD, MPKD);


    rep = 0;
 

    cout << "0%" << flush;
    for(i=0; i<nb_cryp; i++)
    {

        for(j=0; j<N0; j++)
        {
            message[j] = (rand()%2);
        }

        IBE_Encrypt(Ciphertext, message, id, MPKD);
        IBE_Decrypt(decrypted, Ciphertext, SKID, MPKD);

        for(j=0; j<N0; j++)
        {
            if(message[j] != decrypted[j])
            {
                cout << "ERROR : Dec(Enc(m)) != m " << endl;
                rep++;
                break;
            }
        }

        if((i+1)%(nb_cryp/10)==0)
        {
            cout << "..." << (i+1)/(nb_cryp/10) << "0%" << flush;
        }
    }

    cout << endl;
    if(rep == 0)
    {    cout << endl << nb_cryp << " encryptions+decryptions successfully performed!" << endl << endl;    }
    else
    {    cout << endl << rep << " out of " << nb_cryp << " encryptions+decryptions failed miserabily!" << endl << endl;    }
}
