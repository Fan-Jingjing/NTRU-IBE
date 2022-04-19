#include <stdlib.h>
#include <assert.h>
#include <math.h>
#include <complex.h>
#include <time.h>
#include <NTL/ZZ.h>
#include <NTL/ZZX.h>
#include <NTL/mat_ZZ.h>
#include <gmp.h>

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
void Keygen(ZZ_pX& PublicKey0, ZZX* PublicKey1, ZZX& PublicKey2, ZZX* PrivateKey)
{
    ZZ SqNorm;
    ZZX f,g,F,G;
    ZZ_pX h;
    ZZX H[idlength+1];
    ZZX u;

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

    

    for (unsigned int i = 0; i<idlength+1; i++)
    {
       
        H[i] = RandomPolyq(N0-1);
        //cout<<deg(H[i])<<endl;;
        PublicKey1[i] = H[i];
        //cout<<i<<endl;
        
    }
    u = RandomPolyq(N0-1);
    
    PublicKey0 = Quotient(f, g);

    PublicKey2 = u;
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





void GPV(RR_t * v, const RR_t * const c, const RR_t s, const MSK_Data * const MSKD)
{

    int i;
    unsigned j;
    RR_t ci[2*N0], zi, cip, sip, aux;

    for(j=0; j<2*N0;j++)
    {
        ci[j] = c[j];
    }

    //for(j=0; j<2*N0; j++)
    //{

   // }    

    for(i=2*N0-1; i>=0; i--)
    {
        aux = (MSKD->GS_Norms)[i];
        cip = DotProduct(ci, MSKD->Bstar[i])/(aux*aux);
        sip = s/aux;
        zi = Sample4(cip, sip*PiPrime);

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



void CompleteMPK(MPK_Data * MPKD, ZZ_pX MPK0 ,ZZX * MPK1, ZZX MPK2 )
{
    MPKD->h = MPK0;
    ZZXToFFT(MPKD->h_FFT, conv<ZZX>(MPK0));
    for (unsigned i=0; i<idlength+1; i++)
    {
        MPKD->H[i] = MPK1[i];
    }
    MPKD->u = MPK2;
}



void IBE_Extract(ZZX SK_id[2], long id[idlength], const MSK_Data * const MSKD, const MPK_Data * const MPKD)
{
    unsigned int i;
    RR_t c[2*N0], sk[2*N0], sigma;
    ZZX f,g,aux;
    ZZ_pX h;
    ZZX u,h_id,mid;
    ZZX H[idlength+1];
    CC_t mid_FFT[N0], hid_FFT[N0], S_FFT[N0];

    f = MSKD -> PrK[0];
    g = MSKD -> PrK[1];
    sigma = MSKD->sigma;
    SK_id[0].SetLength(N0);
    SK_id[1].SetLength(N0);
    SK_id[2].SetLength(N0);

    //SK_id[2].SetLength(N0);
    /*
    for (i=0; i<N0; i++)
    {
        SK_id[2] = rand()%q0;
    }*/

    
    h = MPKD->h;
    u = MPKD->u;

    for (i=0;i<idlength+1;i++)
    {
        H[i] = (MPKD->H)[i];
    }
    h_id = H[0];
    //cout<<"3"<<endl;
    for (i=1;i<idlength+1;i++)
    {
        h_id += H[i]* id[i-1];
    }

    //c->u-s_3*h_id

    ZZXToFFT(S_FFT,SK_id[2]);

    ZZXToFFT(hid_FFT,h_id);
   
    for(i=0; i<N0; i++)
    {
        mid_FFT[i] = S_FFT[i]*hid_FFT[i];
    }
    FFTToZZX(mid,mid_FFT);
    

    for(i=0;i<N0;i++)
    {
        c[i] = conv<double>(u[i]-mid[i]) ;
        c[i+N0] = 0;
    }
    


    GPV(sk, c, sigma, MSKD);

    //for(i=0; i<N0; i++)
    //{
    //   sk[i] = c[i] - sk[i];
    //    sk[i+N0] = - sk[i+N0];
    //}

    for(i=0; i<N0; i++)
    {
        SK_id[0][i] = sk[i];
        SK_id[1][i] = sk[i+N0];
    }
    
}


unsigned long IBE_Verify_Key(const ZZX SK_id[3], const long id[idlength], const MSK_Data * const MSKD, const MPK_Data * const MPKD)
{
    unsigned int i;
    ZZX f,g,t,aux,u,h_id, mid;
    ZZ_pX h;
    ZZX H[idlength+1];
    CC_t mid_FFT[N0], hid_FFT[N0], S_FFT[N0];
    f = MSKD -> PrK[0];
    g = MSKD -> PrK[1];
    
    //t = conv<ZZX>(id);
    //aux = ((SK_id[0] - t)*f + g*SK_id[1])%phi;
    h = MPKD->h;
    u = MPKD->u;
    for (i=0;i<idlength+1;i++)
    {
        H[i] = (MPKD->H)[i];
    }
    h_id = H[0];
    for (i=1;i<idlength+1;i++)
    {
        h_id += H[i]* id[i-1];
    }

    ZZXToFFT(S_FFT,SK_id[2]);
    ZZXToFFT(hid_FFT,h_id);
    for(i=0; i<N0; i++)
    {
        mid_FFT[i] = S_FFT[i]*hid_FFT[i];
    }
    FFTToZZX(mid,mid_FFT);

    aux = ((SK_id[0] - u)*f + g*SK_id[1]+mid*f)%phi;
    for(i=0; i<N0; i++)
    {
        aux[i] %= q1;
    }

    if( IsZero(aux) != 0)
    {
        cout << "The signature (s1,s2,s3) doesn't verify the required equality [ (s1 - u)*f + g*s2 + s3*h_id*f = 0 ] !\nActually, (s1 - u)*f + g*s2 + s3*h_id*f = " << aux << endl << endl;
    }
    return IsZero(aux);
}


void IBE_Encrypt(long C[3][N0], const long m[N0], const long id0[idlength], const MPK_Data * const MPKD)
{

    unsigned long i;
    long r[N0], e1[N0], e2[N0], e3[N0];
    CC_t r_FFT[N0], t_FFT[N0], aux1_FFT[N0], aux2_FFT[N0], aux3_FFT[N0], hid_FFT[N0], u_FFT[N0];
    ZZX H[idlength+1];
    ZZX h_id,u;

    for(i=0; i<N0; i++)
    {
        e1[i] = (rand()%3) - 1;
        e2[i] = (rand()%3) - 1;
        e3[i] = (rand()%3) - 1;
        r[i] = (rand()%3) - 1;
    }

    MyIntFFT(r_FFT, r);
    //MyIntFFT(t_FFT, id0);
    u = MPKD->u;
    for (i=0;i<idlength+1;i++)
    {
        H[i] = (MPKD->H)[i];
    }
    h_id = H[0];

    for (i=1;i<idlength+1;i++)
    {
        h_id += H[i]* id0[i-1];
    }

    ZZXToFFT(hid_FFT,h_id);
    ZZXToFFT(u_FFT,u);

    for(i=0; i<N0; i++)
    {
        aux1_FFT[i] = r_FFT[i]*((MPKD->h_FFT)[i]);
        aux2_FFT[i] = r_FFT[i]*hid_FFT[i];
        aux3_FFT[i] = r_FFT[i]*u_FFT[i];
    }

    MyIntReverseFFT(C[0], aux1_FFT);
    MyIntReverseFFT(C[1], aux2_FFT);
    MyIntReverseFFT(C[2], aux3_FFT);


    for(i=0; i<N0; i++)
    {
        C[0][i] = (C[0][i] + e1[i]               + q0/2)%q0 - (q0/2);
        C[1][i] = (C[1][i] + e2[i]               + q0/2)%q0 - (q0/2);
        C[2][i] = (C[2][i] + e3[i] + (q0/2)*m[i] + q0/2)%q0 - (q0/2);
    } 

}


void IBE_Decrypt(long message[N0], const long C[3][N0], const CC_t * const SKid1_FFT,const CC_t * const SKid2_FFT)
{
    unsigned int i;
    CC_t c0_FFT[N0], c1_FFT[N0], aux_FFT[N0];

    MyIntFFT(c0_FFT, C[0]);
    MyIntFFT(c1_FFT, C[1]);

    for(i=0; i<N0; i++)
    {
        aux_FFT[i] = c0_FFT[i]*SKid1_FFT[i]+c1_FFT[i]*SKid2_FFT[i];
    }

    MyIntReverseFFT(message, aux_FFT);

    for(i=0; i<N0; i++)
    {
        message[i] = C[2][i] - message[i];
        message[i] = ((unsigned long)(message[i] ))%q0;
        message[i] = (message[i] + (q0>>2) )/(q0>>1);
        message[i] %= 2;
    }

}



//==============================================================================
//==============================================================================
//                             BENCHES AND TESTS
//                   FOR EXTRACTION AND ENCRYPTION/DECRYPTION
//==============================================================================
//==============================================================================


void Extract_Bench(const unsigned int nb_extr, MSK_Data * MSKD, MPK_Data * MPKD)
{
    clock_t t1, t2;
    float diff;
    unsigned int i,j;
    //vec_ZZ id;
    ZZX SK_id[3];
    long int identity[idlength];
    SK_id[2] = RandomPolyq(N0-1);
    t1 = clock();

    cout << "0%" << flush;
    for(i=0; i<nb_extr; i++)
    {
        /*
        for (j = 0; j<idlength; j++)
        {
            id[j] = rand()%2;

        }*/
        for(j=0; j<idlength; j++)
        {
            identity[j] = rand()%2;
        }
        IBE_Extract(SK_id, identity, MSKD, MPKD);
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
    //long id[idlength];
    ZZX SK_id[3], w;
    CC_t SKid_FFT[N0],SKid1_FFT[N0], SKid2_FFT[N0];
    long int message[N0], decrypted[N0];
    long int identity[N0], Ciphertext[3][N0];

/*
    for (i = 0; i<idlength; i++)
    {
        id[i] = rand()%2;
    }*/
    
    for(i=0; i<idlength; i++)
    {
        identity[i] = rand()%2;
    }


    IBE_Extract(SK_id, identity, MSKD, MPKD);
    IBE_Verify_Key(SK_id, identity, MSKD, MPKD);

    ZZXToFFT(SKid1_FFT, SK_id[1]);
    ZZXToFFT(SKid2_FFT, SK_id[2]);

    
    t1 = clock();

    cout << "0%" << flush;
    for(i=0; i<nb_cryp; i++)
    {

        for(j=0; j<N0; j++)
        {
            message[j] = (rand()%2);
        }

        IBE_Encrypt(Ciphertext, message, identity, MPKD);
        IBE_Decrypt(decrypted, Ciphertext, SKid1_FFT, SKid2_FFT);

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


void Extract_Test(const unsigned int nb_extr, MSK_Data * MSKD, MPK_Data * MPKD)
{
    unsigned int i,j , rep;
    //vec_ZZ id;
    ZZX SK_id[3];
    long int identity[idlength];

    rep = 0;

    cout << "0%" << flush;
    for(i=0; i<nb_extr; i++)
    {
        //id = RandomVector();
        /*
        for (j = 0; j<idlength; j++)
        {
            id[j] = rand()%2;
        }*/
        for(j=0; j<idlength; j++)
        {
            identity[j] = rand()%2;
        }
        IBE_Extract(SK_id, identity, MSKD, MPKD);
        rep += IBE_Verify_Key(SK_id, identity, MSKD, MPKD);
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
    //vec_ZZ id;
    ZZX SK_id[3], m, w;
    CC_t SKid_FFT[N0],SKid1_FFT[N0], SKid2_FFT[N0];
    long int id0[N0], Ciphertext[3][N0];
    //long int identity[N0];
    long int message[N0], decrypted[N0];


    //id = RandomVector();
    /*
    for (i = 0; i<idlength; i++)
    {
        id[i] = rand()%2;
    }
*/
    for(i=0; i<N0; i++)
    {
        id0[i] = rand()%2;
    }

    IBE_Extract(SK_id, id0, MSKD, MPKD);
    IBE_Verify_Key(SK_id, id0, MSKD, MPKD);
    ZZXToFFT(SKid1_FFT, SK_id[1]);
    ZZXToFFT(SKid2_FFT, SK_id[2]);


    rep = 0;

    

    cout << "0%" << flush;
    for(i=0; i<nb_cryp; i++)
    {

        for(j=0; j<N0; j++)
        {
            message[j] = (rand()%2);
        }

        IBE_Encrypt(Ciphertext, message, id0, MPKD);
        IBE_Decrypt(decrypted, Ciphertext, SKid1_FFT, SKid2_FFT);
        
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
