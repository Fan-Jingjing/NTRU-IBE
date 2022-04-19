#ifndef LIBE_SCHEME_H
#define LIBE_SCHEME_H

#include "params.h"
#include "Sampling.h"
#include "FFT.h"
#include "Random.h"
#include "Algebra.h"

void Keygen(ZZ_pX& PublicKey0,ZZX* PublicKey1,ZZX& PublicKey2,ZZX* PrivateKey);
void CompletePrivateKey(mat_ZZ& B, const ZZX * const PrivateKey);
void GPV(RR_t * v, const RR_t * const c, const RR_t s, const MSK_Data * const MSKD);
void CompleteMSK(MSK_Data * MSKD, ZZX * MSK);
void CompleteMPK(MPK_Data * MPKD, ZZ_pX MPK0, ZZX * MPK1, ZZX MPK2);
void IBE_Extract(ZZX SK_id[2], vec_ZZ id, const MSK_Data * const MSKD, const MPK_Data * const MPKD);
unsigned long IBE_Verify_Key(const ZZX SK_id[2], const vec_ZZ id, const MSK_Data * const MSKD, const MPK_Data * const MPKD);
void IBE_Encrypt(long C[2][N0], const long m[N0], const long id0[N0], const MPK_Data * const MPKD);
void IBE_Decrypt(long message[N0], const long C[2][N0], const CC_t * const SKid_FFT);
void Extract_Bench(const unsigned int nb_extr, MSK_Data * MSKD, MPK_Data * MPKD);
void Encrypt_Bench(const unsigned int nb_cryp, MPK_Data * MPKD, MSK_Data * MSKD);
void Extract_Test(const unsigned int nb_extr, MSK_Data * MSKD, MPK_Data * MPKD);
void Encrypt_Test(const unsigned int nb_cryp, MPK_Data * MPKD, MSK_Data * MSKD);

#endif
