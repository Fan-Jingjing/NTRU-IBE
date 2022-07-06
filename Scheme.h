#ifndef LIBE_SCHEME_H
#define LIBE_SCHEME_H

#include "params.h"
#include "Sampling.h"
#include "FFT.h"
#include "Random.h"
#include "Algebra.h"

void Keygen(ZZ_pX& PublicKey, ZZX* PrivateKey);
void CompletePrivateKey(mat_ZZ& B, const ZZX * const PrivateKey);
vec_ZZ RandomBinaryVector();
void RandomSmallRingElement(ZZ_pX& e, const unsigned long range);
void Hash1ID(ZZ_pX& h1id, const MPK_Data * const MPKD, vec_ZZ id);
void Hash2ID(ZZ_pX& h2id, const MPK_Data * const MPKD, vec_ZZ id);
void GPV(RR_t * v, const RR_t * const c, const RR_t s, const MSK_Data * const MSKD);
void GPV(RR_t * v, const RR_t * const c, const MSK_Data * const MSKD);
void rng_uint64(uint64_t &t);
void randomShortVec_ZZ(vec_ZZ &s, const uint32_t stdev);
void randomShortVec_ZZ(vec_ZZ &s, const MSK_Data * const MSKD);
void Gram_Schmidt(RR_t Bst[2*N0][2*N0], RR_t B[2*N0][2*N0]);
void SampleD(RR_t * v, const RR_t * const c, const RR_t dev, Basis * Base);
RR_t singleBitGaussian (const uint32_t stdev);
void CompleteMSK(MSK_Data * MSKD, ZZX * MSK);
void CompleteMPK(MPK_Data * MPKD, ZZ_pX MPK);
void IBE_Extract(SKID_Data * SKID, vec_ZZ id, const MSK_Data * const MSKD, const MPK_Data *const MPKD);
unsigned long IBE_Verify_Key(const SKID_Data * SKID, const vec_ZZ id, const MSK_Data * const MSKD, const MPK_Data * const MPKD);
void IBE_Encrypt(long C[5][N0], const long m[N0], const vec_ZZ id, const MPK_Data * const MPKD);
void IBE_Decrypt(long message[N0], const long C[5][N0], SKID_Data * SKID, const MPK_Data * const MPKD);
void Extract_Bench(const unsigned int nb_extr, MSK_Data * MSKD, MPK_Data *MPKD);
void Encrypt_Bench(const unsigned int nb_cryp, MPK_Data * MPKD, MSK_Data * MSKD);
void Extract_Test(const unsigned int nb_extr, MSK_Data * MSKD, MPK_Data *MPKD);
void Encrypt_Test(const unsigned int nb_cryp, MPK_Data * MPKD, MSK_Data * MSKD);

#endif
