#ifndef __SUM_OF_FOUR_SQUARES_DECOMPOSE_H__
#define __SUM_OF_FOUR_SQUARES_DECOMPOSE_H__

#include <gmp.h>

#define FOUR 4
#define REPS 15

void four_squares_init(mpz_t fours[FOUR]);
void four_squares_clear(mpz_t fours[FOUR]);

int special_case_p(mpz_t fours[FOUR], mpz_t input);
int iunit(mpz_t iu, const mpz_t p);
int decompose_prime(mpz_t a, mpz_t b, const mpz_t n);

// this function is guaranteed to return okay
// pre-condition: n >= 0
// returns 1 on success, 0 otherwise
// sorted 4 square roots in incresing order
int decompose(mpz_t fours[FOUR], mpz_t n);

#endif //  __SUM_OF_FOUR_SQUARES_DECOMPOSE_H__
