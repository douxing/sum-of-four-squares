#include <stdio.h>
#include "decompose.h"

int from_args(int argc, char *argv[])
{
  if (argc == 1) {
    printf("usage: %s NUMBER [MORE NUMBERS ...]\n", argv[0]);
    return 0;
  }

  mpz_t t, n;
  mpz_inits(t, n, NULL);
  mpz_t fours[FOUR];
  four_squares_init(fours);

  for (int i = 1; i < argc; ++i) {
    if (gmp_sscanf(argv[i], "%Zd", n)) {
      mpz_set(t, n);
      if (decompose(fours, t)) {
  	gmp_printf("\n%d:\n%Zd =\n%Zd\n%Zd\n%Zd\n%Zd\n",
  		   i, n, fours[0], fours[1], fours[2], fours[3]);
      } else {
	printf("impossible: %Zd");
      }
    } else {
      printf("%4d: failed to parse: %s\n", i, argv[i]);
    }
  }

  mpz_clears(t, n, NULL);
  four_squares_clear(fours);
  return 0;
}

int u32()
{
  unsigned long limit = -1;
  mpz_t t, square, n;
  mpz_inits(t, square, n, NULL);
  mpz_t fours[FOUR];
  four_squares_init(fours);  

  for (unsigned long i = 0; i < 4294967296; ++i) {
    mpz_set_ui(n, i);
    mpz_set_ui(t, i);

    if (decompose(fours, t)) {
      mpz_mul(t, fours[0], fours[0]);
      mpz_mul(square, fours[1], fours[1]);
      mpz_add(t, t, square);
      mpz_mul(square, fours[2], fours[2]);
      mpz_add(t, t, square);
      mpz_mul(square, fours[3], fours[3]);
      mpz_add(t, t, square);

      if (!mpz_cmp(n, t)) {
	gmp_printf("%11d: [%Zd, %Zd, %Zd, %Zd]\n\n",
		   i, fours[0], fours[1], fours[2], fours[3]);
      } else {
	gmp_printf("error for: %d\n", i);
	break;
      }
    } else {
      gmp_printf("impossible for %d\n", i);
      break;
    }
  }

  mpz_clears(t, square, n, NULL);
  four_squares_clear(fours);
  return 0;
}

int main(int argc, char *argv[])
{
  return u32();
}
