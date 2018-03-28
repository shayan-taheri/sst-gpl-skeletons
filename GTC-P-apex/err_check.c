#include <libgen.h>
#include "bench_gtc.h"

static real err_check_FOM=0.0;

int err_check_setFOM( real f1, real f2 ){
  err_check_FOM = f1 * f2;
  return 0;
}

int err_check_print( const char *fname, int Nparticle ){

  real FOM, expected, diff, rel_err;
  char tname[256]="", *bname=NULL;
  const real correct[4] = { 0.0, 2.52889981e+02, 1.26978981e+05, 1.72485862e+09 };
	
  printf(" Checking answers after last time step.... \n");
  strncpy( tname, fname, 256 );
  bname = basename( tname );

  int icorrect = 0;
  if( strcmp(bname,"A.txt")==0 && Nparticle== 30 ) icorrect=1;
  if( strcmp(bname,"E.txt")==0 && Nparticle==200 ) icorrect=2;
  if( strcmp(fname,"E.txt")==0 && Nparticle==400 ) icorrect=3; 
  if( icorrect==0 ){
    printf(" Skipping error check... could not identify problem size. \n");
    return 0;
  }
  expected = correct[icorrect];

  FOM = err_check_FOM;
  diff = fabs(expected-FOM);
  rel_err = fabs( diff / FOM );

  printf(" Expected Answer: %15.8e \n", expected );
  printf(" This Run       : %15.8e \n", FOM      );
  printf(" Difference     : %15.8e \n", diff     );
  printf(" Relative Error : %15.8e %5s\n", 
	 rel_err, ( rel_err<1.0e-6 ? "PASS" : "FAIL" ) );
  printf("\n");

  return 0;
} //end err_check
