// This is a rewritten version of Luscher's test program using the c++ class I adapted from his code

/*******************************************************************************
*
* File testlx.c
*
* Copyright (C) 2005 Martin Luescher
*
* This software is distributed under the terms of the GNU General Public
* License (GPL)
*
* This program checks that ranlxs and ranlxd work correctly
*
*******************************************************************************/



#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include "ranlxd.h"

#define NXD 99

int main(void)
{
   int k,test;
   long *state;
   double base;
   double xd[NXD],yd[NXD],xdn[48];

   base=ldexp(1.0,48);

   Ranlxd rng(1);
   state = new long [ rng.StateSize() ];
   
   rng.Init(32767);


/*******************************************************************************
*
* Check that the correct sequences of random numbers are obtained
*
*******************************************************************************/

  for (k=0;k<20;k++) {
    for (int j=0;j<NXD;++j) { xd[j] = rng.UniformDeviate(); }
   }

   xdn[0]=135665102723086.0;
   xdn[1]=259840970195871.0;
   xdn[2]=110726726657103.0;
   xdn[3]=53972500363809.0;
   xdn[4]=199301297412157.0;
   xdn[5]=63744794353870.0;
   xdn[6]=178745978725904.0;
   xdn[7]=243549380863176.0;
   xdn[8]=244796821836177.0;
   xdn[9]=223788809121855.0;
   xdn[10]=113720856430443.0;
   xdn[11]=124607822268499.0;
   xdn[12]=25705458431399.0;
   xdn[13]=155476863764950.0;
   xdn[14]=195602097736933.0;
   xdn[15]=183038707238950.0;
   xdn[16]=62268883953527.0;
   xdn[17]=157047615112119.0;
   xdn[18]=58134973897037.0;
   xdn[19]=26908869337679.0;
   xdn[20]=259927185454290.0;
   xdn[21]=130534606773507.0;
   xdn[22]=205295065526788.0;
   xdn[23]=40201323262686.0;
   xdn[24]=193822255723177.0;
   xdn[25]=239720285097881.0;
   xdn[26]=54433631586673.0;
   xdn[27]=31313178820772.0;
   xdn[28]=152904879618865.0;
   xdn[29]=256187025780734.0;
   xdn[30]=110292144635528.0;
   xdn[31]=26555117184469.0;
   xdn[32]=228913371644996.0;
   xdn[33]=126837665590799.0;
   xdn[34]=141069100232139.0;
   xdn[35]=96171028602910.0;
   xdn[36]=259271018918511.0;
   xdn[37]=65257892816619.0;
   xdn[38]=14254344610711.0;
   xdn[39]=137794868158301.0;
   xdn[40]=269703238916504.0;
   xdn[41]=35782602710520.0;
   xdn[42]=51447305327263.0;
   xdn[43]=247852246697199.0;
   xdn[44]=65072958134912.0;
   xdn[45]=273325640150591.0;
   xdn[46]=2768714666444.0;
   xdn[47]=173907458721736.0;
   
   test=0;

   for (k=0;k<48;k++)
   {
      if (xdn[k]!=(xd[k+39]*base))
         test=1;
   }

   if (test==1)
   {
      printf("\n");
      printf("Test failed: ranlxd gives incorrect results\n");
      printf("=> do not use ranlxd on this machine\n");
      printf("\n");
   }


/*******************************************************************************
*
* Check of the I/O routines
*
*******************************************************************************/

  rng.GetState( state );

  for (k=0;k<10;k++)
  {
    for (int j=0;j<NXD;++j) { xd[j] = rng.UniformDeviate(); }
  }

  rng.SetState( state );

  for (k=0;k<10;k++)
  {
    for (int j=0;j<NXD;++j) { yd[j] = rng.UniformDeviate(); }
  }

  for (k=0;k<NXD;k++)
  {
     if (xd[k]!=yd[k])
        test=2;
  }

  if (test==2)
  {
     printf("\n");
     printf("Test failed: I/O routines for ranlxd do not work properly\n");
     printf("=> do not use ranlxd on this machine\n");
     printf("\n");
  }


/*******************************************************************************
*
* Success messages
*
*******************************************************************************/

  if (test==0)
  {
     printf("\n");
     printf("All tests on ranlxd passed\n");
     printf("=> ranlxd works correctly on this machine\n");
     printf("\n");
  }

  delete[] state;

  exit(EXIT_SUCCESS);

}

