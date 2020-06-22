// ranlxd.C
// Adapted from Luscher's ranlxd code.

#include <limits.h>
#include <float.h>
#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include "ranlxd.h"

#define BASE 0x1000000
#define MASK 0xffffff

#define STEP(pi,pj) \
        d=(*pj).c1.c1-(*pi).c1.c1-carry.c1; \
        (*pi).c2.c1+=(d<0); \
        d+=BASE; \
        (*pi).c1.c1=d&MASK; \
        d=(*pj).c1.c2-(*pi).c1.c2-carry.c2; \
        (*pi).c2.c2+=(d<0); \
        d+=BASE; \
        (*pi).c1.c2=d&MASK; \
        d=(*pj).c1.c3-(*pi).c1.c3-carry.c3; \
        (*pi).c2.c3+=(d<0); \
        d+=BASE; \
        (*pi).c1.c3=d&MASK; \
        d=(*pj).c1.c4-(*pi).c1.c4-carry.c4; \
        (*pi).c2.c4+=(d<0); \
        d+=BASE; \
        (*pi).c1.c4=d&MASK; \
        d=(*pj).c2.c1-(*pi).c2.c1; \
        carry.c1=(d<0); \
        d+=BASE; \
        (*pi).c2.c1=d&MASK; \
        d=(*pj).c2.c2-(*pi).c2.c2; \
        carry.c2=(d<0); \
        d+=BASE; \
        (*pi).c2.c2=d&MASK; \
        d=(*pj).c2.c3-(*pi).c2.c3; \
        carry.c3=(d<0); \
        d+=BASE; \
        (*pi).c2.c3=d&MASK; \
        d=(*pj).c2.c4-(*pi).c2.c4; \
        carry.c4=(d<0); \
        d+=BASE; \
        (*pi).c2.c4=d&MASK


void Ranlxd::Error(int no)
{
   switch(no)
   {
      case 0:
         printf("Error in rlxd_init\n");
         printf("Arithmetic on this machine is not suitable for ranlxd\n");
         break;
      case 1:
         printf("Error in subroutine rlxd_init\n");
         printf("Bad choice of luxury level (should be 1 or 2)\n");
         break;
      case 2:
         printf("Error in subroutine rlxd_init\n");
         printf("Bad choice of seed (should be between 1 and 2^31-1)\n");
         break;
      case 3:
         printf("Error in rlxd_get\n");
         printf("Undefined state (ranlxd is not initialized)\n");
         break;
      case 4:
         printf("Error in rlxd_reset\n");
         printf("Arithmetic on this machine is not suitable for ranlxd\n");
         break;
      case 5:
         printf("Error in rlxd_reset\n");
         printf("Unexpected input data\n");
         break;
   }
   printf("Program aborted\n");
   exit(0);
}


Ranlxd::Ranlxd(int level)
{
//  init=0;
  LuxLevel(level);
  Init(1); // default seed
}

Ranlxd::~Ranlxd()
{
}

void Ranlxd::Init(long seed)
{

  int i,k,l;
  int ibit,jbit,xbit[31];
  int ix,iy;

  if ((INT_MAX<2147483647)||(FLT_RADIX!=2)||(FLT_MANT_DIG<24)||(DBL_MANT_DIG<48)) {
    Error(0);
  }

  DefineConstants();

  i=seed;

  for (k=0;k<31;k++) {
    xbit[k]=i%2;
    i/=2;
  }

  if ((seed<=0)||(i!=0)) { Error(2); }

  ibit=0;
  jbit=18;

  for (i=0;i<4;i++) {
    for (k=0;k<24;k++) {
      ix=0;
      for (l=0;l<24;l++) {
        iy=xbit[ibit];
        ix=2*ix+iy;

        xbit[ibit]=(xbit[ibit]+xbit[jbit])%2;
        ibit=(ibit+1)%31;
        jbit=(jbit+1)%31;
      }

      if ((k%4)!=i) { ix=16777215-ix; }
      x.num[4*k+i]=ix;
    }
  }

  carry.c1=0;
  carry.c2=0;
  carry.c3=0;
  carry.c4=0;

  ir=0;
  jr=7;
  is=91;
  is_old=0;
  prm=pr%12;

}

double Ranlxd::UniformDeviate()
{
  is=next[is];
  if (is==is_old) { Update(); }
  return one_bit*((double)(x.num[is+4])+one_bit*(double)(x.num[is]));
}

void Ranlxd::LuxLevel(int level)
{

  lux_level = level;
  if (level==1) {
    pr=202;
  } else if (level==2) {
    pr=397;
  } else {
    Error(1);
  }
}

void Ranlxd::Update(void)
{
  int k,kmax,d;
  dble_vec_t *pmin,*pmax,*pi,*pj;

  kmax=pr;
  pmin=&x.vec[0];
  pmax=pmin+12;
  pi=&x.vec[ir];
  pj=&x.vec[jr];
     
  for (k=0;k<kmax;k++) {
    STEP(pi,pj);
    pi+=1;  
    pj+=1;  
    if (pi==pmax) { pi=pmin; }      
    if (pj==pmax) { pj=pmin; }
  }

  ir+=prm;
  jr+=prm;
  if (ir>=12) { ir-=12; } 
  if (jr>=12) { jr-=12; }
  is=8*ir;
  is_old=is;

}

void Ranlxd::DefineConstants(void)
{
  int k;  

  one_bit=ldexp(1.0,-24);
  for (k=0;k<96;k++) {
    next[k]=(k+1)%96;
    if ((k%4)==3) { next[k]=(k+5)%96; }
  }  

}

int Ranlxd::StateSize()
{
  return(105);
}

void Ranlxd::GetState(long* state)
{

  int k;  

  state[0]=StateSize();

  for (k=0;k<96;k++)
     state[k+1]=x.num[k];

  state[97]=carry.c1;
  state[98]=carry.c2;
  state[99]=carry.c3;
  state[100]=carry.c4;

  state[101]=pr;
  state[102]=ir;
  state[103]=jr;
  state[104]=is;

}

void Ranlxd::SetState(long* state)
{

  int k;  

  if ((INT_MAX<2147483647)||(FLT_RADIX!=2)||(FLT_MANT_DIG<24)||(DBL_MANT_DIG<48)) {
    Error(4);
  }

  DefineConstants();

  if (state[0]!=StateSize()) {
    Error(5);
  }

  for (k=0;k<96;k++) {
    if ((state[k+1]<0)||(state[k+1]>=167777216)) {
      Error(5);
    }
    x.num[k]=state[k+1];
  }

  if (((state[97]!=0)&&(state[97]!=1))||((state[98]!=0)&&(state[98]!=1))||
     ((state[99]!=0)&&(state[99]!=1))||((state[100]!=0)&&(state[100]!=1))) {
    Error(5);
  }

  carry.c1=state[97];
  carry.c2=state[98];
  carry.c3=state[99];
  carry.c4=state[100];

  pr=state[101];
  ir=state[102];
  jr=state[103];
  is=state[104];
  is_old=8*ir;
  prm=pr%12;
//  init=1;

  if (((pr!=202)&&(pr!=397))||(ir<0)||(ir>11)||
     (jr<0)||(jr>11)||(jr!=((ir+7)%12))||(is<0)||(is>91)) {
    Error(5);
  }

}


#undef BASE
#undef MASK

double Ranlxd::Normal(double mu, double sigma)
{
  double x, y, r2;

  do {
    x = 2*UniformDeviate() - 1;
    y = 2*UniformDeviate() - 1;;
    r2 = x * x + y * y;
  } while (r2 > 1.0 || r2 == 0);

  return sigma * y * sqrt (-2.0 * log (r2) / r2) + mu; // Box-Muller transform
}

std::complex<double> Ranlxd::Normal(double sigma)
{
  double x, y, r2;
  double z;

  do {
    x = 2*UniformDeviate() - 1;
    y = 2*UniformDeviate() - 1;;
    r2 = x * x + y * y;
  } while (r2 > 1.0 || r2 == 0);

  z = sigma * sqrt (-2.0 * log (r2) / r2);
  return z * std::complex<double>(x,y);
}



Ranlxd rng(2); // Global instance of the random number generator
