#include <iostream>
#include <vector>
#include <algorithm> 
#include<fstream>
#include<mpi.h>

#include "charge.h"

<<<<<<< HEAD
void charge_a(int me, int n, int Np, int &iBeg, int &iEnd)
{
    int q,r;
    q=n/Np;
    r=n-q*Np;
    if (me<r)
    {
        iBeg=me*(q+1);
        iEnd=(me+1)*(q+1)-1;
    }
    else
    {
        iBeg=me*q+r;
        iEnd=iBeg+q-1;
    }

=======

void charge_a(int me,int n,int np,int *ibeg,int *iend)
{
  int r=n%np;
  if(me<r){
    *ibeg=(me)*(n/np+1);
    *iend=(me+1)*(n/np+1)-1;
  }
  else{
    *ibeg=r+me*(n/np);
    *iend=*ibeg+n/np-1;
  }
>>>>>>> c73c8bc (Ajout du dossier Calcul_Parallele_final)
}