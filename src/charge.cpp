#include <iostream>
#include <vector>
#include <algorithm> 
#include<fstream>
#include<mpi.h>

#include "charge.h"


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
}