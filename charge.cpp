#include <iostream>
#include <vector>
#include <algorithm> 
#include<fstream>
#include<mpi.h>

#include "charge.h"

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

}