#include <stdio.h>
#include <stdlib.h>

#include "nrdef.h"
#include "nrutil.h"

#include "vnrdef.h"
#include "vnrutil.h"

#include "mutil.h"

uint8 p1(uint8 M, uint8 I){
	
	if(M > I)
		M -- ;
	else if(M < I)
		M ++ ;	
	return M;  
}

uint8 p2(uint8 M, uint8 I){
	
	uint8 O;
	if(M >= I)
		O = M - I ;
	else if(M < I)
		O = I - M ;	
	return O;  
}

uint8 p3(uint8 V, uint8 O, int N){

	if(V < O*N)
		V ++ ;
	else if(V > O*N)
		V -- ;
	return V;  
}

uint8 p4(uint8 V, uint8 O){
	uint8 E;
	if(O < V)
		E = 0 ;
	else
		E = 255 ;
	return E;  
}


void erosion(uint8** vX, uint8** vY,int nrl,int nrh, int ncl, int nch) {
  uint8 a0, b0, c0;
  uint8 a1, b1, c1;
  uint8 a2, b2, c2;
	
  uint8 row0, row1, row2;
  uint8 y;
  
	
  for (int i = nrl; i < nrh; i++) {
      a0 = vX[i-1][0-1];
      b0 = vX[i-1][0];
      a1 = vX[i][0-1];
      b1 = vX[i][0];
      a2 = vX[i+1][0-1];
      b2 = vX[i+1][0];
      for (int j = ncl; j < nch; j++) {
      	c0 = vX[i-1][j+1];
      	c1 = vX[i][j+1];
      	c2 = vX[i+1][j+1];
      	
      row0 = a0*b0*c0;
      row1 = a1*b1*c1;
      row2 = a2*b2*c2;
      y = row1*row2*row0;
      vY[i][j] = y;

      //rotations de registres
      a0=b0;
      a1=b1;
      a2=b2;
      b0=c0;
      b1=c1;
      b2=c2;
   }
  }


 
}
void dilatation(uint8** vX, uint8** vY,int nrl,int nrh, int ncl, int nch) {
  
  uint8 a0, b0, c0;
  uint8 a1, b1, c1;
  uint8 a2, b2, c2;
  
  uint8 row0, row1, row2;
  uint8 y;
  
  for (int i = nrl; i < nrh; i++) {
      a0 = vX[i-1][0-1];
      b0 = vX[i-1][0];
      a1 = vX[i][0-1];
      b1 = vX[i][0];
      a2 = vX[i+1][0-1];
      b2 = vX[i+1][0];
      for (int j = ncl; j < nch; j++) {
        c0 = vX[i-1][j+1];
        c1 = vX[i][j+1];
        c2 = vX[i+1][j+1];
        
      row0 = a0+b0+c0;
      row1 = a1+b1+c1;
      row2 = a2+b2+c2;
      y = row1+row2+row0;
      vY[i][j] = y;
      //rotations de registres
      a0=b0;
      a1=b1;
      a2=b2;
      b0=c0;
      b1=c1;
      b2=c2;
    }
  }
}

int main(int argc, char * argv[])
{
    int N = 4;
    int i,j,k;


    // nrl, nrh , ncl , nch sont row_low row_high col_low col_high pour matrix uint8matrix
    int nrl, nrh, ncl, nch;
    int row,col;
    char filename[100];
    //on load de l'image du fond pour M 
    k = 3000;
    sprintf(filename,"car3/car_%d.pgm",k);
    uint8 **M = LoadPGM_ui8matrix(filename, &nrl, &nrh, &ncl, &nch);
    row = nrh - nrl + 1;
    col = nch - ncl + 1;


    //on alloue pour IEOV
    uint8 **I = ui8matrix(-1,row,-1,col);
    uint8 **E = ui8matrix(-1,row,-1,col);
    uint8 **O = ui8matrix(-1,row,-1,col);
    uint8 **V = ui8matrix(-1,row,-1,col);
    uint8 **ER = ui8matrix(-1,row,-1,col);
    uint8 **DI = ui8matrix(-1,row,-1,col);
   
   


    //initialiser V à 2
    for(i = nrl; i < nrh; i++){
        for(j = ncl;  j < nch ; j++){
            V[i][j] = 2;
        }
    }


    //pour chaque image, on le load à I
    for(k = 3000; k < 3200; k++)
    {
        sprintf(filename,"car3/car_%d.pgm",k);
        MLoadPGM_ui8matrix(filename,nrl, nrh, ncl, nch,(uint8**)I); 


		for( i = nrl; i < nrh; i ++){
			for( j = ncl; j < nch; j ++){
				uint8 Mt = M[i][j];
				uint8 It = I[i][j];
				uint8 Vt = V[i][j]; 
				uint8 M_t = p1(Mt,It);
				uint8 Ot = p2(Mt,It);
				Vt = p3(Vt,Ot,N);
				uint8 Et = p4(Vt,Ot);
				M[i][j] = M_t;
				O[i][j] = Ot;
				V[i][j] = Vt;
				E[i][j] = Et;
				
 				
			}
		}

    

        sprintf(filename,"result/M_%d.pgm",k);
        SavePGM_ui8matrix((uint8**)M, nrl, nrh, ncl, nch,filename);
        sprintf(filename,"result/O_%d.pgm",k);
        SavePGM_ui8matrix((uint8**)O, nrl, nrh, ncl, nch,filename);
        sprintf(filename,"result/V_%d.pgm",k);
        SavePGM_ui8matrix((uint8**)V, nrl, nrh, ncl, nch,filename);
        sprintf(filename,"result/E_%d.pgm",k);
        SavePGM_ui8matrix((uint8**)E, nrl, nrh, ncl, nch,filename);
		
	//optimisation
	erosion(E,ER,nrl,nrh,ncl,nch);
	dilatation(ER,DI,nrl,nrh,ncl,nch);
	sprintf(filename,"result/ER_%d.pgm",k);
        SavePGM_ui8matrix((uint8**)ER, nrl, nrh, ncl, nch,filename);
	sprintf(filename,"result/DI_%d.pgm",k);
        SavePGM_ui8matrix((uint8**)DI, nrl, nrh, ncl, nch,filename);
	printf("image = %d, files in the %s , N = %d\n", k, filename,N );
    }




    return 0;
}







