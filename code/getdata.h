#include<stdio.h>

void getdata(FILE *fpin,float *tend, float *delt,int *nlys, int *ninf, int *nneut, int *nben, float *rs, float *rl, float *rd, float *ri, float *rt, int *Ninit)

{
  char junk[20];

  fscanf(fpin,"%s %f",junk,tend);
  fscanf(fpin,"%s %f",junk,delt);
  fscanf(fpin,"%s %d",junk,nlys);
  fscanf(fpin,"%s %d",junk,ninf);
  fscanf(fpin,"%s %d",junk,nneut);
  fscanf(fpin,"%s %d",junk,nben);
  fscanf(fpin,"%s %f",junk,rs);     *rs = *rs/((float)*nben);   // s per ben allele 
  fscanf(fpin,"%s %f",junk,rl);
  fscanf(fpin,"%s %f",junk,rd);
  fscanf(fpin,"%s %f",junk,ri);
  fscanf(fpin,"%s %f",junk,rt);
  fscanf(fpin,"%s %d",junk,Ninit);
  fclose(fpin);
}