#include<stdio.h>
#include<math.h>
#include<stdlib.h>
#include<strings.h>
#include"random.h"
#include "getdata.h"

#define MAXNS 31000
#define MAXGENES 40
#define MAXPRINT 2000

long seed=-1;

int main(int argc, char **argv)
{

  float tend;          
  float delt;
  int nsteps;
  int nlys;
  int ninf;
  int nneut;
  int nben;
  int ngenes;
  float rs;
  float rl;
  float rd;
  float ri;
  float rt;
  int Ninit;
  
  short (*prophages)[MAXGENES], (*tmpptr)[MAXGENES];
  short (*newprophages)[MAXGENES];
  short pro1[MAXNS][MAXGENES];
  short pro2[MAXNS][MAXGENES];
  float genemeans[MAXGENES];
  float ismeans[MAXGENES];
  float tempsum[MAXGENES];
  float tempsumt[MAXGENES];
  int occupied[MAXNS];
  int induce[MAXNS];
  int inducible, induciblesum, numben, noinduceflag = 0;
  float ss, css[MAXNS], r, fractiontolose;
  int maxind, i, newi, j, k, ii, jj, lyscapable;
  float ran3(long *);
  void getdata(FILE *,float *, float *,int *, int *, int *, int *, float *, float *, float *, float *, float *, int *);
  int ntoprint, sizeflag=1, sizes[MAXNS], sizehist[MAXGENES+1], ksum=0;

  FILE *fpin, *fpout, *fpout2, *fpout3, *fpout4, *fpout5;

 if (argc>1) seed = -((long)(atof(argv[1])));
 if (argc>2) sizeflag = (int)(atof(argv[2]));

  if ((fpout=fopen("prosim.out","w"))==NULL) {
    fprintf(stderr,"Error opening prosim.out\n");
    printf("\a");
    exit(1);
  }
  if ((fpout2=fopen("genemeans.out","w"))==NULL) {
    fprintf(stderr,"Error opening genemeans.out\n");
    printf("\a");
    exit(1);
  }
  if ((fpin=fopen("prosim.in","r"))==NULL) {
    fprintf(stderr,"Error opening prosim.in\n");
    printf("\a");
    exit(1);
  }
  if (sizeflag)
  if ((fpout3=fopen("sizes.out","w"))==NULL) {
    fprintf(stderr,"Error opening sizes.out\n");
    printf("\a");
    exit(1);
  }
  if ((fpout4=fopen("Ns.out","w"))==NULL) {
    fprintf(stderr,"Error opening Ns.out\n");
    printf("\a");
    exit(1);
  }  
  if ((fpout5=fopen("ismeans.out","w"))==NULL) {
    fprintf(stderr,"Error opening ismeans.out\n");
    printf("\a");
    exit(1);
  }  

  getdata(fpin,&tend,&delt,&nlys,&ninf,&nneut,&nben,&rs,&rl,&rd,&ri,&rt,&Ninit);
  if (Ninit>=MAXNS) {
    fprintf(stderr,"Error, Ninit too large\n");
    exit(1);
  }  
  nsteps = tend/delt;
  ngenes = nlys+ninf+nneut+nben;
  if (ngenes>=MAXGENES) {
    fprintf(stderr,"Error, too many genes\n");
    exit(1);
  }
  rs = rs*delt;
  rl = rl*delt;
  rd = rd*delt;
  ri = ri*delt;
  rt = rt*delt;
  float big = 0.02;
  if ((rs>big) || (rl>big) || (rd>big) || (rs*nben>big))
    fprintf(stdout,"Error: big changes in one timestep\n");
  prophages = pro1;
  newprophages = pro2;
  ntoprint = (int)((float)nsteps/(float)MAXPRINT);
  if (ntoprint == 0) ntoprint = 1;

  //-------------------------------------------------------------------------------
  for (int i = 0; i < Ninit; i++) {
    occupied[i] = 1;
    for (int j = 0; j < ngenes; j++) prophages[i][j] = 1;
    }	
  maxind = Ninit;    // max ind is the maximum possible occupied row in prophages		 
//----------TIME LOOP-----------------------------------------------------------
    for (int istep =0; istep < nsteps; istep++){  	

//------------------INDUCTION---------------------------------------------------------
//
// scan through the prophages, prophage is inducible only if all genes reqd for induction are present
// also keep track that at least something is inducible
    induciblesum = 0;
    for (i = 0; i < maxind; i++) {
      inducible = prophages[i][0]; 
      for (j = 1; j < nlys; j++)
	if (prophages[i][j] != -1) inducible = inducible*prophages[i][j];
      if ( (ran3(&seed) < ri) && (inducible >0)) {
	   occupied[i] = 0; 
           for (j=0;j<ngenes;j++) prophages[i][j]=0; }
      induciblesum += inducible;
    }

    if ((noinduceflag ==0) && (induciblesum == 0)) {
      noinduceflag =1 ;
      fprintf(stdout,"There is no more inducible phage in the population at timestep %d of %d.\n",istep,nsteps);
    } 

//------------------------------------------------------------------------------------
//knockout the genes that have been degraded
    for (i = 0; i < maxind; i++)
      for (j = 0; j < ngenes; j++)
        if (ran3(&seed) < rd)  prophages [i][j] = 0;

    //if all genes from a given prophage have been knocked out replace the corresponding entry at "occupied" by 0.
    for (i = 0; i < maxind; i++) {
      int sumpro = 0;
      for (j = 0; j < ngenes; j++) sumpro += prophages[i][j];
      if (sumpro == 0) occupied[i] = 0;
    }

//IS insertions change the sequence to -1
    for (i = 0; i < maxind; i++)
      for (j = 0; j < ngenes; j++)
        if (ran3(&seed) < rt)  prophages [i][j] = -1;

 
//  add new prophages if they have lysis and infection genes 
  j = 0;   // first possible place to put the new prophage
  for (i=0; i<maxind; i++) {
    lyscapable = 1;
    for (k=0; k<nlys+ninf; k++)
      if ((prophages[i][k]==-1)||(prophages[i][k]==0)) lyscapable = 0;
    if (lyscapable==1)  
      if (ran3(&seed) < rl) {   // make a new copy of prophage[i]
        while (occupied[j]==1) j++; //find the next empty spot
        occupied[j] = 1;
        for (k=0; k<ngenes; k++) prophages[j][k] = prophages[i][k]; 
      }
   }   
   if (j>=maxind) maxind = j+1;
   if (j>MAXNS) { fprintf(stderr,"Error: MAXNS exceeded\n"); istep = nsteps;}          

//----------- SELECTION
// put a copy of each prophage into the next genn w prob rs*(num ben genes)
//  newi will count the number of prophages in the next generation
   newi = -1;
   for (i=0;i<maxind;i++) {
     numben = 0;
     for (j=ngenes-nben;j<ngenes;j++)
       if (prophages[i][j]==1) numben++;
     if (ran3(&seed)<(float)(rs*numben)) {
       newi++;
       for (k=0;k<ngenes;k++) {
         newprophages[newi][k] = prophages[i][k];
       }
       occupied[newi]=1;
       }
   }
/* population size regulation:  Every prophage is copied to the next
   generation with high probability.  If the current population < Ninit,
   every prophage is copied.  If the current population > Ninit, the
   probability is reduced so that on average Ninit are maintained  */
   
   fractiontolose = 1.0-(float)Ninit/(maxind+newi);
   for (i=0;i<maxind;i++) {
     if (ran3(&seed)>fractiontolose) {
       newi++;
       ksum = 0;
       for (k=0;k<ngenes;k++) {
         newprophages[newi][k] = prophages[i][k];
         ksum += prophages[i][k];
       }
       if (ksum>0) occupied[newi]=1;
     }
   }
   for (i=newi+1;i<maxind+1;i++) occupied[i] = 0;   // after newi, unoccupied
   tmpptr = prophages;
   prophages = newprophages;
   newprophages = tmpptr;
   maxind = newi;   // newi is the population size of the next population


   if ((float)istep/ntoprint == (int)istep/ntoprint) {
     for (jj = 0; jj < maxind; jj++) sizes[jj] = 0;  //initialize
     for (ii = 0; ii < ngenes; ii++) {
      tempsum[ii] = 0;
      tempsumt[ii] = 0;   // for transposase genes
      sizehist[ii]=0;     // initialize for later
      for (jj = 0; jj <maxind; jj++) {
        if (prophages[jj][ii] == -1) {
	  tempsumt[ii]++;
	  sizes[jj]++;
	}
	else  {
	  tempsum[ii]=tempsum[ii]+ prophages[jj][ii];
	  sizes[jj] += prophages[jj][ii];
	}
      }	
      tempsum[ii] = (float)(tempsum[ii]/((float)maxind));
      tempsumt[ii] = (float)(tempsumt[ii]/((float)maxind));
    }
    sizehist[ngenes]=0;  // last entry in array didn't get initialized yet
    for (jj=0;jj<maxind;jj++)  sizehist[sizes[jj]]++;
    genemeans[0] = 0; genemeans[1] = 0; genemeans[2] = 0 ; genemeans[3] = 0;
    for (i = 0; i < nlys; i++) genemeans[0] = genemeans[0] + tempsum[i];
    for (i = nlys; i < nlys+ninf; i++) genemeans[1] = genemeans[1] + tempsum[i];
    for (i = nlys+ninf; i < nlys+ninf+nneut; i++) genemeans[2] = genemeans[2] + tempsum[i];
    for (i = nlys+ninf+nneut; i < ngenes; i++) genemeans[3] = genemeans[3] + tempsum[i];
    fprintf(fpout2,"%f %f %f %f %f\n",delt*istep,genemeans[0],genemeans[1],genemeans[2],genemeans[3]);
    ismeans[0] = 0; ismeans[1] = 0; ismeans[2] = 0 ; ismeans[3] = 0;
    for (i = 0; i < nlys; i++) ismeans[0] = ismeans[0] + tempsumt[i];
    for (i = nlys; i < nlys+ninf; i++) ismeans[1] = ismeans[1] + tempsumt[i];
    for (i = nlys+ninf; i < nlys+ninf+nneut; i++) ismeans[2] = ismeans[2] + tempsumt[i];
    for (i = nlys+ninf+nneut; i < ngenes; i++) ismeans[3] = ismeans[3] + tempsumt[i];
    fprintf(fpout5,"%f %f %f %f %f\n",delt*istep,ismeans[0],ismeans[1],ismeans[2],ismeans[3]);

    if (sizeflag) {
      for (i=0;i<=ngenes;i++) fprintf(fpout3,"%d ",sizehist[i]);
      fprintf(fpout3,"\n");
    }
    fprintf(fpout4,"%d\n",maxind);
   } //end toprint if statement
   if (sizehist[0] == maxind) {
       fprintf(stdout,"Prophage population is extinct\n");
       istep = nsteps;
   }    
  }  //end of TIME loop on (istep)
  for (i = 0; i < maxind; i++){
      //fprintf(fpout,"%d ",occupied[i]);
      for (j  =0; j < ngenes; j++) fprintf(fpout," %d ",prophages[i][j]);
      fprintf(fpout,"\n");
  }   
  fclose(fpout2);
  fclose(fpout);
  printf("\a");
}  // end of main
