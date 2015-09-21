// v.1.2.2

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <iostream>
using namespace std;
#include "classes.h"
#define __MAIN
#include "params.h"

#if (defined(GILLESPIE) && defined(FASTER_KMC)) || (defined(GILLESPIE) && defined(NORMAL)) || (defined(NORMAL) && defined(FASTER_KMC))
  #error too many methods defined!
#endif

#if !defined(GILLESPIE) && !defined(FASTER_KMC) && !defined(NORMAL) 
  #error no method defined!
#endif

extern char *NUM ; 
extern int RAND, sample, treatment, max_size ;
extern double tt ;
extern float time_to_treat ;

int sample=0 ;

void save_positions(char *name, float dz) 
{
  FILE *data=fopen(name,"w") ;
  for (int i=0;i<cells.size();i++) {
    Lesion *ll=lesions[cells[i].lesion] ;
    Genotype *g=genotypes[cells[i].gen] ;
    if (abs(int(cells[i].z+ll->r.z))<dz || cells.size()<1e4) fprintf(data,"%d %d %d %u\n",int(cells[i].x+ll->r.x), int(cells[i].y+ll->r.y),int(cells[i].z+ll->r.z),genotypes[cells[i].gen]->index) ;
  }        
  fclose(data) ; 
}

float save_2d_image(char *name)
{
  int i,j,k;
  int minx=1<<20,maxx=-minx,miny=minx,maxy=-minx,minz=minx,maxz=-minx ;  
  for (i=0;i<cells.size();i++) {
    Lesion *ll=lesions[cells[i].lesion] ;
    if (cells[i].x+ll->r.x<minx) minx=int(cells[i].x+ll->r.x) ;
    if (cells[i].x+ll->r.x>maxx) maxx=int(cells[i].x+ll->r.x) ;
    if (cells[i].y+ll->r.y<miny) miny=int(cells[i].y+ll->r.y) ;
    if (cells[i].y+ll->r.y>maxy) maxy=int(cells[i].y+ll->r.y) ;
    if (cells[i].z+ll->r.z<minz) minz=int(cells[i].z+ll->r.z) ;
    if (cells[i].z+ll->r.z>maxz) maxz=int(cells[i].z+ll->r.z) ;
  }
  maxx++ ; maxy++ ; 
  float density=1.*cells.size()/(float(maxx-minx)*float(maxy-miny)*float(maxz-minz)) ;
  if (cells.size()<1e3) return density ;
  float diam=pow(float(maxx-minx)*float(maxy-miny)*float(maxz-minz),1./3) ;
  printf("density=%f\n",density) ;

  int nnn=(maxx-minx)*(maxy-miny) ;
  if (float(maxx-minx)*float(maxy-miny)>2e9) err("(maxx-minx)*(maxy-miny) too large",float(maxx-minx)*float(maxy-miny)) ;
  int *types=new int[nnn] ;
  short int *zbuf=new short int[nnn] ;
  BYTE *br=new BYTE[nnn] ;
  for (i=0;i<nnn;i++) { zbuf[i]=minz ; types[i]=-1 ; br[i]=255 ; }
  printf("%d %d\t %d %d\n",minx,maxx,miny,maxy) ;
  printf("%d x %d = %d\n",maxx-minx,maxy-miny,nnn) ;

  for (i=0;i<cells.size();i++) {
    Lesion *ll=lesions[cells[i].lesion] ;
    int z=int(cells[i].z+ll->r.z) ;
    int adr=int(cells[i].y+ll->r.y-miny)*(maxx-minx)+int(cells[i].x+ll->r.x-minx) ;
    if (adr<0 || adr>=nnn) err("adr",adr) ;
    if (zbuf[adr]<z) { zbuf[adr]=z ; types[adr]=genotypes[cells[i].gen]->index ; }
  }
  
  vecd li(1,-1,-0.3) ;
  normalize(li) ;
  float mult=1.2, dmul=0.93 ; 
  float range=-(0.916291/(density*log(dmul))) ; 
  if (range>diam) { range=diam ; dmul=pow(0.4,1/(density*range)) ; }
  for (i=0;i<maxy-miny;i++) 
    for (j=0;j<maxx-minx;j++) {
      int il,jl,kl ;
      float d=1,o ; 
      k=zbuf[i*(maxx-minx)+j] ;
      for (o=1;o<range;o++) {
        il=int(i-li.y*o) ; jl=int(j-li.x*o) ;  
        if (il>0 && jl>0 && il<maxy-miny && jl<maxx-minx) {
          kl=int(k-li.z*o) ;
          if (types[il*(maxx-minx)+jl]!=-1 && zbuf[il*(maxx-minx)+jl]>kl) d*=dmul ;
        } 
        if (d<0.4) { d=0.4 ; break ; }
      }
      br[i*(maxx-minx)+j]=br[i*(maxx-minx)+j]*d ;
    }

  FILE *data=fopen(name,"w") ;    
  for (i=0;i<maxy-miny;i++) {
    for (j=0;j<maxx-minx;j++) fprintf(data,"%d %d ",types[i*(maxx-minx)+j],br[i*(maxx-minx)+j]) ;
    fprintf(data,"\n") ;
  }
  fclose(data) ; 
  delete [] types ; delete [] zbuf ; delete [] br ;
  return density ;
}
  

void save_genotypes(char *name)
{
  FILE *data=fopen(name,"w") ;
  for (int i=0;i<genotypes.size();i++) {
    Genotype *g=genotypes[i] ;
    if (g!=NULL && g->number>0) {
      fprintf(data,"%d  %d  %d %d  %d\t",i, g->prev_gen,g->no_resistant,g->no_drivers, g->number) ;
      for (int j=0;j<g->sequence.size();j++) fprintf(data," %u",g->sequence[j]) ; 
      fprintf(data,"\n") ;
    } 
  }
  fclose(data) ;  
}

void save_most_abund_gens(char *name, int *most_abund)
{
  FILE *data=fopen(name,"w") ;
  for (int i=0;i<genotypes.size();i++) {
    Genotype *gg=genotypes[i] ;
    if (gg!=NULL && gg->number>0) {
      int r=0,g=0,b=0 ;
      for (int j=0;j<gg->sequence.size();j++) {
        if ((gg->sequence[j]&L_PM)==most_abund[0]) r=1 ;
        if ((gg->sequence[j]&L_PM)==most_abund[1]) g=1 ; 
        if ((gg->sequence[j]&L_PM)==most_abund[2]) b=1 ;
      }
      if (r || g || b) fprintf(data,"%d %d %d\t%d\n",r,g,b,gg->index) ;
    }
  }
  fclose(data) ;  

}

int main(int argc, char *argv[])
{
#if defined(GILLESPIE)   
  cout <<"method: GILLESPIE\n" ;
#endif
#if defined(FASTER_KMC)   
  cout <<"method: FASTER_KMC\n" ;
#endif
#if defined(NORMAL)   
  cout <<"method: NORMAL\n" ;
#endif

  int nsam ;
  if (argc!=4) { err(" Error:: arguments needed: name, no_samples, RAND. Program terminated. \n"); } 
  else { 
    NUM=argv[1] ;
    nsam=atoi(argv[2]) ;
    RAND=atoi(argv[3]) ;
  }
  cout <<NUM<<" "<<" "<<nsam<<" "<<RAND<<endl ;
  _srand48(RAND) ;
  init();
  char name[256] ;
  sprintf(name,"%s/each_run_%d.dat",NUM,max_size) ;
  FILE *er=fopen(name,"w") ; fclose(er) ;
  for (sample=0;sample<nsam;sample++) { 
    reset() ;
#ifdef MAKE_TREATMENT_N
    int s=0 ; while (main_proc(max_size,-1,-1, 10)==1) { s++ ; reset() ; } ; // initial growth until max size is reached, saved every 10 days
    if (s>0) printf("resetted %d times\n",s) ;
    save_data() ; 
    treatment=1 ;       
    double max_time=2*tt ;
    main_proc(1.25*max_size,-1,max_time, 10) ; // treatment
#elif defined MAKE_TREATMENT_T
    int s=0 ; while (main_proc(-1,-1,time_to_treat, 10)==1) { s++ ; reset() ; } ; // initial growth until max time is reached, saved every 10 days
    if (s>0) printf("resetted %d times\n",s) ;
    save_data() ; 
    treatment=1 ;       
    double max_time=2*tt ;
    max_size=cells.size()*1.25 ;
    main_proc(max_size,-1,max_time, 10) ; // treatment

#else    
    int s=0 ; while (main_proc(max_size,2,-1, -1)==1) { s++ ; reset() ; } // initial growth until max size is reached
    if (s>0) printf("resetted %d times\n",s) ;
    fflush(stdout) ;
    save_data() ; 
    
    sprintf(name,"%s/each_run_%d.dat",NUM,max_size) ;
    FILE *er=fopen(name,"a") ;
    fprintf(er,"%d\t%d %lf\n",sample,s,tt) ;
    fclose(er) ;

    // save some more data
    
    int *snp_no=new int[L], *snp_drivers=new int[L] ; // array of SNPs abundances
    for (int i=0;i<L;i++) { snp_no[i]=snp_drivers[i]=0 ; }
    for (int i=0;i<genotypes.size();i++) {
      if (genotypes[i]!=NULL && genotypes[i]->number>0) 
        for (int j=0;j<genotypes[i]->sequence.size();j++) {
          snp_no[((genotypes[i]->sequence[j])&L_PM)]+=genotypes[i]->number ;      
          if (((genotypes[i]->sequence[j])&DRIVER_PM)) snp_drivers[((genotypes[i]->sequence[j])&L_PM)]+=genotypes[i]->number ;
        }
    }

    save_spatial(snp_no) ;

    printf("saving PMs...\n") ;
    int most_abund[100] ;
    sprintf(name,"%s/all_PMs_%d_%d.dat",NUM,RAND,sample) ; save_snps(name,snp_no,max_size,0,most_abund) ;
    if (driver_adv>0 || driver_migr_adv>0) { printf("saving driver PMs...\n") ; sprintf(name,"%s/drv_PMs_%d_%d.dat",NUM,RAND,sample) ; save_snps(name,snp_drivers,max_size,0,NULL) ; }
    delete [] snp_no ; delete [] snp_drivers ;

    if (nsam==1) {  // do this only when making images of tumours
      printf("saving images...\n") ;
      int j=0 ;
      for (int i=0;i<genotypes.size();i++) {
        if (genotypes[i]!=NULL && genotypes[i]->number>0) genotypes[i]->index=j++ ; 
      }       
      sprintf(name,"%s/2d_image_%d.dat",NUM,max_size) ; float density=save_2d_image(name) ;
      sprintf(name,"%s/cells_%d.dat",NUM,max_size) ; save_positions(name,1./density) ; 
      sprintf(name,"%s/genotypes_%d.dat",NUM,max_size) ; save_genotypes(name) ;
      sprintf(name,"%s/most_abund_gens_%d.dat",NUM,max_size) ; save_most_abund_gens(name,most_abund) ;
    }
#endif
  } 
  end() ;
	return 0 ;
}
