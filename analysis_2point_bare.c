// Usage: ./"program_name" (no arguments)
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#define N_file 200
#define N_righe 10
#define N_campi 8

void media_pesata(int N, double *data, double *error, double *media, double *error_media){

  int i;
  double sum_weight=0;

  (*media)=0;
  for(i=0;i<N;i++){
    (*media)+=data[i]/(error[i]*error[i]);
    sum_weight+=1./(error[i]*error[i]);
  }
  (*media)/=sum_weight;
  (*error_media)=sqrt(1./sum_weight); 
}

void media_without_errors(int N, double *data, double *media, double *error_media){

  int i;
  double error=0;

  (*media)=0;
  (*error_media)=0;

  for(i=0;i<N;i++){
    (*media)+=data[i];
    (error)+=data[i]*data[i];
  }
  (*media)/=N;
  (*error_media)=sqrt((error/(double)N-(*media)*(*media))/(N-1.)); 
}



int main(int argc, char *argv[]){

  int i,j,z,N;
  double data[N_campi][N_righe][N_file], media, error;
  FILE *f;
  char s[50];

  for(i=1;i<=N_file;i++)
  {
    sprintf(s,"dati_2pointLL/loop_SG_%d.dat",i);
    f=fopen(s,"r");
    fscanf(f,"# z: 3 h: 0.358 N: %d\n# 1:L 2:con_loop 3:error 4:con_loopAA 5:con_loopAB 6: con_loopAC 7:con_linea 8:error \n",&N);
    for(j=0;j<N_righe;j++)
    {
      for(z=0;z<N_campi;z++)
      {
	fscanf(f,"%lf ",&(data[z][j][i-1]));
      }
      fscanf(f,"\n");
    }
    fclose(f);
  }

  sprintf(s,"gnuplot/analisi_dati_2point/SG_media_bare1.dat");
  f=fopen(s,"w");
  fprintf(f,"# z: 3 h: 0.358 N: %d Nfile: %d \n# 1: L 2:con_linea 3:error \n", N, N_file);
  for(j=0;j<N_righe;j++)
  {
    fprintf(f,"%d ", (int)data[0][j][0]);
    media_pesata(N_file, data[6][j], data[7][j], &media, &error);
    fprintf(f, "%g %g ",media,error);
    fprintf(f,"\n");
  }
  fclose(f);

  return 0;
}
