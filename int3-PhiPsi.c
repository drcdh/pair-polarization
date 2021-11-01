#include <stdio.h>
#include <stdlib.h>
#include <stdint.h>
#include <math.h>
#include <quadmath.h>
#include <string.h>

#include "nr.h"
#include "nrutil.h"
#include "main.h"
#include "cuba.h"

long idum;

int ndim,VERBOSE;
int comp, nregions, neval, fail,failS;
double integral[NCOMP],error[NCOMP], prob[NCOMP];
double psi, phi , omega, Pi,Z;

int Integral(const int *ndim, const double xx[],
  const int *ncomp, double ff[], void *userdata) 
// float fxn(float pt[], float wgt)
{
#define sig ff[0]

       double b,bp,a,ap,a1,a1p,enp,q2,dd,s1,s2,qa,qb;
       double en,t,tp,f;
//       double Pi = 2.*acos(0.); //,omega,psi,phi;


         en  = xx[0]*omega;
         t   = xx[1]*Pi/2.;
         tp  = xx[2]*Pi/2.;
//         phi = pt[4];
//         psi = pt[5];
         enp = omega - en;
   
//        printf ("%e %e %e %e \n",en,t,tp, Pi);    
 
 /*        if ( phi < 0.) {
	  
             printf ("%e %e \n",phi,psi);
         }
 */        
 //        phi = phi-psi;
  
          if (en > 0.511) { 
              b = sqrt(1-pow((0.511/en),2));
            } else { 
              b=0.;
         }                    
              
         if (enp > 0.511) { 
              bp = sqrt(1-pow((0.511/enp),2));
            } else { 
              bp=0.;  
         }
                       
         qa = en*enp*(1.-b*bp*(sin(t)*sin(tp)*cos(phi)+cos(t)*cos(tp))); 
         qb = omega*en*(b*cos(t)-1.)+omega*enp*(bp*cos(tp)-1)+ pow(0.511,2);
         q2 = -2.*(qa+qb);
  
         if ((q2 == 0.)) {
         printf ("en = %f t = %f tp = %f phi = %f psi = %f \n",en,t,tp,phi,psi);
         printf ("w = %f b = %f bp = %f \n",omega,b,bp);
         printf ("q2 = %f qa = %f qb = %f \n",q2,qa,qb);
         }
  
         a = b*sin(t)/(1.-b*cos(t));
         ap = bp*sin(tp)/(1.-bp*cos(tp));  
         a1 = a*cos(phi+psi);
         a1p = ap*cos(psi);
         dd = sin(t)*sin(tp)*b*bp;
         s1 = (en*enp*dd/ pow(q2,2)/pow(omega,3)) * ( 4.*pow((en*a1p+enp*a1),2) - 
               q2* pow((a1p-a1),2)); 

         s2 = a*ap/omega/pow(q2,2)*(pow(en*b*sin(t),2) 
             + pow(enp*bp*sin(tp),2)
             + 2*en*enp*b*bp*cos(phi)*sin(t)*sin(tp)); 
             

        if (Z!=0.) {f = 1./(1.+pow(111.*sqrt(abs(q2))*pow(Z,(-1./3.)),2));} 
           else
          {f = 0.;}
//         write(*,*) q2,f

         sig = (2./pow((2.*Pi),2))*(s2 - s1)*pow((1.-f),2);
         sig = sig*omega*pow(Pi,2)/4.; // normalización de la integración...
         
         if (sig <0.) {
             printf(" %e %e %e \n", phi, psi, omega);    
             printf ("s2 = %e %e %e %e \n", s2, t, a, a1);
             printf ("s1 = %e %e %e %e \n", s1, tp, ap, a1p);
             sig = 0;
         }                            

      return 0;
      
}


typedef int (*funcptr)(const int *ndim, const double xx[], const int *ncomp, double ff[], void *userdata);


int main(int argc, char *argv[])
{
//  double omegav[] = {5.,6., 7., 8., 9., 10., 15., 20., 30., 40., 50, 60., 70., 80., 90., 100., 200.,300., 400., 500.  };      
//  double omegav[] = {2., 3., 4., 5.,6., 7., 8., 9., 10., 11., 12., 13., 14., 15., 16., 17., 18., 19., 20.,21., 22., 23., 24., 25. };      
//  double omegav[] = { 26., 27., 28., 29.,30., 31., 32., 33., 34., 35., 37.5, 40., 45., 50.  }; // 10 
//  double omegav[] = { 36., 37., 38., 39., 41., 42., 43., 44., 46., 47., 48., 49.,   }; // 10 
  double omegav[] = {60., 70., 80., 90., 100., 200.,300., 400., 600., 700.,800., 900., 1000.  }; // 10 
  char nombre [35],buffer1 [35], buffer2 [35], nombre1 [35], cener[6]; // ,bft [35], diagram [5], UV[2], MI[2];
  char MI[2];
  int init, itmax, ncall, nprn,i,j,k,l,ii,paso=200,ener,eini=0, efin;
  double in, chi2a, sd, lsup, linf,temp,psip; //, psi, phi, omega;
  double res[300][300];
  
 int nener = sizeof(omegav)/sizeof(omegav[0]);   
  
  Pi = 2.*acos(0.); 
  
  FILE *fp, *fphi; 
  
  strcpy(MI,argv[1]);
  strcpy(cener,argv[2]);
  
  if (!strcmp(argv[2],"0")) { 
     efin=nener;
    } else {eini=1;
     efin=2;
     omega=atof(argv[2]);}
     Z=atof(argv[3]);
  printf("%s %s %f %f \n", MI,cener, omega,Z);
//  getc(stdin);

  for (ener=eini;ener<efin;ener++)
   { 
    if (!strcmp(argv[2],"0")) { omega = omegav[ener];}
//    omega = 100.;
    int iomega = omega;
//  Z = 18.;
  
  strcpy(nombre ,"res/dat-");
  strcpy(nombre1,"res/mat-");
 
  strcpy(buffer2, "mc2v1.dat");
  sprintf(buffer1, "%d", iomega);
  strcat(nombre, buffer1);
  strcat(nombre, buffer2);

  strcpy(buffer2, "mc2v1.dat");
  sprintf(buffer1, "%d", iomega);
  strcat(nombre1, buffer1);
  strcat(nombre1, buffer2);
  
  printf("%s \n",nombre); 

  fp = fopen ( nombre, "w" );
    
//  omega = 500.;
  
  linf = 0.;
  lsup = 2.9;
 
   for (i=0; i<=paso;i++) {
//     psi = i*2.*Pi/100.;         

    if (i<100) {
      // phi = (lsup -linf)*i/(paso-100) + linf;
      ii = 0;}
     else {
       linf = 2.9;
       lsup = Pi;
       ii = 100; }
           
       phi = (lsup -linf)*(i-ii)/(paso-100) + linf;
//      fprintf (fphi, " %f \n",phi);
     
     for (j=0; j<=paso;j++) {
          

//          phi = j*(Pi-linfan)/30. + linfan;
//          phi = j*2.*Pi/100.;
          
//          phi = (lsup -linf)*j/500. + linf;
          
//          phi = phi-psi;

          psi = j*Pi/paso;   
//          psi = Pi/2;     

VERBOSE = 0;   
//        Cuba(1,NDIM, NCOMP, diagr[idiagr],inV, error, prob);         
    if (!strcmp(MI,"V")) {
	Vegas(NDIM, NCOMP, Integral, USERDATA, NVEC,
	    EPSREL, EPSABS, VERBOSE, SEED,
	    MINEVAL, MAXEVAL, NSTART, NINCREASE, NBATCH,
	    GRIDNO, STATEFILE, SPIN,
	    &neval, &fail, integral, error, prob);
     }
//          printf (" i %d var_li= %e  int= %e \n", i, var_li, in) ; //)/(p1max -p1min), in); 
//          getc(stdin);
    if (!strcmp(MI,"S")) {
       Suave(NDIM, NCOMP, Integral, USERDATA, NVEC,
             EPSREL, EPSABS, VERBOSE | LAST, SEED,
             MINEVAL, MAXEVAL, NNEW, NMIN, FLATNESS,
             STATEFILE, SPIN,
             &nregions, &neval, &fail, integral, error, prob);
    }

    if (!strcmp(MI,"D")) {
       if (omega>200.) {EPSREL = 1e-8;
	EPSABS = 1e-18;
	}
//       printf("Divonne \n");
       Divonne(NDIM, NCOMP, Integral, USERDATA, NVEC,
               EPSREL, EPSABS, VERBOSE, SEED,
               MINEVAL, MAXEVAL, KEY1, KEY2, KEY3, MAXPASS,
               BORDER, MAXCHISQ, MINDEVIATION,
               NGIVEN, LDXGIVEN, NULL, NEXTRA, NULL,
               STATEFILE, SPIN,
               &nregions, &neval, &fail, integral, error, prob);
   }

    if (!strcmp(MI,"C")) {
	Cuhre(NDIM, NCOMP, Integral, USERDATA, NVEC,
	    EPSREL, EPSABS, VERBOSE | LAST,
	    MINEVAL, MAXEVAL, KEY,
	    STATEFILE, SPIN,
	    &nregions, &neval, &fail, integral, error, prob);
    }
/*          write(80,*) i+1,j+1,in
          ang = i*Pi/30.
          if (i.eq.30) write(81,*) phi,in                      
C	    write(81,*) ang, in
     	    if (i.eq.15) write(82,*) phi,in
*/
//          phi = (lsup -linf)*j/100. + linf;
         
	  if ( phi < 0.) {
	  
             printf ("sali de vegas %e %e %e %i %i\n",phi,psi,integral[0], i, j);
             printf("%e %e \n", lsup, linf);
          }
          res[i][j]=integral[0];
//          psip = psi+phi/2.;
//          if (psip>Pi) {psip=psip-Pi;}
          fprintf (fp, " %f %f %f %i %i %f \n",phi,psip,integral[0],i,j,psi);
//          printf ( " %f %f %f %i %i %f \n",phi,psip,integral[0],i,j,omega); //fpsiphi(omega,phi,psi)

     }              

   }
//      fclose ( fphi );
      fclose ( fp );

     fphi = fopen ( nombre1, "w" );
     
//debería redefinir j para producir el desplazamiento en función de phi

     for (i=0; i<=paso;i++) {
      fprintf (fphi, " %f ",i*Pi/paso);
      for (j=0; j<=paso;j++) {   
//    fprintf (fp, " %e",res[j][i]/res[j][99]);
        fprintf (fphi, "%e ",res[j][i]);
     }
        fprintf (fphi, " \n ");
    }
    fclose ( fphi );
}
   return 0;
      
}


