#include <stdio.h>
#include <math.h>
#include <stdlib.h>
#include <string.h>

int n_data = 300;
double bb = 0.2497;
double bd = 5.16;
double ad = 0.3105;
double ah = 64.3;
int G = 1;
int n = 1000;
double M = 2.325E7;
double e = 2.71828;

double chi_sqrt(double *a,double *b);
double *velocidadCirc(double Mb,double Md,double Mh);
double random_Gauss(void);

int main(void){

	double *r = malloc(n_data*sizeof(double));//Radios en parsec
	double *v = malloc(n_data*sizeof(double));//velocidad en km/s
	double *a = malloc(n_data*sizeof(double));
	double *b = malloc(n_data*sizeof(double));
	double *Mb = malloc(n*sizeof(double));
	double *Md = malloc(n*sizeof(double));
	double *Mh = malloc(n*sizeof(double));
	double *like_w = malloc(n*sizeof(double));
	double *vel_pri = malloc(n_data*sizeof(double));
	Mb[0] = 100.0*random_Gauss()+3500.0;
	Md[0] = 100.0*random_Gauss()+3500.0;
 	Mh[0] = 100.0*random_Gauss()+3500.0;
	double *vel_ini = velocidadCirc(Mb[0], Md[0],Mh[0]);

	//IMPORTAR DATOS 
	FILE *data;
	data = fopen("RV.dat", "r");
	for(int i=0; i<300; i++){
		fscanf(data, "%lf %lf\n", &r[i], &v[i]);
		//printf("%lf %lf\n", r[i], v[i]);
	}
	like_w[0] = chi_sqrt(a, vel_ini);

	for(int i=0; i<n; i++){
    		double Mb_pri = 100.*random_Gauss()+3500.;
    		double Md_pri = 100.*random_Gauss()+3500.;
		double Mh_pri = 100.*random_Gauss()+3500.;
		vel_ini = velocidadCirc( Mb[i], Md[i],Mh[i]);
    		vel_pri = velocidadCirc(Mb_pri, Md_pri, Mh_pri);
    		double like_pri = chi_sqrt(v, vel_pri);
    		double like_ini = chi_sqrt(v, vel_ini);
    
    		double alpha = like_pri/like_ini;
    		if(alpha>=1.){
    			Md[i] = Md_pri;
			Mh[i] = Mh_pri;
        		like_w[i] =  like_pri;
		}
    		else{
       			double beta =drand48();
        		if(beta<=alpha){
            			Mb[i] = Mb_pri;
        			Md[i] = Md_pri;
				Mh[i] = Mh_pri;
        			like_w[i] =  like_pri;
			}
			else{
				Mb[i] = Mb[i];
        			Md[i] = Md[i];
				Mh[i] = Mh[i];
        			like_w[i] =  like_ini;		
			} 
       			
		}
	}
	double num = 0;
	int busc = 0;
	for(int i=0; i<n; i++){
		if(like_w[i]>num){
			num = like_w[i];
			busc = i;
		}	
	}
	//printf("%lf %lf %lf\n", Mb[busc] ,Md[busc] ,Mh[busc]);

	//Modelo
	double *mod = malloc(n_data*sizeof(double));
	for(int i=0; i<n_data; i++){
		mod[i] = ((pow(Mb[busc],1./2.)*r[i])/pow((pow(r[i],2)+pow(bb,2)),3./4.))+((pow(Md[busc],1./2.)*r[i])/pow((pow(r[i],2)+pow(bd+ad,2)),3./4.))+((pow(Mh[busc],1./2.))/pow((pow(r[i],2)+pow(ah,2)),1./4.));
		printf("%lf\n", mod[i]);
	}
	return 0;

}

double chi_sqrt(double *a,double *b){
    double chi = 0.;
    for(int i=0; i<n_data; i++){
    	chi = chi+(pow((a[i]-b[i]),2));
    }
    return pow(e,-1./2.*chi);
}


//la formula del ejercicio
double *velocidadCirc(double Mb, double Md, double Mh){
	double *radio = malloc(n_data*sizeof(double));
	double *v_circ = malloc(n_data*sizeof(double));
	FILE *data1;
	data1 = fopen("RV.dat", "r");
	for(int i=0; i<300; i++){
		fscanf(data1, "%lf\n", &radio[i]);}
	fclose(data1);
	for(int i=0; i<n_data; i++){	
		v_circ[i] = ((pow(Mb,1./2.)*radio[i])/pow((pow(radio[i],2)+pow(bb,2)),3./4.))+((pow(Md,1./2.)*radio[i])/pow((pow(radio[i],2)+pow(bd+ad,2)),3./4.))+((pow(Mh,1./2.))/pow((pow(radio[i],2)+pow(ah,2)),1./4.));
	}	     
  	return v_circ;
}

double random_Gauss(){ 
	static double v1, v2, s;
	static int fase = 0;
	double x;
	if(fase == 0) {
		do {
			double u1 = (double)rand()/RAND_MAX;
			double u2 = (double)rand()/RAND_MAX;
			v1 = 2*u1-1;
			v2 = 2*u2-1;
			s = (v1*v1)+(v2*v2);
		} 
		while(s>=1 || s==0);
			x = v1*sqrt(-2.*log(s)/s);
	} else
		x = v2*sqrt(-2.*log(s)/s);
	fase = 1-fase;

	return x;
}



