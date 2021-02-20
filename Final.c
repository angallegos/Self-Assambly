#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <time.h>



float x[1000], y[1000], z[1000]; 			/*Aquí se van a guardar las partículas*/ 
float auxx,auxy,auxz,x2,y2,z2;
int N=11;						/*Número de partículas*/
float d=1;	
float d2=2;	
float d3=0.02;
float R=2;				/*Radio de la esfera*/
float a;
float theta, phi;
int i,j,k,l;
float r;
int parti;
float yaw, pitch;
float metro,MC;
float paso=0.1;
int rechazos;
int rechazosMC=0;
int MCsteps=1000000;
int npics=100;
int nsample;
int cic;
float enew,eold;
float kb=1;
float epsilon=1;
float tmax=1.01;
float tmin=0.01;
float deltat=0.1;
float t;
float rmin,dista;
int contador=0;
float aux1,aux2;
float Rnew;
float sal;
float cut=1E8;
float acept,aceptacion;
float dr=0.1;
float deltaradio=0.01;
double urandom(){
    //RAND_MAX is a constant defined in the rand module
    return  (double)rand() / (double)RAND_MAX;
}

double random_sphere_coord(){
    return  2 * (urandom() - 0.5);
}	
int energia(){
	int energiaf=0;
	for(i=1;i<N;i++){
		for(j=0;j<i;j++){
			r=sqrt(pow(x[i]-x[j],2)+pow(y[i]-y[j],2)+pow(z[i]-z[j],2));
				if(r<d2){
					energiaf=energiaf+1;
				}
		}
	}
	return energiaf;
}							
float posiciones_inciales(){
	double random_sphere_coord();
	a=random_sphere_coord();
	theta=acos(a);		/*Ángulos*/
	phi=drand48()*2*M_PI;	
	x[0]=R*sin(theta)*cos(phi);
	y[0]=R*sin(theta)*sin(phi);
	z[0]=R*cos(theta);
	/*printf("Pudo poner la partícula %i\n",i+1);*/
	for(i=1;i<N;i++){
		double random_sphere_coord();
		lanza: a=random_sphere_coord();
		theta=acos(a);		/*Ángulos*/
		phi=drand48()*2*M_PI;	
		x[i]=R*sin(theta)*cos(phi), y[i]=R*sin(theta)*sin(phi), z[i]=R*cos(theta);
	
		for(j=0;j<i;j++){
			r=sqrt(pow(x[i]-x[j],2)+pow(y[i]-y[j],2)+pow(z[i]-z[j],2));
			if(r<d){
				goto lanza;
			}
		}
		/*printf("Pudo poner la partícula %i\n",i+1);*/
	}	
	printf("La energía de la configuración inicial es %i\n",energia());
	printf("Ya acomodó las partículas\n");
		
}

int mover(){
	parti=rand()%N;
	auxx=x[parti];
	auxy=y[parti];
	auxz=z[parti];	
	yaw=paso*random_sphere_coord();
	pitch=paso*random_sphere_coord();
	x[parti]=((auxx*cos(yaw)*cos(pitch))-(auxy*cos(yaw)*sin(pitch))+(auxz*sin(yaw)));
	y[parti]=((auxx*sin(pitch))+(auxy*cos(pitch)));
	z[parti]=(-(auxx*sin(yaw)*cos(pitch))+(auxy*sin(yaw)*sin(pitch))+(auxz*cos(yaw)));
	for(i=0;i<N;i++){
		if(i==parti){
			continue; 	
		}
		r=sqrt(pow(x[parti]-x[i],2)+pow(y[parti]-y[i],2)+pow(z[parti]-z[i],2));
		if(r<d){
			x[parti]=auxx;
			y[parti]=auxy;
			z[parti]=auxz;
			return 0;
		}

	}

}
float potencial(float x2,float y2, float z2){
	float energy=0;
		for(i=0;i<N;i++){			
			if(i==parti){
				continue; 	
			}
			r=sqrt(pow(x2-x[i],2)+pow(y2-y[i],2)+pow(z2-z[i],2));
			if(r<d2){
				energy++;			
			}
			if(r-d<d3){
				energy=energy+(1/pow(r-d,2));
			}
		}	
	return energy;
}
float distancia_minima(){
	int ii,jj;
	float dist;
	float rmin=2*R;
	for(ii=1;ii<N;ii++){
		for(jj=0;jj<ii;jj++){
			dist=sqrt(pow(x[ii]-x[jj],2)+pow(y[ii]-y[jj],2)+pow(z[ii]-z[jj],2));
			/*printf("La distancia entre la partícula %i y la %i es: %f\n",ii,jj,dist);*/
			if(dist<rmin){
			rmin=dist;
			}
		}
	}
	return rmin;
}
float min(){
	return (deltaradio<aux1?deltaradio:aux1); //Operador ternario
}
int escala(){
	rmin=distancia_minima();
	/*printf("la distancia mínima es %f\n",rmin);*/
	aux1=R*(1-(d/rmin));
	/*printf("deltaradio es %f y el otro número es %f\n",deltaradio,aux1);*/
	aux2=min();
	/*printf("El mínimo entre %f y %f es %f\n",deltaradio,aux1,aux2);*/
	Rnew=R-aux2;
	/*printf("El nuevo radio es %f\n",Rnew);*/
	for(k=0;k<N;k++){
		x[k]=(Rnew/R)*x[k];
		y[k]=(Rnew/R)*y[k];
		z[k]=(Rnew/R)*z[k];
	}
	R=Rnew;
}

float potencial_final(){
	float energiaf=0;
	for(i=1;i<N;i++){
		for(j=0;j<i;j++){
			r=sqrt(pow(x[i]-x[j],2)+pow(y[i]-y[j],2)+pow(z[i]-z[j],2));
				if(r<d2){
					energiaf=energiaf+1;
				}
				if(r-d<d3){
					energiaf=energiaf+(1/pow(r-d,2));
				}
		}
	}
	return energiaf;
}
int imprime_parametros(){
	FILE *pa;
	char salida[]="dat11.dat";
	pa=fopen(salida,"w");
	fprintf(pa,"%i\n",N);
	fprintf(pa,"%f\n",d);
	fprintf(pa,"%f\n",d2);
	fprintf(pa,"%i\n",MCsteps);
	fprintf(pa,"%f\n",R);
	fprintf(pa,"%f\n",tmax);
	fprintf(pa,"%f\n",tmin);
	fprintf(pa,"%f\n",deltat);
	fprintf(pa,"%f\n",paso);
	fclose(pa);
}
int escribe(){
	FILE *pa;
	char salida[]="dat11.dat";
	pa=fopen(salida,"a+");
	fprintf(pa,"%f\n",R);
	//fprintf(pa,"%i\n",energia());
	for(l=0;l<N;l++){
		fprintf(pa,"%f	%f	%f\n",x[l],y[l],z[l]);
	}
	fclose(pa);
}

int main(){
	srand(time(NULL));
	srand48(time(NULL));
	printf("La energía de corte es: %f\n",cut);
	imprime_parametros();
	posiciones_inciales();
	printf("\n");
	escribe();
	/*getchar();*/
	printf("Radio: %f\n",R);
	do {
		contador++;
		rechazos=0;
		for(cic=0;cic<=10;cic++){
			t=tmax-(cic*deltat);	
			for(j=0;j<MCsteps;j++){
				if(mover()==0){
					rechazos++;
					rechazosMC++;
				}	
				eold=potencial(auxx,auxy,auxz)*epsilon;
				/*printf("La energía antes de mover la partícula es %f\n",eold);*/
				enew=potencial(x[parti],y[parti],z[parti])*epsilon;
				/*printf("La energía después de mover la partícula es %f\n",enew);*/
				/*printf("%i\n",j);*/
				metro=drand48();
				MC=exp(-(enew-eold)/(kb*t));
				if(metro>MC){
					x[parti]=auxx;
					y[parti]=auxy;
					z[parti]=auxz;
					rechazos++;
					rechazosMC++;
				}
				nsample=MCsteps/npics;
				if(j % nsample == 0){
					aceptacion=(1.0-((float)rechazosMC/(float)nsample))*100.0;
					/*printf("El paso es %f %f\n",aceptacion,paso);*/
					/*getchar();*/
					if(aceptacion<48){
						paso=0.5*paso;
					}
					if(aceptacion>52){
						paso=2.0*paso;
					}	
					rechazosMC=0;
				}
				/*printf("%i\n",j);*/
			}
			/*printf("%f\n",t);*/
			/*printf("El paso al final de la temperatura %f es %f\n",t,paso);*/
		}
		printf("Reescalamiento %i\n",contador);
		/*printf("La cantidad de rechazos fue: %i\n",rechazos);*/
		acept=(1-((float)rechazos/(MCsteps*cic)))*100;
		printf("El porcentaje de aceptación es: %f%\n", acept);
		/*printf("El 'paso' del usado fue: %f\n",paso);*/
		sal=potencial_final();
		printf("El radio usado fue %f\n",R);
		printf("La energía de la configuración es %i:\n", energia());
		printf("La 'energía' de la configuración es: %f\n",sal);
		escribe();		
		escala();
		printf("La distancia mínima entre partículas es: %f\n",rmin);
		printf("Reducción: %f\n",aux2);
		printf("El radio nuevo es: %f\n",R);
		printf("\n");
		
	} while(sal<cut || rmin>1);
	printf("La energía de la configuración final es: %i\n",energia());
	return 0;
}

