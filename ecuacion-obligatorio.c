//TRABAJO OBLIGATORIO: ECUACION DE SCHRODINGER

# include <stdio.h>
# include <math.h>
# include "complex.h" 


# define PI 3.1415926535
# define N 2000
# define nciclos 80
# define lambda 0.5
# define T 3008



//Definimos las funciones que vamos a utilizar. 
void norma(fcomplex phi[N+1]);
double posicion(fcomplex phi[N+1]);
double momento(fcomplex phi[N+1]);

int main()
{

	int j,k,n;  //Variables auxiliares
    fcomplex x,y,z,c,num,den;  //Variables auxiliares
	double sum;  //Variables auxiliares

	//Definimos varias variables 

	double k0=1.0*2*PI*nciclos/N; //k0 reescalado (k0 tilde)
	double s=1.0/(4*pow(k0,2)); // s reescalado (s tilde)
	fcomplex m=Complex(2,-(2*1.0)/s);  //m=2-2i/s
	fcomplex i=Complex(0,1); 

	fcomplex V[N+1];
	fcomplex alfa[N];
	fcomplex phi[N+1];
	fcomplex beta[N];
	fcomplex chi[N+1];
	
 
	FILE *f1,*f2,*f3,*f4,*f5;

	//Abrimos ficheros para guardar los datos

	f1=fopen("ecuacion.txt","w"); 
	f2=fopen("potencial.txt","w"); 
	f3=fopen("norma.txt","w");
	f4=fopen("posicion.txt","w");
	f5=fopen("momento.txt","w");

	
	
	//Calculamos el potencial V 


	for(j=0;j<N+1;j++)
	{	
		if(j<2*N/5 || j>3*N/5)
		{
			V[j]=Complex(0,0);
		}
		else 
		{
			V[j]=Complex(lambda*pow(k0,2),0);
		}

		
	}

	//Calculamos el vector alfa que no depende del tiempo

	alfa[N-1]=Complex(0,0);	
	
	for(j=0;j<N-1;j++)
	{
		x=Cadd(m,V[j]);
		y=Csub(x,alfa[N-1-j]);
		alfa[N-2-j]=Cdiv(Complex(1,0),y);
	
	}
	

	//Calculamos la función de onda inicial suponiendo las condiciones de contorno

	phi[0]=Complex(0,0);
	phi[N]=Complex(0,0);

	for(j=1;j<N;j++)
	{
		phi[j]=RCmul(exp(-8*pow((4*j-N),2)/(pow(N,2))),Cgauss(k0*j,1));
	}

	norma(phi);  //normalizamos la función de onda con la función "norma"

	for(j=0;j<N+1;j++)
	{
		fprintf(f1,"%i\t%lf\n",j,pow(Cabs(phi[j]),2));   //imprimimos los valores de la función de onda en el fichero "ecuacion.txt"
	}

	fprintf(f1,"\n\n");

	//Calculamos el valor esperado de la posición en el instante inicial y lo guardamos en el fichero "posicion.txt"

	fprintf(f4,"%i\t%lf\n",0,posicion(phi));


	//Calculamos el valor esperado del momento en el instante inicial y lo guardamos en el fichero "momento.txt"
	
	fprintf(f5,"%i\t%lf\n",0,momento(phi));


	
    //Calculamos la evolución temporal de la función de onda.

	for(n=0;n<T;n++)
	{
		//Imponemos las condiciones de contorno para beta y chi

		beta[N-1]=Complex(0,0); 
		chi[0]=Complex(0,0);
		chi[N]=Complex(0,0);


		for(j=0;j<N-1;j++)
		{
			x=RCmul(-4/s,phi[N-1-j]);
			y=Cmul(i,x);
			num=Cadd(beta[N-1-j],y);

			c=Cadd(m,V[j]);
			den=Csub(c,alfa[N-1-j]);
			
			beta[N-2-j]=Cdiv(num,den);	

		}

		for(j=0;j<N-1;j++)
		{
			chi[j+1]=Cadd(Cmul(alfa[j],chi[j]),beta[j]);
		
		}

		for(j=0;j<N+1;j++)
		{
			phi[j]=Csub(chi[j],phi[j]);
		
		}

		norma(phi);

		sum=0;
		
		for(j=0;j<N+1;j++)
		{	
			sum=sum+pow(Cabs(phi[j]),2);   //calculamos la norma
			fprintf(f1,"%i\t%lf\n",j,pow(Cabs(phi[j]),2)); //imprimimos la función de onda en el instante n en el fichero "ecuacion.txt"
			fprintf(f2,"%i\t%lf\n",j,Cabs(V[j]));  //imprimimos el potencial en el fichero "potencial.txt"
		}
		
		fprintf(f3,"%lf\n",sum);   
		fprintf(f1,"\n\n");
		fprintf(f2,"\n\n");
		fprintf(f4,"%i\t%lf\n",n+1,posicion(phi));
		fprintf(f5,"%i\t%lf\n",n+1,momento(phi));
		
		
	}
	
	fclose(f1);	
	fclose(f2);
	fclose(f3);
	fclose(f4);
	fclose(f5);	
	
	return 0;
}

//La función "normal" devuelve la función de onda normalizada. 
 
void norma(fcomplex phi[N+1])
{
	int j;
	double k;

	for(j=0;j<N+1;j++)
	{
		k+=pow(Cabs(phi[j]),2);
	}

	for(j=0;j<N+1;j++)
	{
		phi[j]=RCmul(pow(sqrt(k),-1),phi[j]);

	}

	return;
}

double posicion(fcomplex phi[N+1])

{
	int j;
	double x=0;

	for(j=1;j<N;j++)
	{
		x=x+j*pow(Cabs(phi[j]),2);

	}

	return x;
}

double momento(fcomplex phi[N+1])

{
	int j;
	double p=0;
	fcomplex i=Complex(0,-1);
	fcomplex x1,x2,x3,x4,x5;

	for(j=0;j<N;j++)
	{
		x1=Csub(phi[j+1],phi[j]);
		x2=Cadd(Conjg(phi[j+1]),Conjg(phi[j]));
		x3=RCmul(0.5,x2);
		x4=Cmul(x1,x3);
		x5=Cmul(i,x4);
		p=p+x5.r;

	}

	return p;

}

