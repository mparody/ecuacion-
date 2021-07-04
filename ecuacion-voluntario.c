//VOLUNTARIO: ECUACION DE SCHRODINGER-COEFICIENTE DE TRANSMISIÓN. 

# include <stdio.h>
# include <math.h>
# include "complex.h" 
# include "gsl_rng.h"

# define pi 3.1415926535
# define N 2000  	//Número de puntos en los que discretizamos la función de onda.
# define nciclos 80 //Número de oscilaciones completas que la función de onda tiene sobre la red.
# define lambda 0.5	//lambda: parámetro que influye en la altura de la barrera de potencial.
# define nD 600 	//nD: número de pasos para los que dejamos evolucionar el sistema. 



gsl_rng *tau; 

//Definimos las funciones que vamos a utilizar. 

void potencial(fcomplex V[N+1],double k0);
void f_onda_inicial(fcomplex phi[N+1],double k0);
void Alfa(fcomplex alfa[N],fcomplex V[N+1],double s);
void Beta(fcomplex beta[N],fcomplex alfa[N],fcomplex phi[N+1],fcomplex V[N+1],double s);
void norma(fcomplex phi[N+1]);
double posicion(fcomplex phi[N+1]);
double momento(fcomplex phi[N+1]);


int main()
{

	int j,k,n,cont,t;  //Variables auxiliares.
	int D,I; //D,I: valen 1 ó 0, dependiendo de si la partícula ha sido detectada o no respectivamente por el detector de la derecha/izquierda.
	int m_T=0; //contador de las detecciones a la derecha.
	int m=1000; //número de simulaciones.
	double PD, PI; //PD: probabilidad de que el detector de la derecha detecte la partícula. PI: análogamente con el detector de la izquierda. 
	double K; //coef de transmisión.
	double x; //número aleatorio.
	extern gsl_rng *tau; //Puntero al estado del número aleatorio.
    int semilla=8293401; //Semilla del generador de números aleatorios.

	
	fcomplex V[N+1];
	fcomplex alfa[N];
	fcomplex phi[N+1];
    fcomplex chi[N+1];
    fcomplex beta[N];
	
	double k0=1.0*2*pi*nciclos/N; //k0 reescalado (k0 tilde)
	double s=1.0/(4*pow(k0,2)); // s reescalado (s tilde)
	
	
	FILE *f1,*f2,*f3;  //definición de ficheros.

	//Abrimos ficheros para guardar los datos.

	f1=fopen("posicion.txt","w"); 
	f2=fopen("contador.txt","w");
	f3=fopen("momento.txt","w");

	fprintf(f1,"Tiempo \t <x>\n");
	fprintf(f2,"Nº simulación \t Nº de veces que evoluciona la partícula\n");
	fprintf(f3,"Tiempo \t <p>\n");
	
	
	tau=gsl_rng_alloc(gsl_rng_taus); //Inicializamos el puntero.
    gsl_rng_set(tau,semilla); //Inicializamos la semilla.

	potencial(V,k0); //Calculamos el potencial V mediante la función "potencial". 
	

	Alfa(alfa,V,s); //Calculamos el vector alpha mediante la función "Alfa".


 	//A continuación, colocamos dos detectores de anchura N/5 actuando a derecha y a izquierda del detector. Queremos contar el número de veces que la partícula atraviesa la barrera de potencial, es decir, que es detectada por el detector de la derecha con el fin de calcular el coeficiente de transmisión. 
	
	t=0; //contador para contar el número de veces que hemos calculado el valor esperado de la posición y el momento. 
	

	for(k=0;k<m;k++)
	{
		cont=0;  //contador para contar el número de veces que evoluciona la partícula antes de ser detectada en cada simulación.

		f_onda_inicial(phi,k0); //Calculamos la función de onda inicial normalizada considerando las condiciones de contorno mediante la función "f_onda_inicial".
	
		
		//Calculamos el valor esperado de la posición en el instante inicial mediante la función "posicion" e imprimimos su valor en el archivo "posicion.txt". Solo lo calculamos para la cuarta simulación, por eso escribimos el "if(k==4)". En el informe se detallará por qué se elige la interacción 4. Ídem para el momento. 

		if(k==4)    
		{
			fprintf(f1,"%i\t%lf\n",t,posicion(phi)); 
			fprintf(f3,"%i\t%lf\n",t,momento(phi));

		}	

		D=0;  
		I=0; 
		
		while(D==0 && I==0)    //Este bucle while nos sirve para hacer evolucionar el sistema tantas veces como sea necesario hasta que la partícula sea detectada por 									uno	de los dos detectores.
		{
			
			//Dejamos evolucionar la función de onda un número de pasos nD

            chi[0]=Complex(0,0);  //condición de contorno de chi.
	        chi[N]=Complex(0,0);  //condición de contorno de chi.

	        for(n=0;n<nD;n++)
	        {		

		        Beta(beta,alfa,phi,V,s); //Calculamos beta mediante la función "Beta".
		
		        //Calculamos chi a partir de alpha y beta.

	    	    for(j=0;j<N-1;j++)
		        {
		        	chi[j+1]=Cadd(Cmul(alfa[j],chi[j]),beta[j]);
		
		        }

		       //Calculamos la función de onda en el siguiente instante n a partir de chi. 

		        for(j=0;j<N+1;j++)
		        {
		    	    phi[j]=Csub(chi[j],phi[j]);
		        }
		
		        //Normalizamos la función de onda 
		
		        norma(phi);

		        //Calculamos el valor esperado de la posición en el instante n mediante la función "posicion" e imprimimos su valor en el archivo "posicion.txt". Solo lo calculamos para la cuarta simulación, por eso escribimos el "if(k==4)". Ídem para el momento. 
				
				
		
		        if(k==4)
				{
					t=t+1;
					fprintf(f1,"%i\t%lf\n",t,posicion(phi));
					fprintf(f3,"%i\t%lf\n",t,momento(phi)); 

				}	
	
	        }

			PD=0; //igualamos a 0 la probabilidad de que sea detectada a la derecha.
			PI=0; //igualamos a 0 la probabilidad de que sea detectada a la izquierda.
			
			for(j=1.0*4*N/5;j<N+1;j++)
			{
				PD=PD+pow(Cabs(phi[j]),2);   //calculamos la probabilidad de que sea detectada a la derecha.
			}

			x=gsl_rng_uniform(tau);   //generamos un número aleatorio

			if(x<PD)      //si se da esta condición, la partícula ha sido detectada a la derecha de la barrera, por lo que ha sido transmitida. Sumamos 1 a m_T y hacemos 								D=1 para acabar el bucle while.
			{
				m_T=m_T+1;
				D=1;			
			}

			if(x>=PD)   //si no se da, proyectamos la función de onda y calculamos el valor esperado de la posición y el momento en la cuarta simulación.
			{
				for(j=1.0*4*N/5;j<N+1;j++)
				{
					phi[j]=Complex(0,0);
				}

				norma(phi);
				
				
				if(k==4)
		        {
				   t=t+1;
			       fprintf(f1,"%i\t%lf\n",t,posicion(phi));
				   fprintf(f3,"%i\t%lf\n",t,momento(phi));
    
		        }
                

				for(j=0;j<=(N*1.0)/5;j++)    //calculamos la probabilidad de que sea detectada a la izquierda. 
				{
					PI=PI+pow(Cabs(phi[j]),2);  
				}
				
				x=gsl_rng_uniform(tau);   //generamos un número aleatorio.
			
				if(x<PI)   // si se da esta condición, la partícula ha sido detectada a la izquierda de la barrera. Hacemos I=1 para acabar el bucle while.
				{
					I=1;
				}
	
				if(x>=PI)    //si no se da, proyectamos la función de onda  y calculamos el valor esperado de la posición y el momento en la cuarta simulación. Volvemos 								   a repetir el ciclo while. 
				{
					for(j=0;j<=1.0*N/5;j++)
					{
						phi[j]=Complex(0,0);
					}
		
					norma(phi);
					
					
					if(k==4)
		        	{
						 t=t+1;
			       		 fprintf(f1,"%i\t%lf\n",t,posicion(phi));
						 fprintf(f3,"%i\t%lf\n",t,momento(phi));
    
		        	}
                    				
				}
			
			}
			cont++;
		}

		fprintf(f2,"%i\t%i\n",k,cont);  //imprimimos en el fichero "contador.txt" el valor del contador para cada simulación. 
	}

	K=(m_T*1.0)/m;

	printf("Para los siguientes parámetros:\n N=%i \n n_ciclos=%i \n lambda=%lf \n n_D=%i \n\n El coeficiente de transmisión obtenido ha sido %lf.\n",N,nciclos,lambda,nD,K);

	fclose(f1);	
	fclose(f2);
	fclose(f3);
	
	
	return 0;
}

//La función "potencial" calcula la barrera de potencial.

void potencial(fcomplex V[N+1],double k0)
{
	int j;

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
	return;
}


//La función "norma" devuelve la función de onda normalizada. 
 
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

//La función "f_onda_inicial" calcula la función de onda inicial en el instante n=0.

void f_onda_inicial(fcomplex phi[N+1],double k0)
{
	int j;
	double k=0.;
	

	phi[0]=Complex(0,0);
	phi[N]=Complex(0,0);

	for(j=1;j<N;j++)
	{
		phi[j]=RCmul(exp(-8*pow((4*j-N),2)/(pow(N,2))),Cgauss(k0*j,1));
	}

	norma(phi);

	return;
}

//La función "Alfa" calcula el vector alfa auxiliar. 

void Alfa(fcomplex alfa[N], fcomplex V[N+1],double s)
{	
	int j;
	fcomplex x,y;
	fcomplex m=Complex(2,-(2*1.0)/s);


	alfa[N-1]=Complex(0,0);	
	
	for(j=0;j<N-1;j++)
	{
		x=Cadd(m,V[j]);
		y=Csub(x,alfa[N-1-j]);
		alfa[N-2-j]=Cdiv(Complex(1,0),y);
	
	}
	
	return;
}
	
//La función "beta" calcula el vector beta auxiliar.

void Beta(fcomplex beta[N],fcomplex alfa[N],fcomplex phi[N+1],fcomplex V[N+1],double s)
{
	//Definimos variables auxiliares. 

	int j;
	fcomplex x,y,z,num,den;
	fcomplex m=Complex(2,-(2*1.0)/s); 
	fcomplex i=Complex(0,1);

	beta[N-1]=Complex(0,0);  //condición inicial.

	for(j=0;j<N-1;j++)		
	{
		x=RCmul(-4/s,phi[N-1-j]);
		y=Cmul(i,x);
		num=Cadd(beta[N-1-j],y);
		z=Cadd(m,V[j]);
		den=Csub(z,alfa[N-1-j]);
		beta[N-2-j]=Cdiv(num,den);	

	}

	return;
}

//La función "posicion" calcula el valor esperado del observable "posición".

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




