#include <stdio.h>
#include <math.h>

/*
CODIGO DESARROLLADO EN CLASE DE MODELADO MOLECULAR
DATE/FECHA: OCT 2019

SIMULACION DE DINAMICA MODELUCLAR PARA 64 PARTIDULAS
DE ARGON EN UN CUBO.

SE UTILIZA MECANICA NEWTONEANA Y MECANICA ESTADISTICA CLASICA
*/


/*CONSTANTES*/

/*N es el numero de particulas*/
#define N 64
/*L es la longitud de la caja; proporcional a NL*/
#define L 20.0
/*NL es aproximadamente la raiz cubica de N*/
#define NL 4
/*Nt es el tiempo total de simulacion*/
#define Nt 25000
/*PrintTime es cada cuantos datos se imprimiran las energias*/
#define PrintTime 10
/*Constante de Boltzmann en UE/K*/
#define KB 0.8317


/*V A R I A B L E S*/

/*Temperatura final*/
double T =100.0;
/*Incremento temporal de la simulacion*/
double dt = 0.001; /*Tiempo en picosegundos*/
/*Masa*/
double m = 39.948; /*Unidades en umas*/
/*Sigma*/
double sigma=3.4;/*(3.4E-10)/(1/1E-10); Unidades personales: 1UL = 1.0E-10 m (Unidades de Longitud)*/
/*Epsilon*/
double epsilon=(1.65E-21)*(1/1.66E-23); /*Las unidades son personales: 1UE = 1.66E-23 J (Unidades de Energia)*/
/*Posiciones de las N particulas*/
double rx[N];
double ry[N];
double rz[N];
/*Velocidades de las N particulas*/
double vx[N];
double vy[N];
double vz[N];
/*Nuevas posiciones y velocidades*/


/*LLAMADO DE FUNCIONES*/

/*Potencial de Lennard-Jones*/
double Upot(double r);
/*Funcion de asignacion de valores iniciales*/
void iniciales(void);
/*Funcion de calculo de fuerzas*/
void calculaFij(double r, int i, int j, double Fij[], int ipx, int ipy, int ipz); /*A la funcion le pasaremos dos enteros y un arreglo*/


/*FUNCION PRINCIPAL*/
void main(void)
{
    int i, j, t, ix, iy, iz, ipx, ipy, ipz, tau, icero;
    double Ui, r, Ecin, Epot, Fij[3], Fi[3], Ecinprom, Epotprom;
    double T_t, lambda, doubletau;
    double newrx[N],newry[N],newrz[N];
    double newvx[N],newvy[N],newvz[N];

    tau=40;

    printf("\n\t\t\t\tDINAMICA MOLECULAR, TAU=%i\n\n", tau);
    printf("TIEMPO\t\tENERGIA CINETICA\tENERGIA POTENCIAL\tENERGIA TOTAL\t\tENERGIA FINAL\t\tECINPROM\t\tEPOTPROM\n");

    iniciales(); /*Se llama a la funcion "iniciales"*/
    /*Tau es el tiempo necesario para que el termometro haga efecto*/
    doubletau = (double)tau;
    doubletau = doubletau * dt;
    Ecinprom=0.0; /*La energia cinetica promedio al final de las iteraciones*/
    Epotprom=0.0;/*La energia potencial promedio al final de las iteraciones*/
    /*DINAMICA MOLECULAR PARA EL SISTEMA*/
    for (t=1; t<=Nt; t++)
    {
        Ecin=0.0; /*La energia cinetica tras cada iteracion de tiempo se resetea a cero*/
        Epot=0.0; /*La energia potencial tras cada iteracion de tiempo se resetea a cero*/

        for (i=0; i<N; i++) /*Contador para todas las particulas*/
        {
            /*Calculo de las fuerzas*/
            Fi[0]=0.0;
            Fi[1]=0.0;
            Fi[2]=0.0;
            Ui = 0.0;

            /*Contadores relacionados con las cajas aleda�as a la principal*/
            for (ipx=-1; ipx<=1; ipx++)
            {
                for (ipy=-1; ipy<=1; ipy++)
                {
                    for (ipz=-1; ipz<=1; ipz++)
                    {
                       for(j=0; j<N; j++) /*Contador para todas las particulas menos la i-esima*/
                        {
                            if (j!=i) /*Mientras "j" sea diferente a "i"*/
                            {
                                /*Calculo del potencial en i debido a j*/
                                /*r = (rx[j]-rx[i])*(rx[j]-rx[i])+(ry[j]-ry[i])*(ry[j]-ry[i])+(rz[j]-rz[i])*(rz[j]-rz[i]);
                                r = sqrt(r); Aproximacion a una caja*/
                                /*Aproximacion a 27 cajas*/
                                /*Distancia entre la particula "i" y la "j"*/
                                r = (rx[j] + ipx*L - rx[i])*(rx[j] + ipx*L - rx[i]);
                                r = r + (ry[j] + ipy*L - ry[i])*(ry[j] + ipy*L - ry[i]);
                                r = r + (rz[j] + ipz*L - rz[i])*(rz[j] + ipz*L - rz[i]);
                                r = sqrt(r);
                                /*Calculo de la fuerza j sobre i*/
                                calculaFij(r, i, j, Fij, ipx, ipy, ipz); /*Cuando se pasa un arreglo, solo se nombra*/
                                Fi[0]=Fi[0]+Fij[0]; /*Se suman las fuerzas en "X"*/
                                Fi[1]=Fi[1]+Fij[1]; /*Se suman las fuerzas en "Y"*/
                                Fi[2]=Fi[2]+Fij[2]; /*Se suman las fuerzas en "Z"*/
                                /*Calculo del potencial en i debido a j*/
                                Ui = Ui + Upot(r);
                            }
                        }
                    }
                }
            }

            /*Calculo de las nuevas velocidades*/
            newvx[i]=vx[i]+Fi[0]/m*dt;
            newvy[i]=vy[i]+Fi[1]/m*dt;
            newvz[i]=vz[i]+Fi[2]/m*dt;
            /*Calculo de las nuevas posiciones (con primera aproximacion)
            newrx[i]=rx[i]+vx[i]*dt;
            newry[i]=ry[i]+vy[i]*dt;
            newrz[i]=rz[i]+vz[i]*dt;*/
            /*Calculo de las nuevas posiciones (con segunda aproximacion)*/
            newrx[i]=rx[i]+vx[i]*dt + 0.5*Fi[0]/m*dt*dt;
            newry[i]=ry[i]+vy[i]*dt + 0.5*Fi[1]/m*dt*dt;
            newrz[i]=rz[i]+vz[i]*dt + 0.5*Fi[2]/m*dt*dt;
            /*Se calcula las velocidades actuales*/
            Ecin = Ecin + 0.5*m*(vx[i]*vx[i]+vy[i]*vy[i]+vz[i]*vz[i]);
            Epot = Epot + 0.5*Ui;
            /*printf("x=%lf, y=%lf, z=%lf\n",rx[i],ry[i],rz[i]);*/
            /*printf("Fx=%e, Fy=%e, Fz=%e\n\n", Fi[0], Fi[1], Fi[2]);*/
        }
       /* printf("%lf\n", Ecin);*/
       /* printf(" \t %lf \t %lf \t %lf\n",Ecin, Epot, Ecin+Epot);*/

       /*Calculo de la temperatura en funcion del tiempo*/
       if ((t % tau) == 0)
       {
           T_t = 2.0*Ecin/(3.0*N*KB);
           lambda = 1 + dt/doubletau * (T/T_t-1.0);
           lambda = sqrt(lambda);

           for(i=0; i<N; i++)
           {
               newvx[i] = newvx[i] * lambda;
               newvy[i] = newvy[i] * lambda;
               newvz[i] = newvz[i] * lambda;
           }
       }


        /*Aseguramiento de particulas dentro de la caja*/
        if (newrx[i] < 0.0) newrx[i] = newrx[i] + L;
        if (newry[i] < 0.0) newry[i] = newry[i] + L;
        if (newrz[i] < 0.0) newrz[i] = newrz[i] + L;
        if (newrx[i] > L) newrx[i] = newrx[i] - L;
        if (newry[i] > L) newry[i] = newry[i] - L;
        if (newrz[i] > L) newrz[i] = newrz[i] - L;

        /*Impresion de valores cada PrintTime valores*/
        if ((t)%PrintTime==0)
        {
            printf("%i \t\t %lf \t\t %lf \t\t %lf \t\t %lf \t\t",t,Ecin, Epot, Ecin+Epot, 3.0/2.0*N*KB*T);
            printf("%lf \t\t %lf\n", Ecinprom, Epotprom);
        }



/*Actualizacion de posiciones y velocidades*/
        for (i=0; i<N; i++)
        {
            vx[i]=newvx[i];
            vy[i]=newvy[i];
            vz[i]=newvz[i];

            rx[i]=newrx[i];
            ry[i]=newry[i];
            rz[i]=newrz[i];
        }
        icero = 2000;
        if (t > icero)
        {
            Ecinprom = Ecinprom + Ecin;
            Epotprom = Epotprom + Epot;
        }
    }
    Ecinprom = Ecinprom/((double)(Nt - icero));
    Epotprom = Epotprom/((double)(Nt - icero));
    printf("ENERGIA CINETICA PROMEDIO: %lf\n", Ecinprom);
    printf("ENERGIA POTENCIAL PROMEDIO: %lf\n", Epotprom);
    printf("ENERGIA TOTAL PROMEDIO: %lf\n", Ecinprom+Epotprom);
    printf("TEMPERATURA PROMEDIO: %lf", 2.0/3.0*Ecinprom/(N*KB));
    return(0);
}

/*FUNCION DE POTENCIAL*/
double Upot(double r)
{
    return(4.0*epsilon*(pow(sigma/r,12)-pow(sigma/r,6)));
}

/*CALCULO DE FUERZAS SOBRE LA PARTICULA "i" DEBIDO A LA "j"*/
void calculaFij(double r, int i, int j, double Fij[], int ipx, int ipy, int ipz)
{
    double fij;
    /*Valor absoluto de la fuerza ejercida por "j" sobre "i"*/
    fij=-24.0*(epsilon/(sigma*sigma))*pow(sigma/r,8)*(2.0*pow(sigma/r,6)-1.0);
    /*Componentes de la fuerza sobre "i"*/
    Fij[0]=fij*(rx[j] + ipx*L - rx[i]); /*Fuerza ejercida por "j" sobre "i" en "x"*/
    Fij[1]=fij*(ry[j] + ipy*L - ry[i]); /*Fuerza ejercida por "j" sobre "i" en "y"*/
    Fij[2]=fij*(rz[j] + ipz*L - rz[i]); /*Fuerza ejercida por "j" sobre "i" en "z"*/

    /*printf("r%i%i=%e, f%i%i=%e\n", i+1,j+1,r,i+1,j+1,fij);*/
    /*printf("F%i%ix=%e, F%i%iy=%e, F%i%iz=%e\n", i+1, j+1,Fij[0], i+1, j+1,Fij[1], i+1, j+1,Fij[2]);*/
}

/*POSICIONES INICIALES DE LAS PARTICULAS*/
void iniciales(void)
{
    int i, ix, iy, iz;
    for(i=0; i<N; i++)
    {
        vx[i]=0.0;
        vy[i]=0.0;
        vz[i]=0.0;
    }

    /*Acomodo homogeneo de las particulas*/
    for (iz=0; iz<NL; iz++)
    {
        for (iy=0; iy<NL; iy++)
        {
            for (ix=0; ix<NL; ix++)
            {
            rx[iz*NL*NL+iy*NL+ix] = L /(NL+1.0)*(ix+1);
            ry[iz*NL*NL+iy*NL+ix] = L /(NL+1.0)*(iy+1);
            rz[iz*NL*NL+iy*NL+ix] = L /(NL+1.0)*(iz+1);
            }
        }
    }
}
