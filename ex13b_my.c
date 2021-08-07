//  Linear parabolic equation:
//
//  du/dt = (d/dx)(k(x)du/dx) + f(x,t), xa < x < xb, t>0
//
//  u(x,0) = g0(x), u(xa,t) = g1(t), u(xb,t) = g2(t) 
//
//  Implicit scheme
//
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <unistd.h>
#include <math.h>
#include "mpi.h"
#include "mycom.h"
#include "mynet.h"
#include "myio.h"
#include "myprog.h"

int np, mp, nl, ier, lp;//np-кол.процессов;mp-номер вызывающего
int mp_l, mp_r;
char pname[MPI_MAX_PROCESSOR_NAME];
char vname[10] = "ex13b_my";
char checkname[18] = "check_myparametr";//для проверки
char checknameABC[20] = "check_myparametrACBF";//для проверки 
char checknamey1[20] = "check_myparametry1";//для проверки 

char sname[20];
MPI_Status status;
union_t buf;
double tick, t1, t2, t3;

FILE *Fi = NULL;
FILE *Fo = NULL;
FILE *Fc = NULL;
FILE *Fabc = NULL;
FILE *Fy1 = NULL;

int nx, ntp, ntm, ntv;
double xa, xb, xk, x0, r0, q0, u0, u1;
double k1, k2, tau0, tau1, tmax, epst,epst1,epst2;
double tv, u10, omg0, omg1, gt;

//СВОИ ФУНКЦИИ 
double Krwff(double Swfs_up_loc);
double Krwff(double Swfs_up_loc) {
 //return 0.5;//0.04E-15;
 return 0.03*pow(Swfs_up_loc,2) +0.002*Swfs_up_loc+0.002;
}

double Kroff(double Swfs_up_loc);
double Kroff(double Swfs_up_loc) {
 // return 0.5;//0.6E-15;
 return 7.7*pow(Swfs_up_loc,4)-12.07*pow(Swfs_up_loc,3)+6.9*pow(Swfs_up_loc,2)-1.8*Swfs_up_loc+0.2;
}
//Плотности 
double DofMfunction(double P);//  плотность нефти
double DofMfunction(double P) {
  //return 730.0;?????????????????????
  return 710.0+(5e-6)*(P-101325);// нужно подправить альфа 5e-6
}

double DwfMfunction(double P);//плотность воды
double DwfMfunction(double P) {
 // return 1118.0;
 double alw;
 alw=(1118.0-1000.0)/(25000000.0-101325.0);
  return 1000.0+alw*(P-101325); // нужно подправить альфа 5e-6
}
double DofM_Pfunction(double P);//  (плотность нефти)' по давлению
double DofM_Pfunction(double P) {
  return 5e-6*P;
 
}
double DwfM_Pfunction(double P);//(плотность воды)' по давлению 
double DwfM_Pfunction(double P) {
  return 5e-6*P;
 }
//Dwfwdfff
double Dwfwdfunction(double P);//  плотность  воды с весами 
double Dwfwdfunction(double P) {

  return 1118.0;
  //DwfMfunction(double P1)*d+(1-d)*DwfMfunction(P2) 
 
}
double Dofwdfunction(double P);//  плотность   нефти с весами
double Dofwdfunction(double P) {
  return 730.0;
  //DofMfunction(double P1)*d+(1-d)*DofMfunction(P2)
 
}

//вес d1fff
double d(double P1,double P2);//  вес от пористоcти(давление)  
double d(double P1,double P2) 
{
  return 0.5;
// return sqrt(m(P1))/(sqrt(m(P1))+sqrt(m(P2)));
}
//вес пористость
double m(double P1);//  пористость(давление)
double m(double P1) 
{
  return 0.01;
  
}

int main(int argc, char *argv[])
{
  int i, j, ii, i1, i2, nc, ncm, ncp, ncx,S_count,z;
  double hx, hx2, tau, gam, s0, s1, s2, s3,time;
  double *xx, *y0, *y1, *y2, *y3, *y4, *al;

  double *kfhff,*Dwfwdfff,*DwfMff,*MwfMff,*Dofwdfff,*DofMff,*MofMff,*Swfsff,*Swfw_Dwfw_dfff,*mfw_Dwfw_Pf,*Swf1w_Dofw_dfff,*mfw_Dofw_Pf,*Prfsff;
  double *Apkff,*Cpkff,*Bpkff,*Fpkff;
  double *mfw_Dwfw,*mfw_Dwfw_t0,*mfw_Dofw,*mfw_Dofw_t0,*Prfsff_t1,*Prfsff_time;
  double *hclff,*hknff;

  double Swfs_up_loc,ves;
  MyNetInit(&argc,&argv,&np,&mp,&nl,pname,&tick);

  fprintf(stderr,"Netsize: %d, process: %d, system: %s, tick=%12le\n",np,mp,pname,tick);
  sleep(1);

  sprintf(sname,"%s.p%02d",vname,mp);
  ier = fopen_m(&Fo,sname,"wt");
  
  //sprintf(checkname,".dat");

  fopen_m(&Fc,checkname,"wt");
  fopen_m(&Fy1,checknamey1,"wt");
  fopen_m(&Fabc,checknameABC,"wt");
  if (ier!=0) mpierr("Protocol file not opened",1);

  if (mp==0) {
    sprintf(sname,"%s.d",vname);
    ier = fopen_m(&Fi,sname,"rt");
    if (ier!=0) mpierr("Data file not opened",2);
    i = fscanf(Fi,"xa=%le\n",&xa);
    i = fscanf(Fi,"xb=%le\n",&xb);
    i = fscanf(Fi,"xk=%le\n",&xk);
    i = fscanf(Fi,"x0=%le\n",&x0);
    i = fscanf(Fi,"r0=%le\n",&r0);
    i = fscanf(Fi,"q0=%le\n",&q0);
    i = fscanf(Fi,"u0=%le\n",&u0);
    i = fscanf(Fi,"u1=%le\n",&u1);
    i = fscanf(Fi,"k1=%le\n",&k1);
    i = fscanf(Fi,"k2=%le\n",&k2);
    i = fscanf(Fi,"tau0=%le\n",&tau0);
    i = fscanf(Fi,"tau1=%le\n",&tau1);
    i = fscanf(Fi,"tmax=%le\n",&tmax);
    i = fscanf(Fi,"epst=%le\n",&epst);
    i = fscanf(Fi,"nx=%d\n",&nx);
    i = fscanf(Fi,"ntp=%d\n",&ntp);
    i = fscanf(Fi,"ntm=%d\n",&ntm);
    i = fscanf(Fi,"lp=%d\n",&lp);
    fclose_m(&Fi);
    if (argc>1) sscanf(argv[1],"%d",&nx);
    if (argc>2) sscanf(argv[2],"%d",&ntp);
    if (argc>3) sscanf(argv[3],"%d",&ntm);
  }

  if (np>1) {
    if (mp==0) {
      buf.ddata[0]  = xa;
      buf.ddata[1]  = xb;
      buf.ddata[2]  = xk;
      buf.ddata[3]  = x0;
      buf.ddata[4]  = r0;
      buf.ddata[5]  = q0;
      buf.ddata[6]  = u0;
      buf.ddata[7]  = u1;
      buf.ddata[8]  = k1;
      buf.ddata[9]  = k2;
      buf.ddata[10] = tau0;
      buf.ddata[11] = tau1;
      buf.ddata[12] = tmax;
      buf.ddata[13] = epst;
      buf.idata[28] = nx;
      buf.idata[29] = ntp;
      buf.idata[30] = ntm;
      buf.idata[31] = lp;
    }
    MPI_Bcast(buf.ddata,16,MPI_DOUBLE,0,MPI_COMM_WORLD);
    if (mp>0) {
      xa   = buf.ddata[0];
      xb   = buf.ddata[1];
      xk   = buf.ddata[2];
      x0   = buf.ddata[3];
      r0   = buf.ddata[4];
      q0   = buf.ddata[5];
      u0   = buf.ddata[6];
      u1   = buf.ddata[7];
      k1   = buf.ddata[8];
      k2   = buf.ddata[9];
      tau0 = buf.ddata[10];
      tau1 = buf.ddata[11];
      tmax = buf.ddata[12];
      epst = buf.ddata[13];
      nx   = buf.idata[28];
      ntp  = buf.idata[29];
      ntm  = buf.idata[30];
      lp   = buf.idata[31];
    }
  }

  fprintf(Fo,"Netsize: %d, process: %d, system: %s, tick=%12le\n",
    np,mp,pname,tick);

  fprintf(Fo,"xa=%le xb=%le xk=%le x0=%le r0=%le\n",xa,xb,xk,x0,r0);
  fprintf(Fo,"q0=%le u0=%le u1=%le k1=%le k2=%le\n",q0,u0,u1,k1,k2);
  fprintf(Fo,"tau0=%le tau1=%le tmax=%le epst=%le\n",tau0,tau1,tmax,epst);
  fprintf(Fo,"nx=%d ntp=%d ntm=%d lp=%d\n",nx,ntp,ntm,lp);

  t1 = MPI_Wtime();

//  u10 = u1 - u0; omg0 = 1.0 / tau0; omg1 = 1.0 / tau1;
  hx = (xb-xa)/nx; hx2 = hx * hx;
  //tau = 0.5 * hx / sqrt(dmax(k1,k2));
  //tau = dmin(tau,1.0/q0); gam = tau / hx2;
  //s0 = dmin(tmax/tau,1000000000.0); ntm = imin(ntm,(int)s0);

  fprintf(Fo,"u10=%le omg0=%le omg1=%le\n",u10,omg0,omg1);
  fprintf(Fo,"hx=%le tau=%le ntm=%d\n",hx,tau,ntm);
  if (mp == 0) fprintf(stderr,"nx=%d hx=%le tau=%le ntm=%d\n",nx,hx,tau,ntm);

  if (mp ==    0) mp_l = -1; else mp_l = mp - 1;
  if (mp == np-1) mp_r = -1; else mp_r = mp + 1;

  MyRange(np,mp,0,nx,&i1,&i2,&nc);
  ncm = nc-1; ncp = 2*(np-1); ncx = imax(nc,ncp);

  fprintf(Fo,"i1=%d i2=%d nc=%d\n",i1,i2,nc);

  xx = (double*)(malloc(sizeof(double)*nc));
  y0 = (double*)(malloc(sizeof(double)*nc));
  y1 = (double*)(malloc(sizeof(double)*nc));

  

  Apkff = (double*)(malloc(sizeof(double)*nc));
  Bpkff = (double*)(malloc(sizeof(double)*nc));
  Cpkff = (double*)(malloc(sizeof(double)*nc));
  Fpkff=(double*)(malloc(sizeof(double)*(nc)));
  
  

  
  al = (double*)(malloc(sizeof(double)*ncx));
//свои переменные
  kfhff=(double*)(malloc(sizeof(double)*(nc+1)));
  Dwfwdfff=(double*)(malloc(sizeof(double)*(nc+1)));
  DwfMff=(double*)(malloc(sizeof(double)*(nc+1)));
  MwfMff=(double*)(malloc(sizeof(double)*(nc+1)));
  Dofwdfff=(double*)(malloc(sizeof(double)*(nc+1)));
  DofMff=(double*)(malloc(sizeof(double)*(nc+1)));
  MofMff=(double*)(malloc(sizeof(double)*(nc+1)));
  Swfsff=(double*)(malloc(sizeof(double)*(nc+1)));
  Swfw_Dwfw_dfff=(double*)(malloc(sizeof(double)*(nc+1)));
  mfw_Dwfw_Pf=(double*)(malloc(sizeof(double)*(nc+1)));
  Swf1w_Dofw_dfff=(double*)(malloc(sizeof(double)*(nc+1)));
  mfw_Dofw_Pf=(double*)(malloc(sizeof(double)*(nc+1)));

  Prfsff=(double*)(malloc(sizeof(double)*(nc+1)));
  Prfsff_t1=(double*)(malloc(sizeof(double)*(nc+1)));
  Fpkff=(double*)(malloc(sizeof(double)*(nc+1)));
  mfw_Dwfw=(double*)(malloc(sizeof(double)*(nc+1)));
  mfw_Dwfw_t0=(double*)(malloc(sizeof(double)*(nc+1)));
  mfw_Dofw=(double*)(malloc(sizeof(double)*(nc+1)));
  mfw_Dofw_t0=(double*)(malloc(sizeof(double)*(nc+1)));

  hclff=(double*)(malloc(sizeof(double)*(nc+1)));
  hknff=(double*)(malloc(sizeof(double)*(nc+1)));


  if (np>1) {
    y2 = (double*)(malloc(sizeof(double)*nc));
    y3 = (double*)(malloc(sizeof(double)*nc));
    y4 = (double*)(malloc(sizeof(double)*9*ncp));
  }

fprintf(Fc,"%20s %20s %20s %20s %20s %20s\n","kfhff", "Dwfwdfff","DwfMff","MwfMff","Swfsff","Dofwdfff");

  ntv = 0; tv = 0.0; gt = 1.0;

  for (i=0; i<=nc; i++) {
  xx[i] = xa + hx * (i); 
  y1[i] = 0.0;
  Prfsff[i]=(25E6-(0.1*nc))+0.1;// явный 
  Prfsff[i]=25E6;// явный
  Prfsff[1]=22000000.0;// явный
  Prfsff_t1[i]=25E6;//неявный
  }// initial profile

//шаг
for (i=0; i<=nc+1; i++) {hclff[i]=0.0;}
for (i=1; i<=nc; i++) {hclff[i]=hx;}//*hx
for (i=0; i<=nc; i++) {hknff[i]=0.5*(hclff[i]+hclff[i+1]);}


// Time loop:
time=360;//1*24*60*60
S_count=0;
epst1=0.001;
epst2=0.001;
//time=time*24.0* 60.0*60.0;
tau=0.1;
z=time/tau;
Prfsff_time=(double*)(malloc(sizeof(double)*(z)));
printf("nc=%d\n",nc);
printf("z=%d\n",z);
  do {
    ///
 //  tau=0.1;
    //printf("tau=%20.8le\n",tau);
    ///
    S_count=0;
    z=0;

for (i=0; i<=nc; i++)
{Apkff[i]=0.0;Cpkff[i]=0.0;Bpkff[i]=0.0;Fpkff[i]=0.0;}

for (i=0; i<=nc+1; i++) {
 ves=d(Prfsff[i],Prfsff_t1[i]);
      
   
                
        kfhff[i]=1.0E-12/hknff[i];//проницаемость/шаг 
        Dwfwdfff[i]=1.0/(ves*DwfMfunction(Prfsff[i])+(1.0-ves)*DwfMfunction(Prfsff_t1[i]));//обратная плотность воды c весами
        DwfMff[i]=DwfMfunction(Prfsff[i]);//плотность воды(к-+0.5)
        MwfMff[i]=0.67E-3;// вязкость воды
        Dofwdfff[i]=1.0/(ves*DofMfunction(Prfsff[i])+(1.0-ves)*DofMfunction(Prfsff_t1[i]));//обратная плотность нефти с весами
        DofMff[i]=DofMfunction(Prfsff[i]);//плотность нефти(к-+0.5)
        MofMff[i]=0.86E-3;//вязкость нефти
        Swfsff[i]=i*(1.0/nc);//(0,1)
        Swfw_Dwfw_dfff[i]=0.36*Dwfwdfff[i]; //насыщенность воды*обратная плотность воды
        
        mfw_Dwfw_Pf[i]=hknff[i]*0.01*DwfM_Pfunction(Prfsff[i]);//шаг*пористость* (плотность воды)' 
        Swf1w_Dofw_dfff[i]=0.64*Dofwdfff[i];//насыщенность нефти*обратная плотность нефти
        mfw_Dofw_Pf[i]=hknff[i]*0.01*DofM_Pfunction(Prfsff[i]);//шаг*пористость*(плотность нефти)'
        
        mfw_Dwfw[i]=0.01*DwfMfunction(Prfsff_t1[i]);//1.1e-4//_пористость*плотность воды (неявный слой)
        mfw_Dwfw_t0[i]=0.01*DwfMfunction(Prfsff[i]);//_пористость*плотность воды (явный слой)
        mfw_Dofw[i]=0.01*DofMfunction(Prfsff_t1[i]);//_пористость*плотность нефти (неявный слой)
        mfw_Dofw_t0[i]=0.01*DofMfunction(Prfsff[i]);//_пористость*плотность нефти (явный слой)

        //fprintf(Fc,"%20.13le %20.13le %20.13le %20.13le %20.13le %20.13le\n",kfhff[i],Dwfwdfff[i],DwfMff[i],MwfMff[i],Swfsff[i],Dofwdfff[i]);
  }


//Apkff->
for (i=0; i<=nc; i++) 
  {
    if (i==0) continue;
    if (Prfsff[i]>Prfsff[i-1]) {Swfs_up_loc=Swfsff[i];}
	    else //!(kNicl_loc+1!up)
      {Swfs_up_loc=Swfsff[i-1];}

	 Apkff[i] = tau*kfhff[i]*\
   (Dwfwdfff[i]*(DwfMff[i]/MwfMff[i])*Krwff(Swfs_up_loc)+\
    Dofwdfff[i]*(DofMff[i]/MofMff[i])*Kroff(Swfs_up_loc));
  //fprintf(Fc,"%20.13le %20.13le %20.13le %20.13le %20.13le %20.13le\n",tau,Dwfwdfff[i],DwfMff[i],MwfMff[i],Swfsff[i],Dofwdfff[i]);
  //printf("%20.25le\n",Cpkff[i]);
  }

  //Bpkff->
for (i=0; i<=nc; i++) 
  {
    if (i==nc) continue;
    if (Prfsff[i]>Prfsff[i+1]) {Swfs_up_loc=Swfsff[i];}
	    else //!(kNicl_loc+1!up)
      {Swfs_up_loc=Swfsff[i+1];}

	 Bpkff[i] = tau*kfhff[i+1]*\
   (Dwfwdfff[i]*(DwfMff[i+1]/MwfMff[i+1])*Krwff(Swfs_up_loc)+\
    Dofwdfff[i]*(DofMff[i+1]/MofMff[i+1])*Kroff(Swfs_up_loc));
    
  }
  
   //Cpkff ->
    for (i=0; i<=nc; i++) { 
	

	 if (i!=0) {//then !no_left_bound

	 if (Prfsff[i]>Prfsff[i-1]) //then !(kNicl_loc!up)
	 Swfs_up_loc=Swfsff[i];
	 else //(kNicl_loc-1!up)
	 Swfs_up_loc=Swfsff[i-1];
       

	 Cpkff[i]=Cpkff[i]+\
      tau*kfhff[i]*(\
      Dwfwdfff[i]*(DwfMff[i]/MwfMff[i])*\
      Krwff(Swfs_up_loc)\
      +\
      Dofwdfff[i]*(DofMff[i]/MofMff[i])*\
      Kroff(Swfs_up_loc));
   }


       if (i!=nc)
       { //then !no_rigtt_bound
       if (Prfsff[i]>Prfsff[i+1]) //then !(kNicl_loc!up)
	 Swfs_up_loc=Swfsff[i];
	 else //!(kNicl_loc+1!up)
	 Swfs_up_loc=Swfsff[i+1];
     
	 Cpkff[i]=Cpkff[i]+\
      tau*kfhff[i+1]*(\
      Dwfwdfff[i]*(DwfMff[i+1]/MwfMff[i+1])*\
     Krwff(Swfs_up_loc)\
      +\
     Dofwdfff[i]*(DofMff[i+1]/MofMff[i+1])*\
     Kroff(Swfs_up_loc)\
      );\
       }
	 
       Cpkff[i]=Cpkff[i]+\
      Swfw_Dwfw_dfff[i]*mfw_Dwfw_Pf[i]\
      +\
      Swf1w_Dofw_dfff[i]*mfw_Dofw_Pf[i];  
    	 }
//



    for (i=0; i<=nc; i++) {y0[i] = y1[i];}

    for (i=0; i<=nc; i++) {
    ii = i1 + i;

	 if (i!=0){

	 if (Prfsff[i]>Prfsff[i-1]){
	 Swfs_up_loc=Swfsff[i];}
	 else 
	 {Swfs_up_loc=Swfsff[i-1];}
       

	 Fpkff[i]=Fpkff[i]+\
   tau*kfhff[i]*( \

   Dwfwdfff[i]*(DwfMff[i]/MwfMff[i])*\
   Krwff(Swfs_up_loc)\
    +\
   Dofwdfff[i]*(DofMff[i]/MofMff[i])*\
   Kroff(Swfs_up_loc))\
   *(Prfsff[i]-Prfsff[i-1]); //!(-!darsy)*(-!left flux)=+  нужно исправить
    //printf("%20.25le\n",Fpkff[i]);
    } //!no_left_bound

  if (i!=nc) // !no_rigtt_bound
       {if (Prfsff[i]>Prfsff[i+1]) {
	 Swfs_up_loc=Swfsff[i];}
	 else 
	 {Swfs_up_loc=Swfsff[i+1];}
       //!(kNicl_loc!up.or.down)

	 Fpkff[i]=Fpkff[i]+\
      tau*kfhff[i+1]*(\
      Dwfwdfff[i]*(DwfMff[i+1]/MwfMff[i+1])*\
      Krwff(Swfs_up_loc)\
      +\
      Dofwdfff[i]*(DofMff[i+1]/MofMff[i+1])*\
      Kroff(Swfs_up_loc)\
      )*(-(Prfsff[i+1]-Prfsff[i]));// !(-!darsy)=-
      //fprintf(Fc,"%20.20le %20.13le %20.13le %20.13le %20.13le %20.13le \n",Fpkff[i],tau*kfhff[i+1],Dwfwdfff[i],DwfMff[i+1]/MwfMff[i+1],Krwff(Swfs_up_loc),Dofwdfff[i]*(DofMff[i+1]/MofMff[i+1])*\
      Kroff(Swfs_up_loc));
     }// !no_rigtt_bound

       Fpkff[i]=Fpkff[i]+\
      Swfw_Dwfw_dfff[i]*\
      (mfw_Dwfw[i]-mfw_Dwfw_t0[i]) \
      +\
      Swf1w_Dofw_dfff[i]*\
      (mfw_Dofw[i]-mfw_Dofw_t0[i]);
       //
       Fpkff[i]=-Fpkff[i];// нужно спросить ????????????
    //fprintf(Fabc,"%20.25le %20.25le %20.25le %20.25le\n",Apkff[i],Cpkff[i],Bpkff[i],Fpkff[i]);
    }
   
  //printf("%20.25le\n",Cpkff[0]);
  ier = prog_right(nc,Apkff,Bpkff,Cpkff,Fpkff,al,y1);
   

    if (ier!=0) {printf("%d\n",ier);mpierr("Bad solution",1);};


for (i=0; i<=nc; i++)
{
if(i%100==0)
{

    //fprintf(Fc,"%20.13le %20.13le %20.13le \n",Fpkff[i],Swfw_Dwfw_dfff[i],mfw_Dwfw[i]-mfw_Dwfw_t0[i]);
  //  fprintf(Fabc,"%20.25le %20.25le %20.25le %20.25le\n",Apkff[i],Cpkff[i],Bpkff[i],Fpkff[i]);
    // нужно спросить 
    //fprintf(Fy1,"%20.25le %20.25le %20.25le\n",y1[i],epst1*Prfsff_t1[i]+epst2,(-Apkff[i]*y1[i-1]+Cpkff[i]*y1[i]-Bpkff[i]*y1[i+1])-Fpkff[i]);
}


//Prfsff[i]=Prfsff_t1[i];
//Prfsff_t1[i]=Prfsff_t1[i]+y1[i];
Prfsff_t1[i]=Prfsff[i];
Prfsff[i]=Prfsff[i]+y1[i];

}
Prfsff_time[ntv]=Prfsff_t1[0];
//printf("%20.8le %20.8le\n",Prfsff[0],y1[0]);
if(ntv%100==0){
fprintf(Fc,"%20.13le %20.13le\n",tv,Prfsff_time[ntv]);
}
//printf("%20.8le\n",Prfsff_time[ntv]);
ntv++; tv += tau;


/*
 for (i=0; i<=nc; i++){
    if (y1[i]<(epst1*Prfsff_t1[i]+epst2))
    {
      Prfsff[i]=Prfsff_t1[i];
      Prfsff_t1[i]=Prfsff_t1[i]+y1[i];
      ntv++; tv += tau;
    }
    else{Prfsff[i]=Prfsff_t1[i];
      Prfsff_t1[i]=Prfsff_t1[i]+y1[i];
      ntv++; tv += tau;
    S_count=S_count+1;
    if (S_count>6)
    {
      Prfsff[i]=25E6;
      Prfsff_t1[i]=25E6;
      tau=tau/2.0;
      tv=0;
      ntv=0;
    }
    }

    }
    */
    
    //(y1<epst1*Prfsff_t1[i]+epst2)&&(S_count<6);

    
  }while ((tv<time)); 
  //while ((tv<time)||(ntv<1000));
  

  t1 = MPI_Wtime() - t1;

  sprintf(sname,"%s_%02d.dat",vname,np);
  OutFun1DP(sname,np,mp,nc,xx,Prfsff_t1);

  fprintf(Fo,"ntv=%d tv=%le gt=%le time=%le\n",ntv,tv,gt,t1);
  if (mp == 0) fprintf(stderr,"ntv=%d tv=%le gt=%le tcpu=%le\n",ntv,tv,gt,t1);

  ier = fclose_m(&Fo);
  fclose_m(&Fc);
  fclose_m(&Fabc);
  fclose_m(&Fy1);
  MPI_Finalize();
  return 0;
}
