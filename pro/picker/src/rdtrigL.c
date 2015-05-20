#include <errno.h>
#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <stdarg.h>
#include <time.h>
#include <sys/types.h>
#include <sys/stat.h>
#include <dirent.h>
#include <libmseed.h>
#define VERSION "1.1"
#define PACKAGE "rtpick"
#include "ew_bridge.h"
#include "PickData.h"
#include "FilterPicker5_Memory.h"
#include "FilterPicker5.h"
//-------------------------------------------------
// the following includes will produce some
// graphics that can be used for a web browser
//-------------------------------------------------
//MKS making it so that it will always do graphics to return to normal 
//remove one line below --i.e., remove defint RTQUAKE 1 line
//#define RTQUAKE "1"
#ifdef RTQUAKE
#include "gdfontl.h"
#include "gdfontg.h"
#include "gdfonts.h"
#include "gd.h"
#endif

#define MAXCHA 300
#define MAXPAR 10
#define MAXSMP 200000
#define PACKAGE "rtpick"
//---------------------------------------------
// mseed2askii
//---------------------------------------------
struct listnode {
  char *key;
  char *data;
  struct listnode *next;
};

static void appl(int fnd,int nsamp);
static void design();
static int  buroots();
static void DTE();
static void TIMSEC();
//---------------------------------------------
// mseed2askii
//---------------------------------------------
static int writeascii (MSTrace *mst,int comp_cnt);
static struct listnode *addnode (struct listnode **listroot, void *key, int keylen,
				 void *data, int datalen);
//------------------------------------------------------
static int parameter_proc (int argcount, char **argvec);
static char *getoptval (int argcount, char **argvec, int argopt);
static void usage (void);
/**********************************************************************/
//---------------------------------------------
// mseed2askii
//---------------------------------------------
static int    verbose      = 0;      /* Verbosity level */
static int    reclen       = -1;     /* Record length, -1 = autodetected */
static int    indifile     = 1;      /* Individual file processing flag */
static char  *outputfile   = 0;      /* Output file name for single file output */
static FILE  *ofp          = 0;      /* Output file pointer */
static int    outformat    = 1;      /* Output file format */
static double timetol      = -1.0;   /* Time tolerance for continuous traces */
static double sampratetol  = -1.0;   /* Sample rate tolerance for continuous traces */


/* A list of input files */
struct listnode *filelist = 0;
//---------------------------------------------
int    mail1          =     0;       // 0-no mail, 1-mail
int    mail2          =     0;       // 0-no mail, 1-mail
int    mail3          =     0;       // 0-no mail, 1-mail1
int    mail4          =     0;       // 0-no mail, 1-mail
int    mail5          =     0;       // 0-no mail, 1-mail
int    keep           =     0;       // parameter to save start of existing s-file, 0: make new
int    prt            =     0;       // the higher value, the more printouts
char   sfilename[256];               // s-filename, given as input, complete path
char   wfilename[256];               // waveform filename
char   dbname[256];                  // database name, extracted from s-filename path
int    iterations     =   200;       // number of iterations removing phases giving bad residuals, input
int    locate         =     0;       // 1 = autolocate, 0 = no location, input
int    automag        =     0;       // 1 = automag, 0 = no automag
float  maxres         =   3.0;       // stop iterations when this value is reached, input
int    wavefiledef    =     0;       // switch to indicate if wavefilename is given as input
int    no_stations_trg=     5;       // min number of stations with phase reading needed to locate
int    sendmail=0;
char   mailaddress1[256];            // mailaddress 
char   mailaddress2[256];            // mailaddress 
char   mailaddress3[256];            // mailaddress 
char   mailaddress4[256];            // mailaddress 
char   mailaddress5[256];            // mailaddress 
char   mail_message[1024];           // mail content
char   TIME_START[MAXCHA][80];       // start time each component
char   HYPFILE[MAXCHA][500];         // temporary storage for lines from hyp.out
int    nusmp[MAXCHA];                // number of samples each component
float  vflt[MAXCHA][MAXSMP];         // filtered data components
float  vraw[MAXCHA][MAXSMP];         // raw data components  
char   STATION_NAMES[MAXCHA][80];    // station name: XX_OSL___00_EHZ, where XX is missing network code
char   trig_stations[MAXCHA][256];   // station name code (KONO) is extracted from STATION_NAMES
char   trig_components[MAXCHA][256]; // station component code (HHZ) is extracted from STATION_NAMES
int    srates[500];                  // samplerate each component, samples per second
char   COMP_HDR[MAXCHA][256];        // component header including sr, ns, start time etc, read from miniseed file
int    verdier[MAXCHA][MAXSMP];      // samples all components read from miniseed file
char   TRIGGER_TIMES_P[MAXCHA][256]; // p-phase pick from filterpicker
char   TRIGGER_TIMES_S[MAXCHA][256]; // s-phase pick from filterpicker
char   PICKLINES[MAXCHA][1024];      // result lines from filterpicker
char   PICKTIMES[MAXCHA][500];       // contains number of picks + picktimes extracted from PICKLINES

/**********************************************************************/  

//****************  FILTER  **********************************
/* Type-independant absolute value for C standard types */
#ifndef abs
#define ABS(x) ((x)<0?-(x):(x))
#endif
#define TRUE 1
static  int nsects;
static  float sn[5*3+1],sd[5*3+1];

float   snc[MAXCHA][16];                     // filter coefficients
float   sdc[MAXCHA][16];                     // filter coefficients
int     cset[MAXCHA];
int     cset_pick[MAXCHA];

static  float save_y2[MAXCHA][16];           // overlap between 2 filter buffers
static  float save_y1[MAXCHA][16];
static  float save_x2[MAXCHA][16];
static  float save_x1[MAXCHA][16];
static  float save_y2_pick[MAXCHA][16];      // overlap between 2 filter buffers
static  float save_y1_pick[MAXCHA][16];
static  float save_x2_pick[MAXCHA][16];
static  float save_x1_pick[MAXCHA][16];

float   fl;
float   fh;
float   eqavst;
int     orden;



typedef struct {float r, i;} complex;

/* Table of constant values */
static  float c_b12 = 2.;
static  complex c_b43 = {1.,0.};


/* Local function prototypes */
static  double warp();
static  int chebparm(), lptbp(), lptbr(), lpthp(), bilin2(), lp();
static  int c1roots(), c2roots();
static  int cutoffs(), beroots();
static  complex cmul(), cpowi(), csqrt1(), conjg(), cdiv();

//****************  END FILTER  **********************************


float picks[2][MAXSMP];
float ratios[MAXSMP];
float env[MAXSMP];
float xtap[MAXSMP];
float ltas[MAXSMP];
float ssttaa[MAXSMP];
float h[MAXSMP];
float smratios[MAXSMP];
float pk[2][MAXSMP];  

  float maxsmA;
  float maxsmB;
  int pa;
  int   varighet[MAXCHA];
  int thr[MAXCHA];
  int tmn[MAXCHA];
  int tsc[MAXCHA];
  int tms[MAXCHA];
  int start_trg[MAXCHA];


int printarg=0;
double trg_tim[MAXCHA];             // MSECS for klon n and channel m
int klon;






  int no_lines;
float fseismo[MAXSMP];
char *topdir        = 0;               // top directory SEISAN
char *topdir_rt     = 0;               // top directory RTQUAKE





int PeakA(int npts,float dt,int cmp)

{
  float mx=0.0f;
  int top=0;
  int pp;
  for(pp=0;pp<npts;pp++)                                   // find max peak in filtered ratios
  {
    if(picks[1][pp]>mx)
    {
      mx=picks[1][pp];
      pa=pp;
    }
  }
  maxsmA=mx;
  top=(int)(picks[0][pa]/dt);
  return top;
}
int PeakB(int pix,int npts,float dt,int cmp)
{
  int b1,b2;
  b1 = pa - 20;
  b2 = pa + 20;
  int p;
  float mx=0.0f;
  int top = 0;

  for(p=0;p<npts;p++)                                      // find max peak in filtered ratios
  {
    if(p<b1 || p>b2)
    {
      if (picks[1][p] > mx)
      {
        mx = picks[1][p];
        top = p;
      }
    }
  }
  maxsmB=mx;
  top=(int)(picks[0][top]/dt);
  return top;
}
int MPICKS(float tstart,float dt,int npts,float cutoff)
{
  int npk;
  int idx_min;
  int idx_inflex;
  float curr_slope;
  float prev_slope;
  float curr_inflex;
  float prev_inflex;
  int ii,i;
  float loc_cut;
cutoff=2.0;
  loc_cut = cutoff;


for(;;)
{
  npk = 0;
  idx_min = 1;
  idx_inflex = 0;

  for(ii=0;ii<npts;ii++)
  {
    if(smratios[ii] != 0.0)
      break;
  }
  prev_slope = smratios[ii+1] - smratios[ii];
  prev_inflex= (smratios[ii+2] + smratios[ii])/2 - smratios[ii+1];

  for(i=(ii+1);i<(npts-ii+1);i++)
  {
    curr_slope = smratios[i+1] - smratios[i];
    curr_inflex = (smratios[i+2] + smratios[i])/2 - smratios[i+1];

    if(prev_inflex >= 0.0 && curr_inflex <= 0.0)
      idx_inflex = i+1;

    if(prev_slope >= 0.0 && curr_slope < 0.0 && smratios[i+1] >= loc_cut && idx_inflex != 0)
    {
      pk[0][npk] = tstart + dt * idx_inflex;
      pk[1][npk] = smratios[i];
      npk = npk + 1;
    }
    prev_slope = curr_slope;
    prev_inflex = curr_inflex;
  }
  if(npk<=3)
  {
    break;
  }else{
    loc_cut=loc_cut + 0.1;
//  printf("LOC_CUT: %8.3f  NPK: %3d\n",loc_cut,npk);
  }
}
return(npk);
}
HANFILT_N(int n,int nfil)
{
  float filt[400];
  float pihalf;
  float tfilt;
  float s1;
  int k,ii,i;
  for(i=0;i<90000;i++)
    smratios[i]=0.0;

  pihalf = 3.1415927/2.0;
  tfilt = 0.0;

  for(k=0;k<(2*nfil);k++)
  {
    s1=(float)(k)/(float)(nfil+1);
    filt[k] = cos(s1*pihalf) * cos(s1*pihalf);
    tfilt = tfilt + filt[k]; 
  }
  for(ii=0;ii<n;ii++)
  {
    if(ratios[ii] != 0.0)
      break;
  }
  for(i=(ii+nfil);i<(n-nfil);i++)
  {
    smratios[i] = 0.0;
    for(k=0;k<(2*nfil);k++)
    {
      smratios[i] = smratios[i] + ratios[i+k-nfil] * filt[k];
    }
    smratios[i] = smratios[i]/tfilt;
  }
  for(i=0;i<(nfil+ii);i++)
    smratios[i]=0.0;
  for(i=0;i<nfil;i++)
    smratios[n+1-i]=0.0;
}
STALTA_N(int npts,int n1,int n2)
{
  float lta;
  float sta;
  float lta_tot;
  float sta_tot;
  int i,ii,nz,i1,j;
  for(i=0;i<(n1+n2);i++)
    ratios[i] = 0.0;
  for(ii=0;ii<npts;ii++)
  {
    if(env[ii] != 0.0)
    {
      nz = ii;
      break;
    }
  }
  i1=nz+n1-1;
  lta_tot = 0.0;
  for(j=nz;j<i1;j++)
    lta_tot = lta_tot + (fabs)(env[j]);
  lta = lta_tot/n1;
  
  sta_tot = 0.0;
  for(j=(i1+1);j<(i1+n2);j++)
    sta_tot = sta_tot + (fabs)(env[j]);
  sta = sta_tot/n2;

  ratios[i1+n2] = sta/lta;

  for(i=(i1+1);i<(npts-2);i++)
  {
    lta_tot = lta_tot + (fabs)(env[i]) - (fabs)(env[i-n1]);
    sta_tot = sta_tot + (fabs)(env[i+n2]) - (fabs)(env[i]);
    lta = lta_tot/n1;
    sta = sta_tot/n2;

    ratios[i+n2] = sta/lta;
  }
  for(i=(npts-50);i<npts;i++)
    ratios[i]=0.0;
}
TAPER1_N(int n,float frac)
{
  float xlen;
  float taplen;
  float angle;
  float fact;
  float rest;
  int   hel;
  float fhel;
  float pihalf=3.1415927/2.0;
  int i,i2;
  for(i=0;i<n;i++)
    xtap[i]=h[i];
  if(frac == 0.0)
    return(0);;
  xlen = (float)(n-1);
  taplen = frac * xlen/2.0;
  hel= (int)(taplen);
  fhel = hel;
  rest= taplen-hel;
  if(rest > 0.5)
  {
    i2=hel+1;
  }else{
    i2=hel;
  }
  for(i=0;i<i2;i++)
  {
    angle = pihalf * (float)(i-1)/(i2-1);
    fact = sin(angle) * sin(angle);
    xtap[i] = h[i] * fact;
    xtap[n+1-i] = h[n+1-i] * fact;
  }
}
TILBERT_N(float dt,int npts,int nfil)
{
  float b[MAXSMP];
  int i;
  int j;
  int i1;
  int i2;
  int ii;
  int nfiltot;
  float t;
  float dt0;
  float npts0;
  float pi=3.141592654;

  if(nfil == 0)
  {
    for(i=0;i<npts;i++)
      env[i] = fabs(fseismo[i]);
    return(0);
  }
  j=0;
  for(i = -nfil;i<nfil;i++)
  {
    t=(float)(i) * dt;
    if(i != 0)
    {
      h[j] = -1.0/(pi*t);
    }else{
      h[j] = 0.0;
    }
    j++;
  }
  nfiltot = 2*nfil+1;
  TAPER1_N(nfiltot,1.0);
  dt0 = dt;
  npts0 = npts;
  for(i=0;i<npts;i++)
  {
    b[i]=0.0;
    env[i]=0.0;
  }
  i1=1+nfil;
  i2=npts - nfil;
  for(i=i1;i<i2;i++)
  {
    for(j=0;j<(2*nfil);j++)
    {
      ii=i+nfil-j;
      b[i]=b[i]+fseismo[ii] * xtap[j];
    }
    b[i]=b[i]*dt;
  }
  for(i=i1;i<i2;i++)
    env[i]= sqrt(fseismo[i]*fseismo[i] + b[i]*b[i]);

  for(i=i1;i<i2;i++)
    env[i]= env[i] * env[i];
}

//----------------------------------------------------------------------
// MAKEPICKS makes estimates of phase arrival times
//
// INPUTS:
//
// a      : STA/LTA vector
// tstart : time of first element in a (seconds)
// dt     : time increment (seconds)
// npts   : number of points in a
// cutoff : minimum amplitude of local max of a to be considered a pick
//
// OUTPUT :
//
// p      : array holding pick information
// npk    : number of picks
//--------------------------------------------------------------------------
int MAKEPICKS_N(int npts, float dt,float tstart, float sta,float lta,float cutoff,int nenvfil,float hanfil,int npks)
{

  int ista,ilta,nhanfil,minpts,nonzero,ii;
  int maxpts = 150000;
  ista = (int)(sta/dt);
  ilta = (int)(lta/dt);
  nhanfil = (int)(hanfil/dt);
  minpts = 2*(nenvfil + ista + ilta + nhanfil);
  nonzero = 0;
  for(ii=0;ii<npts;ii++)
  {
    if(fseismo[ii] != 0.0)
      nonzero++;
  }
  if(nonzero < minpts)
  {
    printf("Not non-zero points enough to process\n");
    npks = 0;
    return(-1);
  }
  TILBERT_N(dt,npts,nenvfil);
  STALTA_N(npts,ilta,ista);
  HANFILT_N(npts,nhanfil);
  npks=MPICKS(tstart,dt,npts,cutoff);
  return(npks);
}

//----------------------------------------------------------------------
// MAKEPICKS makes estimates of phase arrival times
//
// INPUTS:
//
// a      : STA/LTA vector
// tstart : time of first element in a (seconds)
// dt     : time increment (seconds)
// npts   : number of points in a
// cutoff : minimum amplitude of local max of a to be considered a pick
//
// OUTPUT :
//
// p      : array holding pick information
// npk    : number of picks
//--------------------------------------------------------------------------
int MAKEPICKS(float a[], float tstart, float dt, int npts,float cutoff,int pindex,int cmp)
{
  int idx_min=0;
  int idx_inflex=0;
  int ii=0;
  int i;
  int logi=0;
  int npkidx=-1;
  int npk;
  float curr_slope=(float)0.0;
  float prev_slope=(float)0.0;
  float curr_inflex=(float)0.0;
  float prev_inflex=(float)0.0;


  idx_min=1;
  idx_inflex=0;
  a[0]=a[1];

  for(ii=pindex;ii<npts;ii++)
  {
    if(a[ii] != (float)0.0)
      break;
  }
//  printf("PICK First non-zero at: %d cutoff: %6.3f\n",ii,cutoff);
  prev_slope=a[ii+1] - a[ii];
  prev_inflex=(a[ii+2]+a[ii])/(float)2.0 - a[ii+1];
//  printf("PR_SL: %6.3f PR_IN: %6.3f\n",prev_slope,prev_inflex);

  for( i=(ii+1);i<(npts);i++)
  {
    curr_slope=a[i+1]-a[i];
    curr_inflex=(a[i+2]+a[i])/(float)2.0 - a[i+1];

    if(prev_inflex >= (float)0.0 && curr_inflex <= (float)0.0)
    {
//  printf("%4d PR_SL: %6.3f PR_IN: %6.3f CU_SL: %6.3f CU_IN: %6.3f\n",i,prev_slope,prev_inflex,curr_slope,curr_inflex);
      idx_inflex=i+1;
//      if(idx_inflex<1)
//        printf("IDX_INFLEX: %6d\n",idx_inflex);
    }
// look for local max
    if(prev_slope >= (float)0.0 && curr_slope < (float)0.0 && a[i+1] >= cutoff && idx_inflex != 0)
    {
      npkidx++;
      picks[0][npkidx]=tstart+dt*(float)idx_inflex;
      picks[1][npkidx]=a[i];
//      printf("NPK : %4d  %6.3f  %6.3f\n",npkidx,picks[0][npkidx],picks[1][npkidx]);
    }
    prev_slope=curr_slope;
    prev_inflex=curr_inflex;

  }
  npk=npkidx+1;
//  printf("Npicks: %d\n",npk);
  npkidx=-1;
  return(npk);
}


void STALTAPHFL(int ne,int ws,int wl)
  {
    float extra=(float)0.0;
    float sta1=(float)0.0;
    float lta1=(float)0.0;
    float k1=(float)0.0;
    float k2=(float)0.0;
    float R=(float)0.0;
    int s,ss;
    sta1=(float)0.0;
    lta1=(float)0.0;
    
//      System.out.println("STALTAPHFL: ne:"+ne);

    for(s=1;s<800;s++)
    {
      extra=(float)fseismo[s]-(float)fseismo[s-1];
      extra=extra * extra;
      k1=fabsf((float)fseismo[s-1]);
      k2=k1 * k1;
      sta1=sta1+(k2+ (float)3.0*extra-sta1)/(float)ws;
      k1=fabsf((float)fseismo[s]);
      k2=k1 * k1;
      lta1=lta1+(k2-lta1)/(float)wl;
    }

    for(ss=1;ss<(ne-1);ss++)
    {
      extra=(float)fseismo[ss]-(float)fseismo[ss-1];
      extra=extra * extra;
      k1=fabsf((float)fseismo[ss-1]);
      k2=k1 * k1;
      sta1=sta1+(k2+ (float)3.0*extra-sta1)/(float)ws;
      ssttaa[ss] = sta1;
      k1=fabsf((float)fseismo[ss]);
      k2=k1 * k1;
      lta1=lta1+(k2-lta1)/(float)wl;
      R=sta1/lta1;
      ratios[ss]=R;
      ltas[ss] = lta1;
    }

  }
void TRIGGER_TIME_S(char timestamp_trg[],int ws,int ss,int kanal,float dt)
{
  int DDAY,MMTH;
  int DY;
  int YR;
  long   yrn;
  long   mon;
  long   day;
  long   hrn;
  long   min;
  long   doy;

  int   iyrn;
  int   imon;
  int   iday;
  int   ihrn;
  int   imin;

  float  sek;
  double MSECS;

  int seconds_cpu;
  char tiden_cpu[80];
  char millis_cpu[4];
  long day_cpu;
  long mon_cpu;
  long yrn_cpu;
  long hrn_cpu;
  long min_cpu;
  float sek_cpu;
  char dm[10];
  double MSECS_CPU;

  char   timestamp[80];
  char dum[50];
  int i,j;
  int klon=0;
  sscanf(timestamp_trg,"%4d%1c%2d%1c%2d%1c%2d%1c%2d%1c%f6.3",&YR,dum,&imon,dum,&iday,dum,&ihrn,dum,  &imin,dum,&sek);
  hrn=ihrn;
  min=imin;
  mon=imon;
  day=iday;
  yrn=YR;
  yrn=yrn-2000;

  iday=day;
  imon=mon;
  iyrn=yrn;
  ihrn=hrn;
  imin=min;

  TIMSEC(yrn,mon,day,hrn,min,sek,&MSECS);
  MSECS=MSECS+(dt * ss)-(dt * ws);
  SECTIM(MSECS,&yrn,&doy,&mon,&day,&hrn,&min,&sek);

  iday=day;
  imon=mon;
  iyrn=yrn;
  ihrn=hrn;
  imin=min;
  iyrn=iyrn-2000;

  sprintf(timestamp,"%2d/%2d/%2d %2d:%2d:%6.3f",iday,imon,iyrn,ihrn,imin,sek);
  if(timestamp[0] == ' ') timestamp[0] = '0';
  if(timestamp[3] == ' ') timestamp[3] = '0';
  if(timestamp[6] == ' ') timestamp[6] = '0';
  if(timestamp[9] == ' ') timestamp[9] = '0';
  if(timestamp[12] == ' ') timestamp[12] = '0';
  if(timestamp[15] == ' ') timestamp[15] = '0';
  if(timestamp[18] == ' ') timestamp[18] = '0';
  if(timestamp[19] == ' ') timestamp[19] = '0';
  if(timestamp[20] == ' ') timestamp[20] = '0';

  sprintf(TRIGGER_TIMES_S[kanal],"%s",timestamp);

  trg_tim[kanal]=MSECS;                      // save trigger time

}
void TRIGGER_TIME_P(char timestamp_trg[],int ws,int pindex,int kanal,float dt)
{
  int DDAY,MMTH;
  int DY;
  int YR;
  long   yrn;
  long   mon;
  long   day;
  long   hrn;
  long   min;
  long   doy;

  int   iyrn;
  int   imon;
  int   iday;
  int   ihrn;
  int   imin;

  float  sek;
  double MSECS;

  int seconds_cpu;
  char tiden_cpu[80];
  char millis_cpu[4];
  long day_cpu;
  long mon_cpu;
  long yrn_cpu;
  long hrn_cpu;
  long min_cpu;
  float sek_cpu;
  char dm[10];
  double MSECS_CPU;

  char   timestamp[80];
  char dum[50];
  int i,j;
  int klon=0;

  sscanf(timestamp_trg,"%4d%1c%2d%1c%2d%1c%2d%1c%2d%1c%f6.3",&YR,dum,&imon,dum,&iday,dum,&ihrn,dum,  &imin,dum,&sek);
  hrn=ihrn;
  min=imin;
  day=iday;
  mon=imon;
  yrn=YR;
  yrn=yrn-2000;

  iday=day;
  imon=mon;
  iyrn=yrn;
  ihrn=hrn;
  imin=min;

  TIMSEC(yrn,mon,day,hrn,min,sek,&MSECS);
  MSECS=MSECS+(dt * pindex)-(dt * ws);
  SECTIM(MSECS,&yrn,&doy,&mon,&day,&hrn,&min,&sek);

  iday=day;
  imon=mon;
  iyrn=yrn;
  ihrn=hrn;
  imin=min;
  iyrn=iyrn-2000;

  sprintf(timestamp,"%2d/%2d/%2d %2d:%2d:%6.3f",iday,imon,iyrn,ihrn,imin,sek);
  if(timestamp[0] == ' ') timestamp[0] = '0';
  if(timestamp[3] == ' ') timestamp[3] = '0';
  if(timestamp[6] == ' ') timestamp[6] = '0';
  if(timestamp[9] == ' ') timestamp[9] = '0';
  if(timestamp[12] == ' ') timestamp[12] = '0';
  if(timestamp[15] == ' ') timestamp[15] = '0';
  if(timestamp[18] == ' ') timestamp[18] = '0';
  if(timestamp[19] == ' ') timestamp[19] = '0';
  if(timestamp[20] == ' ') timestamp[20] = '0';

  sprintf(TRIGGER_TIMES_P[kanal],"%s",timestamp);

  trg_tim[kanal]=MSECS;                      // save trigger time

}
int FIND_P(int tot)
{
  int indx=0;
  int tel=0;
  int k;
  float max;
  max=0.0;
  for(k=0;k<tot;k++)
  {
    if(ratios[k]>max)
    {
      max=ratios[k];
      indx=tel;
    }
    tel++;
  }
  return(indx-5);
}



void Scale_to_1(int no_of_smp,float MxAmp)                 // Scale channel to max 1
{
  int k;

  for(k=0;k<no_of_smp;k++)
  {
    fseismo[k]=fseismo[k]/MxAmp;
  }

}
void Scale_to_100_Raw(int no_of_smp,int cmpkl,float MxAmp)                 // Scale channel to max 1
{
  int k;

  for(k=0;k<no_of_smp;k++)
  {
    vraw[cmpkl][k]=(vraw[cmpkl][k]/MxAmp)*85.0;
  }

}
void Scale_to_100_Flt(int no_of_smp,int cmpkl,float MxAmp)                 // Scale channel to max 1
{
  int k;

  for(k=0;k<no_of_smp;k++)
  {
    vflt[cmpkl][k]=(vflt[cmpkl][k]/MxAmp)*85.0;
  }

}
void MaxAmp(int no_of_smp,float *MxAmp)                    // Find max amplitude channel
{
  float max;
  int k;

  max=(float)0.0;
  for(k=0;k<no_of_smp;k++)
  {
    if(fabsf(fseismo[k])>max)
      max=fabsf(fseismo[k]);
  }
  *MxAmp=max;
}
void MaxAmpFlt(int no_of_smp,int cmpkl,float *MxAmp)                    // Find max amplitude channel
{
  float max;
  int k;

  max=(float)0.0;
  for(k=0;k<no_of_smp;k++)
  {
    if(fabsf(vflt[cmpkl][k])>max)
      max=fabsf(vflt[cmpkl][k]);
  }
  *MxAmp=max;
}
void MaxAmpRaw(int no_of_smp,int cmpkl,float *MxAmp)                    // Find max amplitude channel
{
  float max;
  int k;

  max=(float)0.0;
  for(k=0;k<no_of_smp;k++)
  {
    if(fabsf(vraw[cmpkl][k])>max)
      max=fabsf(vraw[cmpkl][k]);
  }
  *MxAmp=max;
}
void ZERO_START()
{
  int s;
  for (s=0;s<300;s++)
    fseismo[s]=fseismo[s+300];
}
main(int argc,char *argv[])
{
//---------------------------------------------
// mseed2askii
//---------------------------------------------
  MSTraceGroup *mstg = 0;
  MSTrace *mst;
  MSRecord *msr = 0;
  
  struct listnode *flp;
  
  int retcode;
  int totalsamps = 0;
//--------------------------------------

  FILE *sfile;
  int    comp_cnt=0;
  char record[200];
  char yr_trg[200];
  char mo_trg[200];
  char da_trg[200];
  char hr_trg[200];
  char mi_trg[200];
  char se_trg[200];
  int  br;
  char cmd01[200];
  char cmd02[200];
  char cmd03[200];

  char *triggerfile    = 0;               // station list for recording
  char *detfile=0;
  char streamsfull[200]; 
  char filnam[200];
  char fnam[200];

  char MINTRGTID[80];
  char sfilepath[300];


  char det[80];
  char clone[10];

  
  char comp[10];
  int i,k,nn,mm,p,l;
  int tel;
  int nchannels;
  int ret;
  int klon;

  char from[200];
  char dm[200];
  char temptid[200];
  char temptid1[200];
  char temptid2[200];
  char temptid3[200];
  char temptid4[200];
  int  readings;

//-------------------------------------------------------------------------
// Process command line parameters
//-------------------------------------------------------------------------
  if (parameter_proc (argc, argv) < 0)
  {
    fprintf(stderr, "Parameter processing failed\n\n");
    fprintf(stderr, "Try '-h' for detailed help\n");
    return -1;
  }
//------------------------------------------------------------------------
// set path to SEISAN_TOP as default
//------------------------------------------------------------------------
  topdir = (char*)getenv("SEISAN_TOP");
  if(topdir)
  {
    printf("RTPICK: SEISAN_TOP.......................: %s\n", topdir);    
  }else{
    printf("RTPICK: SEISAN_TOP not defined.\n");
    printf("Install SEISAN !\n");
    exit(0);    
  }
//------------------------------------------------------------------------
// set path to RTQUAKE_TOP as default
//------------------------------------------------------------------------
#ifdef RTQUAKE
  topdir_rt = (char*)getenv("RTQUAKE_TOP");
  if (topdir_rt)
  {
    printf("RTPICK: RTQUAKE_TOP......................: %s\n", topdir_rt);
  }else{
    printf("RTPICK: RTQUAKE_TOP not defined\n");
    printf("Include setup_rt.bash or setup_rt.csh in your .bashrc or .csh\n");
    exit(0);      
  } 
#endif
//------------------------------------------------------------------------
//  if exist, otherwise use defaults
//------------------------------------------------------------------------
Read_Parameters();

//-------------------------------------------------------------------
// sfilename given as input. Complete path and filename
//--------------------------------------------------------------------
  sprintf(sfilepath,"%s",sfilename);
  if(prt >= 1)  
    printf("rtpick.........................: %s\n",sfilepath);
  sprintf(cmd03,"cp %s s_orginal.sav",sfilepath);                  // save original s-file to current directory
  if(prt >= 1)
    printf("Save original s-file to s_original.sav to current directory\n");
  ret=system(cmd03);
//-------------------------------------------------------------------
// extract database name from s-filepath
//-------------------------------------------------------------------
  for(i=0;i<256;i++)
  {
    if(sfilepath[i] == 'R' && sfilepath[i+1] == 'E' & sfilepath[i+2] == 'A')
    {
      dbname[0] = sfilepath[i+4];
      dbname[1] = sfilepath[i+5];
      dbname[2] = sfilepath[i+6];
      dbname[3] = sfilepath[i+7];
      dbname[4] = sfilepath[i+8];
      dbname[5] = '\0';
    }
  }
  if(prt >= 1)
    printf("rtpick: dbname.................: %s\n",dbname);
//-------------------------------------------------------------------
// open s-file and read the '1' line to get the time for trigger
// continue to read line '6' to get the wavefile name
//-------------------------------------------------------------------
  if ((sfile = fopen (sfilepath, "rb")) == NULL)
  {
    printf("Can't open sfile file: %s\n",sfilepath);
    exit(0);
  }
  for(i=0;i<50;i++)
  {
    br = les_linje(sfile,record);  // read first line in new station, and save samplerate + timestamp
    if(prt >= 4)
      printf("%s\n",record);
    if(record[79]=='1')
    {
      sscanf(record,"%s",yr_trg);
      dm[0]=record[6];
      if(dm[0]==' ') dm[0]='0';
      dm[1]=record[7];
      dm[2]='\0';
      sscanf(dm,"%s",mo_trg);
      dm[0]=record[8];

      if(dm[0]==' ') dm[0]='0';
      dm[1]=record[9];
      dm[2]='\0';

      sscanf(dm,"%s",da_trg);
      dm[0]=record[11];
      if(dm[0]==' ') dm[0]='0';
      dm[1]=record[12];
      dm[2]='\0';
      sscanf(dm,"%s",hr_trg);
      dm[0]=record[13];
      if(dm[0]==' ') dm[0]='0';
      dm[1]=record[14];
      dm[2]='\0';
      sscanf(dm,"%s",mi_trg);
      dm[0]=record[16];
      if(dm[0]==' ') dm[0]='0';
      dm[1]=record[17];
      dm[2]=record[18];
      dm[3]=record[19];
      dm[4]='\0';
      sscanf(dm,"%s",se_trg);
//      printf("%s %s %s %s %s %s\n",yr_trg,mo_trg,da_trg,hr_trg,mi_trg,se_trg);
      sprintf(MINTRGTID,"%s/%s/%s0%s:%s:%s",da_trg,mo_trg,yr_trg,hr_trg,mi_trg,se_trg);
      printf("RTPICK:..................................: MINTRGTID: %s\n",MINTRGTID);
    }
    if(record[79]=='6')
      break;
  }
  clone[0]='0';
  clone[1]='\0';
  klon=0;
  sprintf(fnam,"%s",record);
  for(i=0;i<100;i++)
  {
    fnam[i]=fnam[i+1];
    if(fnam[i]==' ')
    {
      fnam[i]='\0';
      break;
    }
  }
  if(prt >= 1)
  {
    printf("Extract waveform filename from original s-file\n");
    printf("rtpick.........................: %s\n",fnam);
  }
  fclose(sfile);

//--------------------------------------------------------------------------------------
// Temporary solution:
// Wavefile is in miniseed format and is converted to askii for reading 
// This should be a proper miniseed read routine
//--------------------------------------------------------------------------------------
  if(wavefiledef == 1)
  {
      /* Add the file name to the intput file list */
    if ( ! addnode (&filelist, NULL, 0, wfilename, strlen(wfilename)+1) )
    {
      fprintf (stderr, "Error adding file name to list\n");
    }
  }else{
  sprintf(cmd02,"%s/WAV/%s/%s/%s/%s",topdir,dbname,yr_trg,mo_trg,fnam);
  printf("CMD02: %s\n",cmd02);
      /* Add the file name to the intput file list */
    if ( ! addnode (&filelist, NULL, 0, cmd02, strlen(cmd02)+1) )
    {
      fprintf (stderr, "Error adding file name to list\n");
    }
  }


  if(prt >= 4)
    printf(" CH      STATION                 TIMESTAMP   SR     DT  NSAMP PINDEX     PP    PN    SS    SN  CODA   DIFF-P   DIFF-S\n");
//---------------------------------------------------------------------------------------
// Read wavefile that was converted from miniseed to askii: trigger.ask
// Each component is read and processed for phases
//---------------------------------------------------------------------------------------
  comp_cnt=0;
  /* Init MSTraceGroup */
  mstg = mst_initgroup (mstg);
  
  /* Read input miniSEED files into MSTraceGroup */
  flp = filelist;

      if ( verbose )
        fprintf (stderr, "Reading %s\n", flp->data);
      
      while ( (retcode = ms_readmsr(&msr, flp->data, reclen, NULL, NULL,
				    1, 1, verbose-1)) == MS_NOERROR )
	{
	  if ( verbose > 1)
	    msr_print (msr, verbose - 2);
	  
	  mst_addmsrtogroup (mstg, msr, 1, timetol, sampratetol);
	  

	  totalsamps += msr->samplecnt;
	}
      
      if ( retcode != MS_ENDOFFILE )
	fprintf (stderr, "Error reading %s: %s\n", flp->data, ms_errorstr(retcode));
      
      /* Make sure everything is cleaned up */
      ms_readmsr (&msr, NULL, 0, NULL, NULL, 0, 0, 0);
      
      /* If processing each file individually, write ASCII and reset */
      if ( indifile )
	{
	  mst = mstg->traces;
	  while ( mst )
	    {
	      writeascii (mst,comp_cnt);
              comp_cnt++;
	      mst = mst->next;
	    }
	   mstg = mst_initgroup (mstg);
	}
      
  
  /* Make sure everything is cleaned up */
  mst_freegroup (&mstg);
  
  if ( ofp )
    fclose (ofp);
  if(prt >= 1)
    printf("Components found reading miniseed file. COMP_CNT: %d\n",comp_cnt);

  nchannels=read_trigger(comp_cnt);

  for(i=0;i<nchannels;i++)
  {
    
    nn=0;
    tel=0;
//printf("STATION_NAMES: %s\n",STATION_NAMES[i]);
    for(k=3;k<10;k++)
    {
      if(STATION_NAMES[i][k] == '_')
      {
        break;
      }else{
        tel++;
      }
    }
    nn=tel;                      // station name length
//printf("nn: %d\n",nn);    
    for(k=0;k<10;k++)
      comp[k]='\0';
    for(k=0;k<nn;k++)
      comp[k]= STATION_NAMES[i][k+3];
    sprintf(trig_stations[i],"%s",comp);
//printf("trig_stations: %s\n",trig_stations[i]);
    for(k=0;k<10;k++)
      comp[k]='\0';
    switch(nn)
    {
      case 3:
      for(k=0;k<3;k++)
      	comp[k]= STATION_NAMES[i][k+12];
      sprintf(trig_components[i],"%s",comp);
      break;
      case 4:
      for(k=0;k<3;k++)
        comp[k]= STATION_NAMES[i][k+12];
      sprintf(trig_components[i],"%s",comp);
      break;
      case 5:
      for(k=0;k<3;k++)
        comp[k]= STATION_NAMES[i][k+12];
      sprintf(trig_components[i],"%s",comp);	
      break;
    }
//printf("trig_components: %s\n",trig_components[i]);
  } // i

for(k=0;k<nchannels;k++)
{
  for(i=0;i<200;i++)
  {
    temptid1[i]='\0';
    temptid2[i]='\0';
    temptid3[i]='\0';
    temptid4[i]='\0';
  }
  sscanf(PICKTIMES[k],"%d",&readings);

  switch(readings)
  {
    case 0:               // 0 phase readings
    if(prt >= 2)
      printf("Number of readings: %s Comp: %s  ",PICKTIMES[k],trig_components[k]);
    sprintf(TRIGGER_TIMES_P[k],"                         ");
    sprintf(TRIGGER_TIMES_S[k],"                         ");
    if(prt >= 2)
      printf("%3d %2d %s %s\n",readings,k,trig_stations[k],trig_components[k]);
    break;

    case 1:               // 1 phase reading
    if(prt >= 2)    
      printf("Number of readings: %2d %s Comp: %s\n",readings,PICKTIMES[k],trig_components[k]);
    if(trig_components[k][2] == 'Z')
    {
      sscanf(PICKTIMES[k],"%4c%21c",dm,temptid);
      temptid[21]='\0';
      sprintf(TRIGGER_TIMES_P[k],"%s",temptid);
      sprintf(TRIGGER_TIMES_S[k],"                         ");
      if(prt >= 2)
        printf("%3d %s %s %s +\n",k,trig_stations[k],trig_components[k],temptid);
    }else{                
      if(prt >= 2)
        printf("%3d %s %s Not valid on horizontal\n",k,trig_stations[k],trig_components[k]);  

      sprintf(TRIGGER_TIMES_P[k],"                         ");
      sprintf(TRIGGER_TIMES_S[k],"                         ");
    }
    break;

    case 2:
    if(prt >= 2)
      printf("Number of readings: %2d %s Comp: %s\n",readings,PICKTIMES[k],trig_components[k]);
    if(trig_components[k][2] == 'Z')
    {
      sscanf(PICKTIMES[k],"%4c%21c",dm,temptid1);
      temptid1[21]='\0';
      sprintf(TRIGGER_TIMES_P[k],"%s",temptid1);
      sprintf(TRIGGER_TIMES_S[k],"                         ");
      if(prt >= 2)
        printf("%3d %s %s %s +\n",k,trig_stations[k],trig_components[k],temptid1);
    }else{
      sscanf(PICKTIMES[k],"%26c%21c",dm,temptid2);
      temptid2[21]='\0';
      sprintf(TRIGGER_TIMES_P[k],"                         ");
      sprintf(TRIGGER_TIMES_S[k],"%s",temptid2);
      if(prt >= 2)
        printf("%3d %s %s %s +\n",k,trig_stations[k],trig_components[k],temptid2);
    }
    break;

    case 3:
    if(prt >= 2)
      printf("Number of readings: %2d %s Comp: %s\n",readings,PICKTIMES[k],trig_components[k]);      
    if(trig_components[k][2] == 'Z')
    {
      sscanf(PICKTIMES[k],"%4c%21c",dm,temptid1);
      temptid1[21]='\0';
      sprintf(TRIGGER_TIMES_P[k],"%s",temptid1);
      sprintf(TRIGGER_TIMES_S[k],"                         ");
      if(prt >= 2)
        printf("%3d %s %s %s +\n",k,trig_stations[k],trig_components[k],temptid1);
    }else{
      sscanf(PICKTIMES[k],"%48c%21c",dm,temptid3);
      temptid3[21]='\0';
      sprintf(TRIGGER_TIMES_P[k],"                         ");
      sprintf(TRIGGER_TIMES_S[k],"%s",temptid3);
      if(prt >= 2)
        printf("%3d %s %s %s +\n",k,trig_stations[k],trig_components[k],temptid3);
    }
    break;
    case 4:
    if(prt >= 2)
      printf("Number of readings: %2d %s Comp: %s\n",readings,PICKTIMES[k],trig_components[k]);      
    if(trig_components[k][2] == 'Z')
    {
      sscanf(PICKTIMES[k],"%4c%21c",dm,temptid1);
      temptid1[21]='\0';
      sprintf(TRIGGER_TIMES_P[k],"%s",temptid1);
      sprintf(TRIGGER_TIMES_S[k],"                         ");
      if(prt >= 2)
        printf("%3d %s %s %s +\n",k,trig_stations[k],trig_components[k],temptid1);
    }else{
      sscanf(PICKTIMES[k],"%48c%21c",dm,temptid3);
      temptid3[21]='\0';
      sprintf(TRIGGER_TIMES_P[k],"                         ");
      sprintf(TRIGGER_TIMES_S[k],"%s",temptid3);
      if(prt >= 2)
        printf("%3d %s %s %s +\n",k,trig_stations[k],trig_components[k],temptid3);
    }
    break;
    case 5:
    if(prt >= 2)
      printf("Number of readings: %2d %s Comp: %s\n",readings,PICKTIMES[k],trig_components[k]);      
    if(trig_components[k][2] == 'Z')
    {
      sscanf(PICKTIMES[k],"%4c%21c",dm,temptid1);
      temptid1[21]='\0';
      sprintf(TRIGGER_TIMES_P[k],"%s",temptid1);
      sprintf(TRIGGER_TIMES_S[k],"                         ");
      if(prt >= 2)
        printf("%3d %s %s %s +\n",k,trig_stations[k],trig_components[k],temptid1);
    }else{
      sscanf(PICKTIMES[k],"%48c%21c",dm,temptid3);
      temptid3[21]='\0';
      sprintf(TRIGGER_TIMES_P[k],"                         ");
      sprintf(TRIGGER_TIMES_S[k],"%s",temptid3);
      if(prt >= 2)
        printf("%3d %s %s %s +\n",k,trig_stations[k],trig_components[k],temptid3);
    }
    break;
  }

}

  Create_Sfile(fnam,sfilepath,dbname,nchannels,klon,locate,yr_trg,mo_trg,da_trg,hr_trg,mi_trg,se_trg,sfilename);
//  printf("mks just called Create_Sfile \n");
  

  
#ifdef RTQUAKE
/*------------------------------------------------------------------------*/
/*                Create png plots for web-page                           */
/*------------------------------------------------------------------------*/
  PNG_PLOTS(nchannels);
//  printf("mks just called PNG_PLOTS \n");
/*------------------------------------------------------------------------*/
/*                Send mail to specified users                            */
/*------------------------------------------------------------------------*/  
  /*  mks removing sendmail for testing.  Delete line below to allow it again */
  /*   Allowed--now just uncomment for testing */
//  sendmail = 0;
  printf("Sendmail..: %2d\n",sendmail);
  if(sendmail ==1)
  {
    if(mail1==1)
    {
      printf("Sent mail to recipient 1\n");
      sprintf(mail_message,"echo %cDETECTION WITHOUT LOCATION \n %s%c | mutt -s %cTRIGGER%c -a %chyp.out%c -a %c%s/loc/ALL.png%c -- %s &",0x22,sfilepath,0x22,0x22,0x22,0x22,0x22,0x22,\
		topdir_rt,0x22,mailaddress1);

       printf("%s\n",mail_message);
      ret=system(mail_message);
    }
    if(mail2==1)
    {
      printf("Sent mail to recipient 2\n");
      printf("%s\n",mail_message);    
      ret=system(mail_message);
    }
    if(mail3==1)
    {
      printf("Sent mail to recipient 3\n");
      printf("%s\n",mail_message);    
      ret=system(mail_message);
    }
    if(mail4==1)
    {
      printf("Sent mail to recipient 4\n");
      printf("%s\n",mail_message);    
      ret=system(mail_message);
    }
    if(mail5==1)
    {
      printf("Sent mail to recipient 5\n");
      printf("%s\n",mail_message);    
      ret=system(mail_message);
    }  
  }else if(sendmail ==2)
  {
    if(mail1==1)
    {
      printf("Sent mail to recipient 1\n");
/* sprintf(mail_message,"echo %cDETECTION WITH LOCATION%c | mutt -s %cLOCATED EVENT%c -a %chyp.out%c -a %c%s/loc/ALL.png%c -- %s &",0x22,0x22,0x22,0x22,0x22,0x22,0x22,topdir_rt,0x22,mailaddress1); */
/* sprintf(mail_message,"echo %chttp://maps.googleapis.com/maps/api/staticmap?center=%f,%f&zoom=7&size=900x1000&maptype=hybrid&markers=color:red_7Ccolor:red_7Clabel:Q_7C%f,%f&sensor=false%c | mutt  -a %c%s/loc/ALL.png%c -s %cTRIGGER%c -- %s &",0x22,lat,lon,lat,lon,0x22,0x22,topdir_rt,0x22,0x22,0x22,mailaddress1);        */
	printf("%s\n",mail_message);
      ret=system(mail_message);
    }
    if(mail2==1)
    {
      printf("Sent mail to recipient 2\n");
      printf("%s\n",mail_message);    
      ret=system(mail_message);
    }
    if(mail3==1)
    {
      printf("Sent mail to recipient 3\n");
      printf("%s\n",mail_message);    
      ret=system(mail_message);
    }
    if(mail4==1)
    {
      printf("Sent mail to recipient 4\n");
      printf("%s\n",mail_message);    
      ret=system(mail_message);
    }
    if(mail5==1)
    {
      printf("Sent mail to recipient 5\n");
      printf("%s\n",mail_message);    
      ret=system(mail_message);
    }  
  }
#endif  
}
#ifdef RTQUAKE
PNG_PLOTS(int nchannels)
{
  FILE   *out;
  int    i;
  int    ix;
  int    kar;
  int    kar2;
  int    displ;
  int    ret2;
  int    kan;
  int    kl;
  int    nemin;
  int    nl;
  int    koda;
  int    rate;
  int    yincr;  
  int    avstand;
  int    bg,sg1,sg2,sg3,sg4,tl,tt,tk,t5,dg;
  int    ix_old,ix_new,iy_old,iy_new;
  int    cmpix;
  int    cmpix2;
  int    ret;
  int    stasjon;
  int    ncomp[MAXCHA];
  int    xsiz,ysiz;


  int    YR;
  int    DY;
  int    DDAY;
  int    MMTH;
  int   iyrn;
  int   imon;
  int   iday;
  int   ihrn;
  int   imin;
  int   mark10;
  int   mark5;
  
  double PPH[MAXCHA];
  double SPH[MAXCHA];
  double COD[MAXCHA];  
  double MSECS;
  double MINI;
  double MSECS_START[MAXCHA];  
  
  char   last_station[200];
  char   fil[200];
  char   fil2[200];
  char   timestring[256];
  char   dum[256];
  char   NAME_ONLY[MAXCHA][100];
  char   COMP_NAMES[MAXCHA][100];
  char   ph[100];
  char   navnihyp[200];
  char   comp[100];
  char   sec10[100];
  char   sec5[100];
  char   pngfile[200];
  char  nam[256];
  
  long   yrn;
  long   mon;
  long   day;
  long   hrn;
  long   min;
  long   doy;

  float  sek;
  float  diff;
  float  idx;
  float  incr;
  float  MAXAMP_RAW;
  float  MAXAMP_FLT;

  
 // printf("mks in PNG_PLOTS \n");
  for(i=0;i<nchannels;i++)                                 // initialize p-phase and s-phase table
  {
    PPH[i]=-1.0;
    SPH[i]=-1.0;
  }
  gdImagePtr im;                                           // Declare the image 


  kar=0;
  yincr=86;                                                // Set distance increment in pixels between traces
  avstand=yincr;

  for(i=0;i<200;i++)                                        // initialize last station name
    last_station[i]='\0';  

//---------------------------------------
// FIND NUMBER OF COMPONENTS EACH STATION
//---------------------------------------
  stasjon=0;
  for(cmpix = 0; cmpix < nchannels; cmpix++)
  {
    for(i=0;i<20;i++)
      fil[i]='\0';
    for(i=0;i<5;i++)
    {
      if(STATION_NAMES[cmpix][i+3] != '_')
      {
        fil[i]=STATION_NAMES[cmpix][i+3];
        kar++;
      }
    }
    for(cmpix2 = 0; cmpix2 < nchannels; cmpix2++)
    {
      for(i=0;i<20;i++)
        fil2[i]='\0';
      for(i=0;i<5;i++)
      {
        if(STATION_NAMES[cmpix2][i+3] != '_')
        {
          fil2[i]=STATION_NAMES[cmpix2][i+3];
          kar2++;
        }
      }
//      printf("%s  %s\n",fil,fil2);
      ret=strncmp(fil,fil2,kar);                           // compare current station with previous station
      if(ret == 0)                                         // New station ?
        stasjon++;
    }
    ncomp[cmpix]=stasjon;
    stasjon=0;
    kar=0;
    kar2=0;
  }
 
//------------------------------------------
// convert all start times to seconds
//------------------------------------------
  for(cmpix = 0; cmpix < nchannels; cmpix++)
  {
    sprintf(timestring,"%s",TIME_START[cmpix]);

    sscanf(timestring,"%4d%1c%2d%1c%2d%1c%2d%1c%2d%1c%f6.3",&YR,dum,&imon,dum,&iday,dum,&ihrn,dum,  &imin,dum,&sek);
    hrn=ihrn;
    min=imin;
    mon=imon;
    day=iday;
    yrn=YR;
    yrn=yrn-2000;
    iday=day;
    imon=mon;
    iyrn=yrn;
    ihrn=hrn;
    imin=min;
    TIMSEC(yrn,mon,day,hrn,min,sek,&MSECS);
    MSECS_START[cmpix]=MSECS;
  }

//---------------------------------------
// component name all channels
//---------------------------------------
  for(cmpix = 0; cmpix < nchannels; cmpix++)
  {
    dum[0]=STATION_NAMES[cmpix][12];
    dum[1]=STATION_NAMES[cmpix][14];
    dum[2]='\0';
    sprintf(COMP_NAMES[cmpix],"%s",dum);
  }
  sprintf(dum,"hyp.out"); 

//--------------------------------------
// read s-file
//--------------------------------------
  if(prt >= 4)
    printf("CALLING CHECK_SS %s\n",dum);
  ret=Check_SS(dum);
  if(prt >= 4)
    printf("RET: %d\n",ret);

//---------------------------------------
// extract the station name only  
//---------------------------------------
  for(cmpix = 0; cmpix < nchannels; cmpix++)
  {
    for(i=0;i<20;i++)
      fil[i]='\0';
    for(i=0;i<5;i++)
    {
      if(STATION_NAMES[cmpix][i+3] != '_')
      {
        fil[i]=STATION_NAMES[cmpix][i+3];
      }
    }
    sprintf(NAME_ONLY[cmpix],"%s",fil);
//    printf("NAME_ONLY: %s\n",NAME_ONLY[cmpix]);
  } 
  
//-------------------------------------------------------------------
// insert all readings of phases from s-file into corresponding
// place for each station component. (for plotting of phases)
//-----------------------------------------------------------------
  for(nl = 0; nl < no_lines; nl++)
  {
    sprintf(timestring,"%s",HYPFILE[nl]);
//    printf("timestring: %s\n",timestring);
    for(i=0;i<20;i++)
      navnihyp[i]='\0';
    strncpy(navnihyp,timestring,6);
//    printf("navnihyp  : %s\n",navnihyp);

    for(i=0;i<5;i++)
      navnihyp[i]=navnihyp[i+1];
    navnihyp[5]='\0';
//    printf("navnihyp  : %s\n",navnihyp);
//    getchar();    
    for(i=0;i<5;i++)
    {
      if(navnihyp[i]==' ')
        navnihyp[i]='\0';
    }
    sscanf(timestring,"%6c%s",dum,comp);
    sscanf(timestring,"%8c%s",dum,ph);
    dum[0]=timestring[18];
    dum[1]=timestring[19];
    dum[2]='\0';
    sscanf(dum,"%d",&ihrn);
    dum[0]=timestring[20];
    dum[1]=timestring[21];
    dum[2]='\0';
    sscanf(dum,"%d",&imin);
    sscanf(timestring,"%23c%f5.2",dum,&sek);
    koda=0;
    if(timestring[32] != ' ')
      sscanf(timestring,"%28c%d",dum,&koda);

    hrn=ihrn;
    min=imin;
    TIMSEC(yrn,mon,day,hrn,min,sek,&MSECS);
      
    for(cmpix = 0; cmpix < nchannels; cmpix++)
    {

      ret = strcmp(navnihyp,NAME_ONLY[cmpix]);
      ret2= strcmp(comp,COMP_NAMES[cmpix]);
      
//      printf("ret: %d  %s    %s    ret2: %d  %s   %s  Z%s\n",ret,navnihyp,NAME_ONLY[cmpix],ret2,comp,COMP_NAMES[cmpix],timestring); 
//      getchar();     
      if(ret == 0 && ret2 == 0)       // found station and component ?
      {
//      printf("%s  %s\n",navnihyp,NAME_ONLY[cmpix]);

	if(ph[1]=='P')
	{
	  PPH[cmpix]=MSECS;
//	  printf("%s %s P\n",navnihyp,comp);
	}
	if(ph[1]=='S')
	{
	  SPH[cmpix]=MSECS;
//	  printf("%s %s S\n",navnihyp,comp);	  
	}
	if(koda!=0)
	{
	  COD[cmpix]=MSECS+koda;
//	  printf("%s %s %d\n",navnihyp,comp,koda);	  
	}
      }
    }
  }        

//------------------------------------  
// GENERATE PLOT EACH STATION
//------------------------------------
for(cmpix = 0; cmpix < nchannels; cmpix++)
{
  for(i=0;i<20;i++)
    fil[i]='\0';
  if(STATION_NAMES[cmpix][0] != '_')
  {
     for(i=0;i<5;i++)
     {
       if(STATION_NAMES[cmpix][i+3] != '_')
       {
         fil[i]=STATION_NAMES[cmpix][i+3];
         kar++;
       }
     }
  }else{
     for(i=0;i<5;i++)
     {
       if(STATION_NAMES[cmpix][i+1] != '_')
       {
         fil[i]=STATION_NAMES[cmpix][i+1];
         kar++;
       }
     }	  
  }

  ret=strncmp(fil,last_station,kar);                       // compare current station with previous station
  if(ret != 0)                                             // New station ?
  {
    if(prt >= 1)
      printf("New station single station %s\n",fil);
    sprintf(pngfile,"%s/loc/%s.png",topdir_rt,fil);   // YES, file name format: BER.png or PB10.png
    sprintf(last_station,"%s",fil);


    xsiz = 795;
    ysiz = 345;
    
    xsiz = 845;
    ysiz = 345; 
    
    yincr=86;                                                // Set distance increment in pixels between traces
    avstand=yincr;
    
    
    im  = gdImageCreate(xsiz, ysiz);                        // Allocate image: tot_wd pixels width by tot_hg pixels high

    bg  = gdImageColorAllocate(im, 255, 255, 255); // WHITE first color definition defines the background color
    tl  = gdImageColorAllocate(im, 200, 200, 200); // GRAY  color for timelines
    tk  = gdImageColorAllocate(im,   0,   0,   0); // BLACK color for text
    tt  = gdImageColorAllocate(im,   0,   0, 255); // BLUE color for title
    t5  = gdImageColorAllocate(im,   0,   0,   0); // BLACK color for 5 minute marks
    sg1 = gdImageColorAllocate(im, 255,   0,   0); // RED   color seismogram
    sg2 = gdImageColorAllocate(im,   0, 255,   0); // GREEN color seismogram
    sg3 = gdImageColorAllocate(im,   0,   0, 255); // BLUE color seismogram
    sg4 = gdImageColorAllocate(im,   0,   0,   0); // BLACK color seismogram     

    ix_old=idx;
    iy_old=avstand;
    gdImageLine(im, 0, 0, xsiz-1, 0, sg4);                     // draw frame  
    gdImageLine(im, xsiz-1, 0, xsiz-1, ysiz-1, sg4);                 // draw frame
    gdImageLine(im, xsiz-1, ysiz-1, 0, ysiz-1, sg4);                 // draw frame
    gdImageLine(im, 0, ysiz-1, 0, 0, sg4);                     // draw frame
//------------------------------------------------- 
// find minimum start time in seconds this station       
//-------------------------------------------------
    MINI=MSECS_START[cmpix];
    kan=cmpix;

    for(kl=0;kl<ncomp[cmpix];kl++)
    {      
      if(MSECS_START[cmpix+kl] <= MINI)
      {
	MINI=MSECS_START[cmpix+kl];
	kan=cmpix+kl;
	nemin=nusmp[cmpix+kl];
	rate=srates[cmpix+kl];
      }
    }

//--------------------------------------------------
// draw station name and start time
//---------------------------------------------------   
    gdImageString(im,gdFontGetLarge(),10,15,fil,sg3);        
    gdImageString(im,gdFontGetLarge(),80,15,TIME_START[kan],tk);      
    sprintf(timestring,"%s",TIME_START[kan]);
    sscanf(timestring,"%4d%1c%3d%1c%2d%1c%2d%1c%f6.3",&YR,dum,&DY,dum,&ihrn,dum,&imin,dum,&sek);
    for(kl=0;kl<7;kl++)
    {
      if((kl * 10) > sek)
      {
	diff = (kl * 10) - sek;
	break;
      }
    }
    ret=kl*10;
    if(ret==60)
      ret=0;

    sprintf(sec10,"%2d",ret);
    incr = (float)(xsiz-1)/((float)rate * ((float)nemin/rate));
  
    ix_old = 0;
    ix_new = 0;
    iy_old = 335;
    iy_new = 335;
    idx    = 0;
    mark10= diff * rate;
    for(i=0;i<nusmp[cmpix];i++)
    {
      idx = idx + incr;
      ix_new = idx;
	
      gdImageLine(im, ix_old, iy_old, ix_new, iy_new, tl); // draw timescale
      if(i==mark10)
      {
        gdImageLine(im, ix_new, 40, ix_new, 335, tl);      // draw 10 seconds marks		    
        if(sec10[0]==' ')
	  gdImageLine(im, ix_new, 320, ix_new, 335, sg4);
        gdImageString(im,gdFontGetLarge(),ix_new-8,313,sec10,tk); // draw in black minute marks
        mark10 = mark10 + 10.0 * rate;
        ret=ret+10;
	if(ret == 60)
          ret = 0;
        sprintf(sec10,"%2d",ret);
      }
      ix_old=ix_new;
    }
//--------------------------------------------------------
// draw station components
//--------------------------------------------------------
    avstand=yincr;
    iy_old=yincr;
    iy_new=yincr;
   
    for(kl=0;kl<ncomp[cmpix];kl++)
    {
//printf("%2d ncomp: %d\n",kl,ncomp[kl]); 
      ix_old=incr * srates[cmpix] * (MSECS_START[cmpix+kl]-MINI );
      idx = ix_old;
      ix=nusmp[cmpix];
      MaxAmpRaw(ix,cmpix+kl,&MAXAMP_RAW);
      Scale_to_100_Raw(ix,cmpix+kl,MAXAMP_RAW);
      MaxAmpFlt(ix,cmpix+kl,&MAXAMP_FLT);
      Scale_to_100_Flt(ix,cmpix+kl,MAXAMP_FLT);
//printf("MAX_RAW: %10.3f  MAX_FLT: %10.3f\n",MAXAMP_RAW,MAXAMP_FLT);
      for(i=0;i<ix;i++)                          // plot this buffer
      {
        ix_new=(int)idx;
//          iy_new=(int)vraw[cmpix+kl][i]*1.0+avstand;
        iy_new=(int)vflt[cmpix+kl][i]*1.01+avstand;
        gdImageLine(im, ix_old, iy_old, ix_new, iy_new, sg4);  // draw from old to new point
        ix_old=ix_new;
        iy_old=iy_new;
        idx=idx+incr;
      }
//-----------------------------------------------
// draw phases this component
//-----------------------------------------------
      if(PPH[cmpix+kl] > 0)
      {

	gdImageSetThickness(im,2);
        sprintf(dum,"P");
        displ=incr * srates[cmpix] * (PPH[cmpix+kl] - MINI);
        gdImageLine(im, displ, avstand+30, displ, avstand-30, sg1); 
        gdImageLine(im, displ, avstand, displ-30, avstand-30, sg2);
        gdImageString(im,gdFontGetLarge(),displ-10,avstand-23,dum,sg3);
	gdImageSetThickness(im,1);	
      }
      if(SPH[cmpix+kl] > 0)
      {	
	gdImageSetThickness(im,2);
	sprintf(dum,"S");
        displ=incr * srates[cmpix] * (SPH[cmpix+kl] - MINI);
        gdImageLine(im, displ, avstand+30, displ, avstand-30, sg1);
        gdImageLine(im, displ, avstand, displ-30, avstand-30, sg2);
        gdImageString(im,gdFontGetLarge(),displ-10,avstand-23,dum,sg3); 
	gdImageSetThickness(im,1);	
      }
      if(COD[cmpix+kl] > 0)
      {	
	gdImageSetThickness(im,2);
	sprintf(dum,"C");
        displ=incr * srates[cmpix] * (COD[cmpix+kl] - MINI);
        gdImageLine(im, displ, avstand+30, displ, avstand-30, sg1);
        gdImageLine(im, displ, avstand, displ-30, avstand-30, sg2);
        gdImageString(im,gdFontGetLarge(),displ-10,avstand-23,dum,sg3); 
	gdImageSetThickness(im,1);	
      }
      idx=0;                                               // update counter before next component
      ix_old=0;
      avstand = avstand + yincr;
      iy_old=avstand;
      iy_new=avstand;
    }
    avstand=yincr;
//--------------------------------------------------
// write to png file
//--------------------------------------------------
    out = fopen (pngfile, "wb");
    gdImagePng (im, out);
    fclose (out);
    gdImageDestroy(im);                   /* Destroy the image in memory. */
  } 
}

  if(prt >= 1)
    printf("PLOT ALL CHANNELS\n");  
//-----------------------------------------------
// find minimum start time all components
//-----------------------------------------------
  MINI=MSECS_START[0];
  kan=cmpix;
  for(kl=0;kl<nchannels;kl++)
  {      
    if(MSECS_START[kl] < MINI)
    {
      MINI=MSECS_START[kl];
      kan=kl;
      nemin=nusmp[kl];
      rate=srates[kl];
    }
  }

//----------------------------------
// plot start time
//----------------------------------
  gdImageString(im,gdFontGetLarge(),80,15,TIME_START[kan],tk);      

  sprintf(timestring,"%s",TIME_START[kan]);
  sscanf(timestring,"%4d%1c%3d%1c%2d%1c%2d%1c%f6.3",&YR,dum,&DY,dum,&ihrn,dum,&imin,dum,&sek);
  for(kl=0;kl<7;kl++)
  {
    if((kl * 10) > sek)
    {
      diff = (kl * 10) - sek;
      break;
    }
  }

 
ret=kl*10;
if(ret==60)
  ret=0;
sprintf(sec10,"%2d",ret);
incr = (float)845/((float)rate * ((float)nemin/rate));
  if(prt >= 1)
    printf("RATE: %d NEMIN: %d INCR: %10.5f\n",rate,nemin,incr);


im  = gdImageCreate(xsiz, 3000);                 // Allocate the image: tot_wd pixels width by tot_hg pixels high
bg  = gdImageColorAllocate(im, 255, 255, 255); // WHITE first color definition defines the background color
tl  = gdImageColorAllocate(im, 200, 200, 200); // GRAY  color for timelines
tk  = gdImageColorAllocate(im,   0,   0,   0); // BLACK color for text
tt  = gdImageColorAllocate(im,   0,   0, 255); // BLUE color for title
t5  = gdImageColorAllocate(im,   0,   0,   0); // BLACK color for 5 minute marks
sg1 = gdImageColorAllocate(im, 255,   0,   0); // RED   color seismogram
sg2 = gdImageColorAllocate(im,   0, 255,   0); // GREEN color seismogram
sg3 = gdImageColorAllocate(im,   0,   0, 255); // BLUE color seismogram
sg4 = gdImageColorAllocate(im,   0,   0,   0); // BLACK color seismogram  
//------------------------------------  
// GENERATE PLOT ALL COMPONENTS
//------------------------------------
yincr = 50;
avstand=yincr;
sprintf(pngfile,"%s/loc/ALL.png",topdir_rt);     // file to write png

for(cmpix = 0; cmpix < nchannels; cmpix++)
{

  dg=sg4;      

  ix_old=0;
  iy_old=avstand;
  ix_new=0;
  iy_new=avstand;
 
  gdImageLine(im, 0, 0, xsiz-1, 0, dg);                       // plot frame
  gdImageLine(im, xsiz-1, 0, xsiz-1, 2950, dg);                   // plot frame
  gdImageLine(im, xsiz-1, 2950, 0, 2950, dg);                   // plot frame
  gdImageLine(im, 0, 2950, 0, 0, dg);                       // plot frame

  iy_old=avstand;
  iy_new=avstand;

  ix_old=incr * srates[cmpix] * (MSECS_START[cmpix]-MINI );

  idx = ix_old;
  /*  MKS  changing ix to be 60 seconds instead of given by sample number */
  //ix=nusmp[cmpix];
  ix=60*srates[cmpix];
  sprintf(dum,"%s",STATION_NAMES[cmpix]);
  strncpy(nam,dum,12);

  gdImageString(im,gdFontGetLarge(),10,iy_new-15,dum,sg3);
//--------------------------------------------------
// draw start time
//---------------------------------------------------   
//    gdImageString(im,gdFontGetLarge(),10,15,fil,sg3);        
    gdImageString(im,gdFontGetLarge(),80,15,TIME_START[kan],tk);      

   //---------------------
  // Draw 5 second markers  (MKS 29 Oct 2014)
  //-------------------
  sprintf(timestring,"%s",TIME_START[kan]);
  sscanf(timestring,"%4d%1c%3d%1c%2d%1c%2d%1c%f6.3",&YR,dum,&DY,dum,&ihrn,dum,&imin,dum,&sek);
//  printf("sek is %f \n",sek);
  for(kl=0;kl<7;kl++)
  {
    if((kl * 10) > sek)
    {
     	diff = (kl * 10) - sek;
	break;
    }
  }
  ret=kl*10;
  if(ret==60)
    ret=0;

  sprintf(sec10,"%2d",ret);
  rate=srates[cmpix];
  //Make it go 45 seconds past trigger, which is at 35 seconds so total is 45+35=80
  incr = (float)794/((float)rate * 80.);
//    incr = (float)(xsiz-1)/((float)rate * ((float)nemin/rate));  Old way
//  printf("cmpix: %d sta: %s  rate: %d ix: %d nusmp: %d incr: %f are \n",cmpix,STATION_NAMES[cmpix],rate,ix, nusmp[cmpix],incr);
  ix_old = 0;
  ix_new = 0;
  iy_old = 335; 
//  iy_old=avstand;
  iy_new = 335;
 // iy_new=avstand; 
//  idx    = incr*20/rate;  Trying to start earlier but gave up.
  idx = 0;
  mark5= diff * rate / 2;
  for(i=0;i<nusmp[cmpix];i++)
  {
    idx = idx + incr;
    ix_new = idx;
    gdImageLine(im, ix_old, iy_old, ix_new, iy_new, tl); // draw timescale
    if(i==mark5)
    {
//      gdImageLine(im, ix_new, 60, ix_new, 335, tl);      // draw 5 seconds marks
      gdImageLine(im, ix_new, avstand, ix_new, avstand+yincr, tl);      // draw 5 seconds marks
//      printf("sec10 is %s\n",sec10);
//      if(sec10[0]==' ')
      if((sec10[0]==' ' || sec10[0] == '3') && sec10[1]=='0')
      gdImageLine(im, ix_new, avstand+yincr/2, ix_new, avstand+yincr, sg4); //draw in black 30 sec mark
//      gdImageLine(im, ix_new, 320, ix_new, 335, sg4);
//      gdImageString(im,gdFontGetLarge(),ix_new-8,313,sec10,tk); // draw in black minute marks
      if(sec10[1]=='0')
      gdImageString(im,gdFontGetLarge(),ix_new-8,avstand,sec10,tk); // draw in labels every 10 sec
      mark5 = mark5 + 5.0 * rate;
      ret=ret+5;
      if(ret == 60)
        ret = 0;
      sprintf(sec10,"%2d",ret);
    }
    ix_old=ix_new;
    iy_new=iy_old+yincr;
    iy_old=iy_old+yincr;
  }

//ix_old=0;
  ix_old=incr * srates[cmpix] * (MSECS_START[cmpix]-MINI );
  idx = ix_old;
  ix=nusmp[cmpix];
  iy_old=avstand;
  ix_new=0;
  iy_new=avstand;

  for(i=0;i<ix;i++)                          // plot this buffer
  {
    ix_new=(int)idx;
//    iy_new=(int)vraw[cmpix+kl][i]*0.08+avstand;
    iy_new=(int)vflt[cmpix][i]*1.001+avstand;
    gdImageLine(im, ix_old, iy_old, ix_new, iy_new, dg);  // draw from old to new point
    ix_old=ix_new;
    iy_old=iy_new;
    idx=idx+incr;
  }

//-----------------------------------------------
// draw phases this component
//-----------------------------------------------
  if(PPH[cmpix] > 0)
  {
    gdImageSetThickness(im,2);    
    sprintf(dum,"P");
    displ=incr * srates[cmpix] * (PPH[cmpix] - MINI);
    gdImageLine(im, displ, avstand+30, displ, avstand-30, sg1); 
    gdImageLine(im, displ, avstand, displ-30, avstand-30, sg2);
    gdImageString(im,gdFontGetLarge(),displ-10,avstand-23,dum,sg3); 
    gdImageSetThickness(im,1);    
  }
  if(SPH[cmpix] > 0)
  {
    gdImageSetThickness(im,2);    
    sprintf(dum,"S");
    displ=incr * srates[cmpix] * (SPH[cmpix] - MINI);
    gdImageLine(im, displ, avstand+30, displ, avstand-30, sg1);
    gdImageLine(im, displ, avstand, displ-30, avstand-30, sg2);
    gdImageString(im,gdFontGetLarge(),displ-10,avstand-23,dum,sg3); 
    gdImageSetThickness(im,1);    
  }


  idx=0;
  ix_old=0;
  avstand = avstand + yincr;
}


out = fopen (pngfile, "wb");
gdImagePng (im, out);
fclose (out);
gdImageDestroy(im);                   /* Destroy the image in memory. */









}
#endif
Check_SS(char sfilepath[])
{
  FILE *sf;
  FILE *out;
  int no,i,br,l,n;
  int totlines;
  char linjer[500];
  char params1[MAXCHA][500];
  char toplot1[MAXCHA][500];
  int index[500];
  char navn[256];

  int goon;
  int cr;
  int readings_left;
  float residual;
  int ant=0;
  int idx;
  int bytes_read;
  char dummy[500];
  float res[256];
  float res_save[256];
  float max_res;
  sf=fopen(sfilepath,"r");

  if(sf == NULL)
  {
    printf("rtpng:CHECK_S:Can't open S-file: %s\n",sfilepath);
    exit(0);
  }else{
    if(prt >= 1)
      printf("rtpng:Found....................: %s\n",sfilepath);

    br=0;
    cr=0;
    no = 200;
    totlines=0;
    sprintf(navn,"CAZ7");
    while ( fgets ( linjer, sizeof linjer, sf ) != NULL ) /* read a line from file */
    {
//      printf("%s",linjer);
      if (strstr(linjer, navn) )
      {
//	printf("%s", linjer);
        break;
      }
    }

    ant=0;
    while ( fgets ( linjer, sizeof linjer, sf ) != NULL ) /* read a line from file */
    {
     if(prt >= 1)   
        printf("%s",linjer);
      sprintf(HYPFILE[ant],"%s",linjer);
      ant++;
    }
    fclose(sf);
}
no_lines=ant;
if(prt >= 1)
  printf("NO_LINES: %d\n",no_lines); 
return(no_lines);

}
Read_Parameters()
{
  FILE *par;
  char  dummy[256];
  char  linjer[256];
  char  keyword[256];
  float val;
  int   valint;
  int   br = 0;
  
  par = fopen("rtquake.par","r");
  if(par == NULL)
  {
    printf("RTPICK: rtquake.par does not exist in current directory.\n");

    sprintf(dummy,"%s/com/rtquake.par",topdir_rt);
    printf("RTPICK: Look in %s\n",dummy);    
    par = fopen(dummy,"r");
    if(par == NULL)
    {
      printf("RTPICK: rtquake.par does not exist in %s/com. Use program defaults.\n",topdir_rt);
      br=1;      
    }
  }
  if(br == 0)
  {
    printf("RTPICK: \n");
    sprintf(keyword,"LOCATION");
    while ( fgets ( linjer, sizeof linjer, par ) != NULL ) /* read a line from file */
    {
      if (strstr(linjer, keyword) )
      {
	if(prt > 0)
	  printf("%s", linjer);
	sscanf(linjer,"%40c%f",dummy,&val);
//	printf("val: %5.1f\n",val);
	locate = (int)val;
        break;
      }
    }
    rewind(par);
    sprintf(keyword,"AUTOMAG");
    while ( fgets ( linjer, sizeof linjer, par ) != NULL ) /* read a line from file */
    {
  
      if (strstr(linjer, keyword) )
      {
	if(prt > 0)
	  printf("%s", linjer);
	sscanf(linjer,"%40c%f",dummy,&val);
	if(prt > 0)
          printf("val: %5.1f\n",val);
	automag = (int)val;
        break;
      }
    }   
    rewind(par);    
    sprintf(keyword,"KEEP");
    while ( fgets ( linjer, sizeof linjer, par ) != NULL ) /* read a line from file */
    {
  
      if (strstr(linjer, keyword) )
      {
	if(prt > 0)
	  printf("%s", linjer);
	sscanf(linjer,"%40c%f",dummy,&val);
//        printf("val: %5.1f\n",val);
	keep = (int)val;
        break;
      }
    }   
    rewind(par);    
    sprintf(keyword,"ITERATION");
    while ( fgets ( linjer, sizeof linjer, par ) != NULL ) /* read a line from file */
    {
  
      if (strstr(linjer, keyword) )
      {
	if(prt > 0)
	  printf("%s", linjer);
	sscanf(linjer,"%40c%f",dummy,&val);
//        printf("val: %5.1f\n",val);
	iterations = (int)val;
        break;
      }
    }   
    rewind(par);    
    sprintf(keyword,"PRINTING");
    while ( fgets ( linjer, sizeof linjer, par ) != NULL ) /* read a line from file */
    {
  
      if (strstr(linjer, keyword) )
      {
	if(prt > 0)
	  printf("%s", linjer);
	sscanf(linjer,"%40c%f",dummy,&val);
//        printf("val: %5.1f\n",val);
//	prt = (int)val;	
        break;
      }
    }   
    rewind(par);    
    sprintf(keyword,"MAX_RESIDUAL");
    while ( fgets ( linjer, sizeof linjer, par ) != NULL ) /* read a line from file */
    {
  
      if (strstr(linjer, keyword) )
      {
	if(prt > 0)
	  printf("%s", linjer);
	sscanf(linjer,"%40c%f",dummy,&val);
//        printf("val: %5.1f\n",val);
	maxres = val;
//	printf("MAX_RES: %5.2f\n",maxres);
//	exit(0);
        break;
      }
    }       
    rewind(par);
    sprintf(keyword,"MINSTALOC");
    while ( fgets ( linjer, sizeof linjer, par ) != NULL ) /* read a line from file */
    {
  
      if (strstr(linjer, keyword) )
      {
	if(prt > 0)
	  printf("%s", linjer);
	sscanf(linjer,"%40c%d",dummy,&valint);
//        printf("val: %5.1f\n",val);
	no_stations_trg = valint;
//	printf("MINSTALOC  no_stations_trg: %d\n",no_stations_trg);
        break;
      }
    }       
    rewind(par);     
    sprintf(keyword,"MAIL1");
    while ( fgets ( linjer, sizeof linjer, par ) != NULL ) /* read a line from file */
    {
  
      if (strstr(linjer, keyword) )
      {
	if(prt > 0)
	  printf("%s", linjer);
	sscanf(linjer,"%50c%s",dummy,mailaddress1);
	sscanf(linjer,"%40c%f",dummy,&val);	
//        printf("val: %5.1f  %s\n",val,mailaddress1);
	mail1 = (int)val;	
        break;
      }
    }       
    rewind(par);    
    sprintf(keyword,"MAIL2");
    while ( fgets ( linjer, sizeof linjer, par ) != NULL ) /* read a line from file */
    {
  
      if (strstr(linjer, keyword) )
      {
	if(prt > 0)
	  printf("%s", linjer);
	sscanf(linjer,"%50c%s",dummy,mailaddress2);
	sscanf(linjer,"%40c%f",dummy,&val);	
//        printf("val: %5.1f  %s\n",val,mailaddress1);
	mail2 = (int)val;	
        break;
      }
    }      
    rewind(par);    
    sprintf(keyword,"MAIL3");
    while ( fgets ( linjer, sizeof linjer, par ) != NULL ) /* read a line from file */
    {
  
      if (strstr(linjer, keyword) )
      {
	if(prt > 0)
	  printf("%s", linjer);
	sscanf(linjer,"%50c%s",dummy,mailaddress3);
	sscanf(linjer,"%40c%f",dummy,&val);	
//        printf("val: %5.1f  %s\n",val,mailaddress1);
	mail3 = (int)val;	
        break;
      }
    }      
    
    rewind(par);    
    sprintf(keyword,"MAIL4");
    while ( fgets ( linjer, sizeof linjer, par ) != NULL ) /* read a line from file */
    {
  
      if (strstr(linjer, keyword) )
      {
	if(prt > 0)
	  printf("%s", linjer);
	sscanf(linjer,"%50c%s",dummy,mailaddress4);
	sscanf(linjer,"%40c%f",dummy,&val);	
//        printf("val: %5.1f  %s\n",val,mailaddress1);
	mail4 = (int)val;	
        break;
      }
    }        
    
    rewind(par);    
    sprintf(keyword,"MAIL5");
    while ( fgets ( linjer, sizeof linjer, par ) != NULL ) /* read a line from file */
    {
  
      if (strstr(linjer, keyword) )
      {
	if(prt > 0)
	  printf("%s", linjer);
	sscanf(linjer,"%50c%s",dummy,mailaddress5);
	sscanf(linjer,"%40c%f",dummy,&val);	
//        printf("val: %5.1f  %s\n",val,mailaddress1);
	mail5 = (int)val;	
        break;
      }
    }      
    
    
    
    
    
    
    
    
    
    
    
    
    
    
  fclose(par);    
    
  }
}

Check_S(int *ne,float *resid,int *nst,float res_lim,int lp)
{
  FILE *sf;
  FILE *out;
  int no,i,br,l,n;
  int totlines;
  char linjer[200];
  char params1[MAXCHA][50];
  char toplot1[MAXCHA][80];
  int index[200];
  char navn[25];
  char seq[2];
  int goon;
  int cr;
  int readings_left;
  float residual;
  int ant=0;
  int idx;
  int bytes_read;
  char dummy[50];
  float res[200];
  float res_save[200];
  float max_res;
//  sf=fopen(sfilepath,"r");
  sf=fopen("hyp.out","r");
  if(sf == NULL)
  {
    printf("RTPICK :CHECK_S:Can't open S-file: hyp.out\n");
    exit(0);
  }else{
    if(lp==0)                                              // write only once
      printf("RTPICK: Found............................: hyp.out\n");

    br=0;
    cr=0;
    no = 200;
    totlines=0;

    for(i=0;i<10;i++)
    {
      for(l=0;l<90;l++)
      {
        linjer[l]=fgetc(sf);
        if(linjer[l] == 0xa)
          cr++;
        linjer[l]='\0';

        if(cr>=5)
        {
          br=1;
          break;
        }
      }
      if(br==1)
        break;
    }

    br=0;
    cr=0;
    seq[0]=' ';
    seq[1]=' ';
    for(i=0;i<no;i++)
    {
      br=0;
      for(l=0;l<90;l++)
      {
        linjer[l]=fgetc(sf);
        if(linjer[l] == EOF)
        {
          totlines = i-1;
          ant = i-1;
          fclose(sf);
          cr=1;
          break;
        }else{
          if(linjer[l] == 0xa )
          {
            br=1;
            linjer[l]='\0';

            for(n=0;n<20;n++)
              navn[n]='\0';
            for(n=0;n<6;n++)
              navn[n]=linjer[n];
            sprintf(toplot1[i], "%s",navn);
            sprintf(HYPFILE[i],"%s",linjer);

            sscanf(linjer,"%63c %f",dummy,&res[i]);

            res_save[i]=res[i];
            ant++;
          }        // 0xa
        }
        if(br == 1)
          break;
      }    // l=0-90
      if(cr == 1)
        break;
    }  // i=0-no
  }
    fclose(sf);
    no_lines=ant;

  residual=0;
  
  for(l=0;l<ant;l++)
    residual=residual+fabs(res_save[l]);
  residual=residual/ant;
  printf("RTPICK: readings left....................: %2d Avg.res: %8.2f  ",ant,residual);
  *resid = residual;
  *nst   = ant;
  for(l=0;l<200;l++)
    index[l]=0;

  for(l=0;l<1;l++)
  {
    max_res = fabs(res[0]);
    idx=0;
    for(i=1;i<ant;i++)
    {
      if(fabs(res[i]) > max_res)
      {
        max_res = fabs(res[i]);
        idx = i;
      }
    }
    res[idx]=0.0;
    index[idx]=1;
    
  }

  residual=0;
  sf=fopen("hyp.out","r");
  if(sf == NULL)
  {
    printf("RTPICK: CHECK_S: Can't open hyp.out\n");
    exit(0);
  }else{
//    printf("Found the: hyp.out\n");
  }
  out=fopen("hyp_new.out","w");
  if(out == NULL)
  {
    printf("RTPICK: CHECK_S: Can't open hyp_new.out\n");
    exit(0);
  }else{

    br=0;
    cr=0;
    readings_left = 0;
    for(i=0;i<(ant+5);i++)
    {
      for(l=0;l<90;l++)
      {
        linjer[l]=fgetc(sf);
        if(linjer[l] == 0xa)
        {
          cr=1;
          br++;
        }
        if(cr == 1)
        {
          if(br <= 5)
          {
            fprintf(out,"%s",linjer);
            cr=0;
            break;
          }else{
            if(index[i-5] == 0)
            {
              fprintf(out,"%s",linjer);
              readings_left++;
              residual=residual+fabs(res_save[i-5]);
            }
            cr=0;
            break;
          }
        }
      }

    }
    fclose(out);
    fclose(sf);

    ant=ant-1;

    printf("phases left: %2d Avg.residual in HYP_NEW: %6.2f\n",readings_left,residual/ant);
    *ne=ant;

    if( (residual/ant) > res_lim)
    {
      goon=1;  // Continue
    }else{
      printf("RTPICK: STOP iterations. Residual below..: %4.2f\n",res_lim);
      goon=0;  // Stop

    }
    return(goon);




  }

}
/*****************************************************************************/
/* Read copy of 'hyp.out' ('hyp.tmp') and 'hyp_save.out'                     */
/* All NOT USED readings in 'hyp_save.out' are added to 'hyp.tmp' to form    */
/* new complete S-file.                                                      */
/*****************************************************************************/
Check_hyp()
{
  FILE *sf;
  FILE *old;
  int i,br,l,n;
  int totlines;
  char linjer[200];
  char stat[MAXCHA*3][20];
  char comp[MAXCHA*3][20];
  char ph[MAXCHA*3][20];
  char stat0[20];
  char comp0[20];
  char ph0[20];
  char comm[200];
  int r1,r2,r3,ret;
  int ant=0;

  sf=fopen("hyp.tmp","rw+");
  if(sf == NULL)
  {
    printf("CHECK_HYP: Can't open: hyp.tmp\n");
    exit(0);
  }else{
    if(prt >= 1)
      printf("rtpick:Found the................: hyp.tmp file\n");
printf("--------------------------------- s-file ----------------------------------------\n");
    br=0;
    totlines=0;
    for(i=0;i<100;i++)
    {
      for(l=0;l<90;l++)
      {
        linjer[l]=fgetc(sf);
        if(linjer[l] == 0xa)                               // check for end-of-line
        {
          linjer[l] = '\0';
          printf("%s\n",linjer);                           // print line
//	  if(linjer[79] == '1')
//	  {
//            if(mail == 1)
//            {
//	      sprintf(comm,"echo %c%s%c | sendmail %s",0x22,linjer,0x22,mailaddress1);
//	      ret=system(comm);
//            }
//	  }
          if(linjer[79] == 0x37)                           // test for last header line
            br=1;
          break;
        }
      }
      if(br==1)                                            // last line in header....."7" line
        break;
    }
    br=0;
    for(i=0;i<1000;i++)
    {
      for(l=0;l<90;l++)
      {
        linjer[l]=fgetc(sf);
        if(linjer[l] == 0xa)                               // check for end-of-line
        {
          linjer[l] = '\0';
          if(linjer[1]!=' ' && linjer[2]!=' ')             // last line ??
          {
            printf("%s\n",linjer);                         // print line
            totlines++;
// save stations and components used for location
            for(n=0;n<5;n++)
              stat[i][n]=linjer[n+1];
            stat[i][5]='\0';
            comp[i][0]=linjer[6];
            comp[i][1]=linjer[7];
            comp[i][2]='\0';
            ph[i][0]=linjer[9];
            ph[i][1]=linjer[10];
            ph[i][2]='\0';
            break;
          }
        }
      }
    }
printf("--------------------------------- s-file ----------------------------------------\n");
//    for(n=0;n<totlines;n++)
//      printf("%s %s %s\n",stat[n],comp[n],ph[n]);

    old=fopen("hyp_all.out","r");                         // read original hyp.out (hyp_all.out)
    if(old == NULL)
    {
printf("rtpick:CHECK_HYP................: Can't open: hyp_all.out\n");
      exit(0);
    }else{
printf("rtpick:Found the................: hyp_all.out file\n");
      br=0;
      for(i=0;i<100;i++)
      {
        for(l=0;l<90;l++)
        {
          linjer[l]=fgetc(old);
          if(linjer[l] == 0xa)                             // check for end-of-line
          {
            linjer[l] = '\0';
//            printf("%s\n",linjer);                         // print line

            if(linjer[79] == 0x37)                         // test for last header line
              br=1;
            break;
          }
        }
        if(br==1)                                          // last line in header....."7" line
          break;
      }
//      printf("last line in header hyp_save.out\n");
/**************************************************************************/
/* Include to avoid blank line in s-file after readings used for location */
/*       fseek(sf,-81,SEEK_END);                                          */
/**************************************************************************/
      br=0;
      for(i=0;i<1000;i++)
      {
//        printf("read line %d\n",i);
        for(l=0;l<90;l++)
        {
          linjer[l]=fgetc(old);
          if(linjer[l] == 0xa)                             // check for end-of-line
          {
            linjer[l] = '\0';
            if(linjer[1]!=' ' && linjer[2]!=' ')           // last line ?
            {
//              printf("%s\n",linjer);                       // check against the other

// save station and component in original s-file
              for(n=0;n<5;n++)
                stat0[n]=linjer[n+1];
              stat0[5]='\0';
              comp0[0]=linjer[6];
              comp0[1]=linjer[7];
              comp0[2]='\0';
              ph0[0]=linjer[9];
              ph0[1]=linjer[10];
              ph0[2]='\0';

              for(n=0;n<totlines;n++)                      // write readings not used in location at end of new s-file
              {
                r1=strncmp(stat0,stat[n],5);
                r2=strncmp(comp0,comp[n],2);
                r3=strncmp(ph0,ph[n],2);
//                printf("%s %s  %s %s r1:%5d r2:%5d\n",stat0,comp0,stat[n],comp[n],r1,r2);
                if(r1 == 0 && r2 == 0 && r3 == 0)          // reading exists, do not write
                {
                  ant=1;
                  break;
                }
//                getchar();
              }
              if(ant==1)
              {
//                printf("do not write\n");
                ant=0;
              }else{
                ant=0;
//                printf("%s\n",linjer);
                linjer[14]='4';                     // 0 % weight
                ret=fprintf(sf,"%s\n",linjer);
//                printf("ret= %d\n",ret);
              }
              break;
            }
          }
        }
      }



    }
  }
fclose(sf);
}

Update_map_file(char snam[],int ant,float resid,int nst,int klon, char sfilepath[])
{
  FILE *sf;
  FILE *loc;
  FILE *out;
  FILE *gog;
  FILE *llc;
  FILE *gmt;
  int i,l,cr,br,tel;
  int karakterer;
  int dagen;
  int timen;
  int minuttet;
  int sekundet;
  int aret;
  int maned;
  float lat,lon;
  char linjer[200];
  char dummy[500];
  char mapname[200];
  char locname[200];
  char mapout[200];
  char gmtmap[200];
  char name[200];
  char copyhyp[200];
  char lastline[200];
  float res[200];
  float mag;
  long size;
  float avg;
  int ret;
  tel=0;
  avg=0;
  if(prt >= 1)  
    printf("Update_map_file: %s\n",snam);
  for(i=0;i<256;i++)
  {
    if(snam[i] == '.' && snam[i+1] =='S')
      break;
  }

  for(l=0;l<20;l++)
    snam[l]=snam[i-11+l];
  
  sf=fopen("hyp.out","r");
  if(sf == NULL)
  {
    printf("Update_map_file: Can't open hyp.out\n");
    exit(0);
  }else{
//    printf("Found the: hyp.out\n");
  }
  
  sprintf(mapout,"%s/map/map.out",topdir_rt); 
  out=fopen(mapout,"r+");
  if(out == NULL)
  {
    printf("rtpick:Update_map_file: Can't open: %s\n",mapout);
    exit(0);
  }else{
    br=0;
    cr=0;
    for(i=0;i<(ant+5);i++)
    {
      for(l=0;l<90;l++)
      {
        linjer[l]=fgetc(sf);
        if(linjer[l] == 0xa)
        {
          cr=1;
          br++;
          if(br==1)
          {
            sscanf(linjer,"%23c %f %f",dummy,&lat,&lon);
            sscanf(linjer,"%56c%f3.1",dummy,&mag);
          }
          if(br==3)
          {
            sscanf(linjer,"%s",name);
          }         
        }
        if(cr == 1)
        {
          if(br <= 5)
          {
            cr=0;
            break;
          }else{
            sscanf(linjer,"%65c %f",dummy,&res[tel]);

            avg=avg+fabs(res[tel]);
            tel++;
            cr=0;
            break;
          }
        }
      }

    }
    fseek(out,0,SEEK_END);
    size=ftell(out);
    fseek(out,size-80,SEEK_SET);
    for(l=0;l<80;l++)
    {
      linjer[l]=fgetc(out);
      if(linjer[l]==0xa)
      {
        fseek(out,0,SEEK_CUR);
        fprintf(out,"circle,red,%s %6.2f %2d,%8.3f,%8.3f\n",snam,resid,nst,lat,lon);
        fprintf(out,"google,yellow,%s %6.2f %2d,%8.3f,%8.3f\n",snam,resid,nst,lat,lon);
        break;
      }
    }
    fclose(out);
    fclose(sf);
  }
  sprintf(mapname,"%s/map/ALL_EPI%d.txt",topdir_rt,klon);
  gog=fopen(mapname,"a+");
  if(gog == NULL)
  {
    printf("rtpick:Update_map_file: Can't open %s\n",mapname);
    exit(0);
  }else{
        printf("RTPICK: MAG..............................: %3.1f\n",mag);
        fprintf(gog,"%s | %6.2f | %8.3f | %8.3f | %2d | %4.1f | text\n",snam,resid,lat,lon,nst,mag);
        fclose(gog);
  }
// write cdate, coordinates,depth and magnitude for gmt map 
  sscanf(snam,"%2d%1c%2d%2d%1c%2d%3c%4d%2d",&dagen,dummy,&timen,&minuttet,dummy,&sekundet,dummy,&aret,&maned);
  sprintf(gmtmap,"%s/map/GMT_EPI%d.txt",topdir_rt,klon);
  gmt=fopen(gmtmap,"a+");
  if(gmt == NULL)
  {
    printf("rtpick:Update_map_file: Can't open %s\n",gmtmap);
    exit(0);
  }else{
        fprintf(gmt,"%4d %2d %2d %2d %2d %2d %6.2f  %8.3f %2d %4.1f\n",aret,maned,dagen,timen,minuttet,sekundet,lon,lat,nst,mag);
        fclose(gmt);
  }
  
  
    sprintf(mapname,"%s/loc/%s.html",topdir_rt,snam);
    loc=fopen(mapname,"w");
    fprintf(loc,"<!DOCTYPE html>\n");
    fprintf(loc,"<html>\n");
    fprintf(loc,"<body>\n");
    fprintf(loc," <img src=%chttp://maps.googleapis.com/maps/api/staticmap?center=%8.3f,%8.3f&zoom=7&size=600x600&maptype=roadmap&markers=color:blue|label:S|40.702147,-74.015794&markers=color:red|color:red|label:C|%8.3f,%8.3f&sensor=false%c>\n",0x22,lat,lon,lat,lon,0x22);
    fprintf(loc,"</body>\n");
    fprintf(loc,"</html>\n");
    fclose(loc);
    printf("RTPICK: mail1............................: %d\n",mail1);
    
    sprintf(copyhyp,"cp hyp.out hyp.txt");
    ret=system(copyhyp);
      for(l=0;l<500;l++)
        mail_message[l]='\0';

//sprintf(mail_message,"echo %chttp://maps.googleapis.com/maps/api/staticmap?center=%f,%f&zoom=7&size=900x1000&maptype=hybrid&markers=color:red_7Ccolor:red_7Clabel:Q_7C%f,%f&sensor=false%c | mutt  -a %c%s/loc/ALL.png%c -s %cTRIGGER%c -- %s &",0x22,lat,lon,lat,lon,0x22,0x22,topdir_rt,0x22,0x22,0x22,mailaddress1);
sprintf(mail_message,"echo %chttp://maps.googleapis.com/maps/api/staticmap?center=%f,%f&zoom=7&size=900x1000&maptype=hybrid&markers=color:red_7Ccolor:red_7Clabel:Q_7C%f,%f&sensor=false \n %s%c\
		| mutt  -a %c%s/loc/ALL.png%c -s %cLOCATED EVENT: M %.1f%c -a %chyp.out%c -- %s &",0x22,lat,lon,lat,lon,sfilepath,0x22,0x22,topdir_rt,0x22,0x22,mag,0x22,0x22,0x22,mailaddress1);

//printf("%s\n",mail_message);
      for(l=0;l<500;l++)
      {
        if(mail_message[l]=='_')
          mail_message[l]='%';
      }

    
/*    
    if(mail1 == 1)
    {
      for(l=0;l<500;l++)
        mail_message[l]='\0';

//sprintf(mail_message,"echo %chttp://maps.googleapis.com/maps/api/staticmap?center=%f,%f&zoom=7&size=900x1000&maptype=hybrid&markers=color:red_7Ccolor:red_7Clabel:Q_7C%f,%f&sensor=false%c | mutt  -a %c%s/loc/ALL.png%c -s %cTRIGGER%c -- %s &",0x22,lat,lon,lat,lon,0x22,0x22,topdir_rt,0x22,0x22,0x22,mailaddress1);
sprintf(mail_message,"echo %chttp://maps.googleapis.com/maps/api/staticmap?center=%f,%f&zoom=7&size=900x1000&maptype=hybrid&markers=color:red_7Ccolor:red_7Clabel:Q_7C%f,%f&sensor=false%c | mutt  -a %c%s/loc/ALL.png%c -s %cTRIGGER%c -a %chyp.txt%c -- %s &",0x22,lat,lon,lat,lon,0x22,0x22,topdir_rt,0x22,0x22,0x22,0x22,0x22,mailaddress1);
printf("%s\n",mail_message);
      for(l=0;l<500;l++)
      {
        if(mail_message[l]=='_')
          mail_message[l]='%';
      }

    }
    if(mail2 == 1)
    {
      for(l=0;l<500;l++)
        mail_message[l]='\0';
          sprintf(mail_message,"echo %cQUAKE%c | mutt -a %c%s/rtquake/map/hyp.out%c -s %cTRIGGER%c -- %s &",0x22,0x22,0x22,topdir_rt,0x22,0x22,0x22,mailaddress1);

sprintf(mail_message,"echo %chttp://maps.googleapis.com/maps/api/staticmap?center=%f,%f&zoom=7&size=900x1000&maptype=hybrid&markers=color:red_7Ccolor:red_7Clabel:Q_7C%f,%f&sensor=false%c | mutt  -a %c%s/loc/ALL.png%c -s %cTRIGGER%c -- %s &",0x22,lat,lon,lat,lon,0x22,0x22,topdir_rt,0x22,0x22,0x22,mailaddress2);
printf("%s\n",mail_message);
      for(l=0;l<500;l++)
      {
        if(mail_message[l]=='_')
          mail_message[l]='%';
      }

    }    
*/    
    
    
    
  sprintf(locname,"%s/map/LAST_LOC.txt",topdir_rt);
  llc=fopen(locname,"w");
  if(llc == NULL)
  {
    printf("rtpick:Update_map_file: Can't open %s\n",locname);
    exit(0);
  }else{
    fprintf(llc,"%s | %6.2f | %8.3f | %8.3f | %2d | %4.1f | text\n",snam,resid,lat,lon,nst,mag);
    fclose(llc);
  }

  

}



Create_Sfile(char filename[],char sfilepath[],char dbname[],int nchannels,int klon,int locate,char yr_trg[],char mo_trg[],char da_trg[],char hr_trg[],char mi_trg[],char se_trg[],char sfilename[])
{
  FILE *sf;
  FILE *so1;
  FILE *so2;
  int retur;
  int i;
  char fullpath[200];
  char syscmd[200];
  char comm0[200];
  char comm1[200];
  char comm2[200];
  char comm3[200];
  char comm4[200];
  char comm5[200];
  char comm6[200];
  char comm7[200];
  char comm8[200];
  char comm9[200];
  char com10[200];
  char com11[200];
  char com12[200];
  char com13[200];  
  char com14[200];
  char com15[200];  
  
  char com20[200];
  char com21[200];
  
  char snam[200];
  char dummy[200];
  char fpath[200];
  char year[5];
  char month[5];
  char day[5];
  char hour[5];
  char minute[5];
  char sec[5];
  
  char record[256];

  char seconds[10];
  char tempo[80];
  int ret;
  int lp;
  int kart;
  int ne;
  int rettur;
  int br;

  float resid;
  float res_lim;
  int nst;
  for(i=0;i<200;i++)
  {
    snam[i]='\0';
    fullpath[i]='\0';
  }


  printf("RTPICK: Path+s-filename..................: %s\n",sfilepath);

  S_REC(yr_trg,mo_trg,da_trg,hr_trg,mi_trg, se_trg, dummy,sfilepath,filename,nchannels,snam);
//---------------------------------------------
// check for auto locate
//---------------------------------------------

  if(keep == 1 && locate == 1)
  {
    printf("RTPICK: Create_Sfile.....................: Locate + new s-file.\n");

    sprintf(com13,"automag %s",sfilepath);
    sprintf(com14,"hyp automag.out >> hyptemp.txt");
    sprintf(com15,"cp automag.out %s",sfilepath);

    sprintf(comm0,"rm hyptemp.txt");
    sprintf(comm1,"hyp %s >> hyptemp.txt",sfilepath);
    sprintf(comm2,"cp hyp.out %s", sfilepath);
    sprintf(comm3,"cp hyp_new.out %s", sfilepath);
    sprintf(comm4,"cp %s hyp_save.out",sfilepath);
    sprintf(comm5,"cp hyp_save.out %s",sfilepath);
    sprintf(comm6,"cp hyp.out hyp.tmp");
    sprintf(comm7,"cp hyp.tmp %s",sfilepath);
    sprintf(comm8,"cp hyp.out hyp_all.out");
    sprintf(comm9,"cp hyp_save.out %s",sfilepath);
    sprintf(com10,"cp %s s_org.out",sfilepath);
#ifdef RTQUAKE
    sprintf(com11,"cp hyp.out %s/map",topdir_rt);
#endif    
    printf("RTPICK: comm0............................: %s\n",comm0);
    ret=system(comm0);                                     // remove old hyptemp.txt
 
    printf("RTPICK: com10............................: %s\n",com10);
    ret=system(com10);                                     // copy original s-file to s_org.out

    printf("RTPICK: comm1............................: %s\n",comm1);
    ret=system(comm1);                                     // run hyp, input s-file, this produces "hyp.out"

    printf("RTPICK: comm2............................: %s\n",comm2);
    ret=system(comm2);                                     // copy hyp.out new s-file

    printf("RTPICK: comm8............................: %s\n",comm8);
    ret=system(comm8);                                     // copy hyp.out to hyp_all.out

    res_lim = maxres;
    kart=0;
    rettur=0;
    for(lp=0;lp<iterations;lp++)
    {
      ret=Check_S(&ne,&resid,&nst,res_lim,lp);       // read 'hyp.out', create 'hyp_new.out'
      
//printf("ret: %d ne: %d resid: %f6.2 nst: %d\n",ret,ne,resid,nst);      
      
      if(ret==1)
      {
        ret=system(comm3);                                 // copy hyp_new to new s-file 
        ret=system(comm0);                                 // remove old hyptemp.txt
        ret=system(comm1);                                 // run hyp, input s-file
        ret=system(comm2);                                 // copy hyp.out new s-file
      }else{
        printf("RTPICK: Average residual.................: %f\n",resid);
        printf("RTPICK: No more iterations...............: Number of stations: %d Avg: res.: %8.3f\n",nst,resid);
        if(  nst < no_stations_trg)
        {
          printf("RTPICK: Too few stations (<%d)............: No location\n",no_stations_trg);
	  sendmail=1;
          ret=system(comm5);
          kart=1;
          rettur = -1;
        }else if(  nst >= no_stations_trg)
	{
	  sendmail=2;
	}
        break;
      }
    }
    if(rettur==0)
    {
      ret=system(comm6);                                   // copy last hyp.out to hyp.tmp
      printf("RTPICK: comm6............................: %s\n",comm6);
#ifdef RTQUAKE
      ret=system(com11);
      printf("RTPICK: com11............................: %s\n",com11);
#endif
      printf("RTPICK: comm2............................: %s\n",comm2);
      ret=system(comm2);
    }else{
      ret=system(comm9);
      printf("RTPICK: comm9............................: %s\n",comm9);
    }
/*------------------------------------------------------------------*/
//    sprintf(com20,"gnuplot -e %cfilename='res.out'%c /home/seismo/resplot2",0x22,0x22);   
//    sprintf(com21,"resiplt");
    
//    ret=system(com21);
//    ret=system(com20);    
    sprintf(com20,"respng");
    ret= system(com20);    
    
/*-------------------------------------------------------------------*/    
    if(kart==0)
    {
      int rss=resid*1000;
      printf("RTPICK: RSS..............................: %d\n",rss);
#ifdef RTQUAKE      
      if(rss != 0)
      {
        printf("RTPICK:..................................: Update map\n");
        Update_map_file(sfilename,ne,resid,nst,klon,sfilepath);

	if(automag == 1)
	{
//--------------AUTOMAG------------------------------------------------- 
          printf("RTPICK: AM:comm0.........................: %s\n",comm0);
          ret=system(comm0);                                     // remove old hyptemp.txt
          printf("RTPICK: AM:com10.........................: %s\n",com10);
          ret=system(com10);                                     // copy original s-file to s_org.out
          printf("RTPICK: AM:com13.........................: %s\n",com13);    
          ret=system(com13);                                     // run automag, creates automag.outfile
          printf("RTPICK: AM:com15.........................: %s\n",com15);
          ret= system(com15);                                    // copy automag.out new s-file

          printf("RTPICK: AM:com14.........................: %s\n",com14);
          ret=system(com14);                                     // run hyp, input automag.out, this produces "hyp.out"       
          printf("RTPICK: AM:comm2.........................: %s\n",comm2);
          ret=system(comm2);                                     // copy hyp.out new s-file
          printf("RTPICK: AM:comm8.........................: %s\n",comm8);
          ret=system(comm8);                                     // copy hyp.out to hyp_all.out
//--------------END AUTOMAG---------------------------------------------   
	}
      }

#else
      printf("Plotting static map is not activated\n");
#endif      
    }
  }else{
    
/*    
    
    
    if(locate == 1 && keep == 1)
  {
    printf("Create_Sfile: Locate + keep old s-file + add phases after n iterations.\n");
    sprintf(comm0,"rm hyptemp.txt");
    sprintf(comm1,"hyp %s >> hyptemp.txt",sfilepath);
    sprintf(comm2,"cp hyp.out %s", sfilepath);
    sprintf(comm3,"cp hyp_new.out %s", sfilepath);
    sprintf(comm4,"cp %s hyp_save.out",sfilepath);
    sprintf(comm5,"cp hyp_save.out %s",sfilepath);
    sprintf(comm6,"cp hyp.out hyp.tmp");
    sprintf(comm7,"cp hyp.tmp %s",sfilepath);
    sprintf(comm8,"cp hyp.out hyp_all.out");
    sprintf(comm9,"cp hyp_save.out %s",sfilepath);
    sprintf(com10,"cp %s s_org.out",sfilepath);
#ifdef RTQUAKE
    sprintf(com11,"cp hyp.out %s/map",topdir_rt);
#endif    
    printf("rtpick:comm0....................: %s\n",comm0);
    ret=system(comm0);                                     // remove old hyptemp.txt
 
    printf("rtpick:com10....................: %s\n",com10);
    ret=system(com10);                                     // copy original s-file to s-sfile

    printf("rtpick:comm1....................: %s\n",comm1);
    ret=system(comm1);                                     // run hyp, input s-file, this produces "hyp.out"

    printf("rtpick:comm2....................: %s\n",comm2);
    ret=system(comm2);                                     // copy hyp.out new s-file

    printf("rtpick:comm8....................: %s\n",comm8);
    ret=system(comm8);                                     // copy hyp.out to hyp_all.out

    res_lim = maxres;
    kart=0;
    rettur=0;
    for(lp=0;lp<iterations;lp++)
    {
      ret=Check_S(&ne,&resid,&nst,res_lim,lp);       // read 'hyp.out', create 'hyp_new.out'
      if(ret==1)
      {
        ret=system(comm3);                                 // copy hyp_new to new s-file 
        ret=system(comm0);                                 // remove old hyptemp.txt
        ret=system(comm1);                                 // run hyp, input s-file
        ret=system(comm2);                                 // copy hyp.out new s-file
      }else{
        printf("rtpick:Average residual.........: %f\n",resid);
        printf("rtpick:No more iterations.......: Number of stations: %d Avg: res.: %8.3f\n",nst,resid);
        if(  nst < 5)
        {
          printf("rtpick:Too few stations (<5)....: No location\n");
          ret=system(comm5);
          kart=1;
          rettur = -1;
        }
        break;
      }
    }
    
    printf("RETTUR: %d\n",rettur);
    if(rettur == 0)
    {
      so1=fopen("hyp_save.out","r");       // open original s-file
      if(so1 == NULL)
      {
        printf("Create_Sfile: Can't open original s-file: hyp_save.out\n");
        exit(0);
      }
      so2=fopen("hyp.out","r");       // open original s-file
      if(so2 == NULL)
      {
        printf("Create_Sfile: Can't open last s-file: hyp.out\n");
        exit(0);
      }   
      
      sf=fopen("comb_sfile","w");
      if(sf == NULL)
      {
        printf("Create_Sfile: Can't open combined S-file: comb_sfile\n");
        exit(0);
      }else{
        for(i=0;i<200;i++)
        {
          br = les_linje(so1,record);  // read first line in new station, and save samplerate + timestamp
          printf("%s\n",record);
//    if(record[79] == '7')
          if(br == 1)
          {
//      fprintf(sf,"%s\n",record);
            fclose(so1);
            break;
          }else{
            fprintf(sf,"%s\n",record);
          }
        }
        for(i=0;i<200;i++)
        {
          br = les_linje(so2,record);  // read first line in new station, and save samplerate + timestamp
          printf("%s\n",record);
          if(record[79] == '7')
	    break;
        }
        for(i=0;i<200;i++)
        {
          br = les_linje(so2,record);  // read first line in new station, and save samplerate + timestamp
          printf("%s\n",record);
//    if(record[79] == '7')
          if(br == 1)
          {
//      fprintf(sf,"%s\n",record);
            fclose(so2);
            break;
          }else{
            fprintf(sf,"%s\n",record);
          }
        }
        fclose(sf);
        sprintf(com12,"cp -f comb_sfile %s",sfilepath);
	printf("COM12: %s\n",com12);
        ret=system(com12);        
      
        
        
      }   // else
    }else{
      ret=system(comm9);
      printf("RTPICK:comm9........................: %s\n",comm9);
    }

    if(kart==0)
    {
      int rss=resid*1000;
      printf("RTPICK: RSS..........................: %d\n",rss);
#ifdef RTQUAKE      
      if(rss != 0)
      {
        printf("RTPICK: Update map\n");
        Update_map_file(sfilename,ne,resid,nst,klon);
//        if(mail == 1)
//        {
//          sprintf(dummy,"echo %cQUAKE%c | mutt -a %c%s/rtquake/map/hyp.out%c -s %cTRIGGER%c -- %s &",0x22,0x22,0x22,topdir_rt,0x22,0x22,0x22,mailaddress1);
//	  printf("%s\n",dummy);
//          ret=system(dummy);
//        }	
      }
#else
      printf("Plotting static map is not activated\n");
#endif      
    }
  }
  
  
*/  
  
  
 }

}
//----------------------------------------------------------------------
// create directory given in 'buffer' if not present
//----------------------------------------------------------------------
CreateDirectory(char* buffer) 
{ 
  int ret;
  int i;
  int siz;
  char path2[80];
  char message[100];
  DIR  *dirp; 
  struct dirent *direntp; 

  siz=strlen(buffer);
  buffer[siz]='/';
  buffer[siz+1]='\0';
  for(i=0;i<siz+1;i++)
  {
    if(buffer[i]=='/' && i > 0)
    {
      strncpy(path2,buffer,i);
      path2[i]='\0';
      dirp = opendir(path2);              // check if dir exists 
      if(dirp == NULL) 
      { 
//--------------------log message-------------------------------
//        sprintf(message,"%d %s does not exist. Create it !",klon,path2);
//        ErrMsg();
        ret=mkdir(path2,S_IRWXU | S_IRGRP | S_IXGRP | S_IROTH | S_IXOTH);
//        if(processor == 0)
//          changemod(path2);
//      }else{
//--------------------log message-------------------------------
//        sprintf(message,"CreateDirectory: '%s' exist.",buffer);
//        ErrMsg();
      } 
      closedir(dirp);
    }
  }
} 

S_REC(year,month,day,hour,minute,seconds,dummy,fullpath,filename,nchannels,snam)
char year[];
char month[];
char day[];
char hour[];
char minute[];
char seconds[];
char dummy[];
char fullpath[];
char filename[];
int nchannels;
char snam[];
{
  FILE    *sf;
  FILE    *so;
  char buf[100];
  char fname[80];
  char def[80];
  char comp[20];
  char record[200];
  int hr1;
  int mn1;
  int sc1;
  int ms1;
  char dm1[50];
  char dm2[50];
  int i;
  int ch; 
  int cnt;
  int br;
  int YR,MTH,DAY,HR,MIN;
  float SEC;
  int AR;

//printf("SECONDS: %s\n",seconds);
  sscanf(year,"%d",&YR);
//printf("%s %d\n",year,YR);
  sscanf(month,"%d",&MTH);
//printf("%s %d\n",month,MTH);
  sscanf(day,"%d",&DAY);
//printf("%s %d\n",day,DAY);
  sscanf(hour,"%d",&HR);
//printf("%s %d\n",hour,HR);
  sscanf(minute,"%d",&MIN);
//printf("%s %d\n",minute,MIN);
  sscanf(seconds,"%f",&SEC);
//printf("%s %f\n",seconds,SEC);

  printf("RTPICK: S_REC: fullpath..................: %s\n",fullpath);

   
/*
  if(keep == 1)
  {
    printf("RTPICK: S_REC: Read old readings from: s_original.sav\n");    
    so=fopen("s_orginal.sav","r");       // open original s-file
    if(so == NULL)
    {
      printf("RTPICK: S_REC: Can't open original s-file: s_original.sav.\n");
      exit(0);
    }
  }
  printf("RTPICK: S_REC: Create s-file.............: %s\n",fullpath);
  sf=fopen(fullpath,"w");
  if(sf == NULL)
  {
    printf("RTPICK: S_REC: Can't open S-file: %s\n",fullpath);
    for(i=0;i<52;i++)
      printf("%2d %2x %c\n",i,fullpath[i],fullpath[i]);
    exit(0);
  }
  if(keep==1)
  {
    printf("RTPICK: S_REC..............: Copy s_original.sav on top in new s-file.\n");
//----------------------------------------------------------------------
// read old s-file and copy first lines to new s-file
//----------------------------------------------------------------------
    for(i=0;i<200;i++)
    {
      br = les_linje(so,record);  // read first line in new station, and save samplerate + timestamp
      printf("%s\n",record);
//    if(record[79] == '7')
      if(br == 1)
      {
//      fprintf(sf,"%s\n",record);
        fclose(so);
        break;
      }else{
        fprintf(sf,"%s\n",record);
      }
    }
  }
*/  
  if(keep == 1)
  {
    printf("RTPICK: S_REC............................: Write new s-file header to s-file.\n");
    
    printf("RTPICK: S_REC: Create s-file.............: %s\n",fullpath);
    sf=fopen(fullpath,"w");
    if(sf == NULL)
    {
      printf("RTPICK: S_REC: Can't open S-file: %s\n",fullpath);
      for(i=0;i<52;i++)
        printf("%2d %2x %c\n",i,fullpath[i],fullpath[i]);
      exit(0);
    }    
    
    for(i=0;i<100;i++)
    buf[i]=' ';
//  printf(" %4d %2d%2d %2d%2d %4.1f\n",YR,MTH,DAY,HR,MIN,SEC);
    sprintf(buf," %4d %2d%2d %2d%2d %4.1f",YR,MTH,DAY,HR,MIN,SEC);
    for(i=0;i<100;i++)
      if(buf[i]=='\0')buf[i]=' ';
    buf[21]='L';
    buf[22]='M';
    buf[45]='T';   // extract Agency from short network name
    buf[46]='S';
    buf[47]='T';
    buf[79]='1';
    buf[80]='\0';
    fprintf(sf,"%s\n",buf);                     // write '1' line sfile
//    printf("%s\n",buf);
    for(i=0;i<100;i++)
      buf[i]=' ';
//printf("In S_REC: FILENAME: %s\n",filename);
    sprintf(buf," %s",filename);
    for(i=0;i<100;i++)
      if(buf[i]=='\0')buf[i]=' ';
    buf[79]='6';
    buf[80]='\0';
    fprintf(sf,"%s\n",buf);

    sprintf(dm1,"%4d",YR);
    sscanf(dm1,"%2c%d",dm2,&AR); 

    for(i=0;i<100;i++)
      buf[i]=' ';
    sprintf(buf," ACTION:NEW %2d-%2d-%2d %2d:%2d OP:SEIS STATUS:               ID:%4d%2d%2d%2d%2d%2d",AR,MTH,DAY,HR,MIN,YR,MTH,DAY,HR,MIN,(int)SEC);

    if(buf[12]==' ')buf[12]='0';
    if(buf[15]==' ')buf[15]='0';
    if(buf[18]==' ')buf[18]='0';
    if(buf[21]==' ')buf[21]='0';
    if(buf[24]==' ')buf[24]='0';
  
    for(i=0;i<100;i++)
      if(buf[i]=='\0')buf[i]=' ';
    for(i=59;i<73;i++)
      if(buf[i]==' ')buf[i]='0';
    buf[79]='I';
    buf[80]='\0';
    fprintf(sf,"%s\n",buf);
//    printf("%s\n",buf);                 
    sprintf(buf," STAT SP IPHASW D HRMM SECON CODA AMPLIT PERI AZIMU VELO SNR AR TRES W  DIS CAZ7");
    fprintf(sf,"%s\n",buf);
//    printf("%s\n",buf);
// STAT SP IPHASW D HRMM SECON CODA AMPLIT PERI AZIMU VELO AIN AR TRES W
// BER  BZ IP       2223  6.00    4
  }
  printf("S_REC: Write new phases to s-file.\n");
  for(ch=0;ch<nchannels;ch++)
  {
//printf("%s\n",trig_components[ch]);
    comp[0]=trig_components[ch][0];
    comp[1]=trig_components[ch][2];
    comp[2]='\0';

    cnt=0;
    for(i=0;i<6;i++)
    {
      if(trig_stations[ch][i] != '\0')
      {
        cnt++;
      }else{
        break;
      }
    }
//printf("%s  CNT: %d\n",trig_stations[ch],cnt);
//      printf("%3d TRIGGER_TIMES: %s %s\n",ch,TRIGGER_TIMES_P[ch],TRIGGER_TIMES_S[ch]);
    hr1=0;
    mn1=0;
    sc1=0;
    ms1=0;
    sscanf(TRIGGER_TIMES_P[ch],"%9c%2d",dm1,&hr1);
    sscanf(TRIGGER_TIMES_P[ch],"%12c%2d",dm1,&mn1);
    sscanf(TRIGGER_TIMES_P[ch],"%15c%2d",dm1,&sc1);
    sscanf(TRIGGER_TIMES_P[ch],"%18c%2d",dm1,&ms1);
    if(TRIGGER_TIMES_P[ch][3] != ' ')
    {
//printf("%3d %3d %3d %3d %3d\n",ch,hr1,mn1,sc1,ms1);
      switch(cnt)
      {
        case 2:
        fprintf(sf," %s   %s IP    A  %2d%2d %2d.%2d %4d\n",trig_stations[ch],comp,hr1,mn1,sc1,ms1,varighet[ch]);
        printf(" %s   %s IP    A  %2d%2d %2d.%2d %4d\n",trig_stations[ch],comp,hr1,mn1,sc1,ms1,varighet[ch]);
        break;
        case 3:
        fprintf(sf," %s  %s IP    A  %2d%2d %2d.%2d %4d\n",trig_stations[ch],comp,hr1,mn1,sc1,ms1,varighet[ch]);
        printf(" %s  %s IP    A  %2d%2d %2d.%2d %4d\n",trig_stations[ch],comp,hr1,mn1,sc1,ms1,varighet[ch]);
        break;
        case 4:
        fprintf(sf," %s %s IP    A  %2d%2d %2d.%2d %4d\n",trig_stations[ch],comp,hr1,mn1,sc1,ms1,varighet[ch]);
        printf(" %s %s IP    A  %2d%2d %2d.%2d %4d\n",trig_stations[ch],comp,hr1,mn1,sc1,ms1,varighet[ch]);
        break;
        case 5:
        fprintf(sf," %s%s IP    A  %2d%2d %2d.%2d %4d\n",trig_stations[ch],comp,hr1,mn1,sc1,ms1,varighet[ch]);
        printf(" %s%s IP    A  %2d%2d %2d.%2d %4d\n",trig_stations[ch],comp,hr1,mn1,sc1,ms1,varighet[ch]);
        break;

//        break;
      }
    }
  }

  for(ch=0;ch<nchannels;ch++)
  {
    comp[0]=trig_components[ch][0];
    comp[1]=trig_components[ch][2];
    comp[2]='\0';

    cnt=0;
    for(i=0;i<6;i++)
    {
      if(trig_stations[ch][i] != '\0')
      {
        cnt++;
      }else{
        break;
      }
    }
    hr1=0;
    mn1=0;
    sc1=0;
    ms1=0;
//      printf("TRIGGER_TIMES_S: %s\n",TRIGGER_TIMES_S[ch]);
    sscanf(TRIGGER_TIMES_S[ch],"%9c%2d",dm1,&hr1);
    sscanf(TRIGGER_TIMES_S[ch],"%12c%2d",dm1,&mn1);
    sscanf(TRIGGER_TIMES_S[ch],"%15c%2d",dm1,&sc1);
    sscanf(TRIGGER_TIMES_S[ch],"%18c%2d",dm1,&ms1);
    if(TRIGGER_TIMES_S[ch][3] != ' ')
    {
//printf("%3d %3d %3d %3d %3d\n",ch,hr1,mn1,sc1,ms1);
      switch(cnt)
      {
        case 2:
        fprintf(sf," %s   %s IS   3A  %2d%2d %2d.%2d\n",trig_stations[ch],comp,hr1,mn1,sc1,ms1);
        printf(" %s   %s IS   3A  %2d%2d %2d.%2d\n",trig_stations[ch],comp,hr1,mn1,sc1,ms1);
        break;
        case 3:
        fprintf(sf," %s  %s IS   3A  %2d%2d %2d.%2d\n",trig_stations[ch],comp,hr1,mn1,sc1,ms1);
        printf(" %s  %s IS   3A  %2d%2d %2d.%2d\n",trig_stations[ch],comp,hr1,mn1,sc1,ms1);
        break;
        case 4:
        fprintf(sf," %s %s IS   3A  %2d%2d %2d.%2d\n",trig_stations[ch],comp,hr1,mn1,sc1,ms1);
        printf(" %s %s IS   3A  %2d%2d %2d.%2d\n",trig_stations[ch],comp,hr1,mn1,sc1,ms1);
        break;
        case 5:
        fprintf(sf," %s%s IS   3A  %2d%2d %2d.%2d\n",trig_stations[ch],comp,hr1,mn1,sc1,ms1);
        printf(" %s%s IS   3A  %2d%2d %2d.%2d\n",trig_stations[ch],comp,hr1,mn1,sc1,ms1);
        break;

        break;
      }
    }
  }

  fclose(sf);
  for(ch=0;ch < nchannels;ch++)
  {
    for(i=0;i<256;i++)
    {
      TRIGGER_TIMES_P[ch][i]='\0';
      TRIGGER_TIMES_S[ch][i]='\0';
    }
  }
}






DC(int nsamp,float value[])
{
  int k;
  float dc;
  dc=0;
  for(k=0;k<nsamp;k++)      // add all samples in buffer
    dc = dc + value[k];
  dc = dc/(float)nsamp;    // dc of nsamp samples
  return(dc);
}

LTA_WND(int pindex,float ltas[],int nsamp)
{
  int i;
  int cod;
  float lta_normal;
//  for(i=0;i<nsamp;i++)
//    printf("%5d %8.3f\n",i,ltas[i]*1000.0);
  lta_normal= ltas[pindex+50];
  printf("LTA_NORMAL: %6.3f\n",lta_normal);
  for(i=(pindex+50);i<nsamp;i++)
  {
    printf("%5d %6.3f\n",i,ltas[i]);
    if(ltas[i] < lta_normal*0.8)
    {
      cod=i;
      return(cod);
    }
  }
  return(0);
}
CODA(int pindex,int ss,float ltas[],int nsamp)
{
  int i;
  int cod;
  float lta_normal;
//  for(i=0;i<nsamp;i++)
//    printf("%5d %8.3f\n",i,ltas[i]*1000.0);
  lta_normal= ltas[pindex-400];
//  printf("LTA_NORMAL: %8.4f %8.4f\n",lta_normal,lta_normal*1.5);
  for(i=(ss+30);i<nsamp;i++)
  {
    if(ltas[i] <= lta_normal*5.0)
    {
      cod=i;
      return(cod);
    }
  }
  return(0);
}

/** date functions */

static char daytab[2][13] = {
    {0, 31, 28, 31, 30, 31, 30, 31, 31, 30, 31, 30, 31},
    {0, 31, 29, 31, 30, 31, 30, 31, 31, 30, 31, 30, 31}
};

/** function to set month / day from day of year */

void MonthDay(int year, int yearday, int* pmonth, int* pday) {
    int i, leap;

    leap = (year % 4 == 0 && year % 100 != 0) || year % 400 == 0;
    for (i = 1; yearday > daytab[leap][i]; i++)
        yearday -= daytab[leap][i];
    *pmonth = i;
    *pday = yearday;

}
//----------------------------------------------------------------------------------
// process each component for phases
//----------------------------------------------------------------------------------
process_component(int cmp,char station[],int nsamp,int sr,char timestamp[],float value[])
{
  int i,j;
  float dc;
  float dt;
  float srate;



  char  stname[30];
  float fl;
  float fh;
  int ch;
  int cod;
  int pp;
  int ss;
  float cutoff=5.0;   // start value to register a pick
  int nenvfil=50;
  float hanfil=0.33;
  float sta=0.1;
  float lta=1.5;
  int   npks;
  float tstart = 0.0;
  float MxAmp;
    int pixA;
    int pixB;
  int idif;
  int minval;
  int maxval;
  int pn;
  int sx;
  int wind;
  int pindex;
  float diff;
  float diffs;
  int nep;
  int ws=5;                  // sta in samples
  int wl=200;                // lta in samples
  ch=0;
  fl=2.0;
  fh=8.0;
/*************************************************************************************/
/*     PICKER                                                                        */
/*************************************************************************************/
    int YEAR;
    int DOY;
    int HOUR;
    int MINU;
    int SEC;
    int MSEC;
    char dummy[200];
    int n;
    int space;
    
    BOOLEAN_INT useMemory = TRUE_INT; 

    double longTermWindow = 10.0; // NOTE: auto set below
    double threshold1 = 10.0;
    double threshold2 = 10.0;
    double tUpEvent = 0.5; // NOTE: auto set below
    double filterWindow = 4.0; // NOTE: auto set below
    //
    // auto set values
    // get dt
    double dtt = 1.0/(float)sr;
    filterWindow = 300.0 * dtt;
    long iFilterWindow = (long) (0.5 + filterWindow * 1000.0);
    if (iFilterWindow > 1)
        filterWindow = (double) iFilterWindow / 1000.0;
    //
    longTermWindow = 500.0 * dtt; // seconds
    long ilongTermWindow = (long) (0.5 + longTermWindow * 1000.0);
    if (ilongTermWindow > 1)
        longTermWindow = (double) ilongTermWindow / 1000.0;
    //
    //tUpEvent = 10.0 * dt;   // time window to take integral of charFunct version
    tUpEvent = 20.0 * dtt; // AJL20090522
    long itUpEvent = (long) (0.5 + tUpEvent * 1000.0);
    if (itUpEvent > 1)
        tUpEvent = (double) itUpEvent / 1000.0;
/*
    printf("dt            : %f\n",dtt);
    printf("filterWindow  : %f\n",filterWindow);
    printf("longTermWindow: %f\n",longTermWindow);
    printf("threshold1    : %f\n",threshold1);
    printf("threshold2    : %f\n",threshold2);
    printf("tUpEvent      : %f\n",tUpEvent);
*/    
    // do picker function test
    // definitive pick data
    PickData** pick_list_definative = NULL;
    int num_picks_definative = 0;
    // persistent memory
    FilterPicker5_Memory* mem = NULL;

    int proc_samples = 0;
    int read_samples = nsamp;
//    printf("read_samples: %4d\n",read_samples);
        // temporary data
    PickData** pick_list = NULL; // array of num_picks ptrs to PickData structures/objects containing returned picks
        int num_picks = 0;

        Pick(
                dtt,
                value,
                read_samples,
                filterWindow,
                longTermWindow,
                threshold1,
                threshold2,
                tUpEvent,
                &mem,
                useMemory,
                &pick_list,
                &num_picks,
                "TEST"
                );
//        printf("num_picks: %d\n",num_picks);    

        // save pick data
        for (n = 0; n < num_picks; n++) {
            PickData* pick = *(pick_list + n);
            pick->indices[0] += proc_samples; // pick indices returned are relative to start of packet
            pick->indices[1] += proc_samples;
            addPickToPickList(pick, &pick_list_definative, &num_picks_definative);
        }
        // clean up temporary data
        free(pick_list); // do not use free_PickList() since we want to keep PickData objects

        proc_samples += read_samples;

//printf("timestamp: %s\n",timestamp);

    int month, day;

//                2013-03-02T14:51:51.419539
sscanf(timestamp,"%4d%1c%2d%1c%2d%1c%2d%1c%2d%1c%2d%1c%3d",&YEAR,dummy,&month,dummy,&day,dummy,&HOUR,dummy,&MINU,dummy,&SEC,dummy,&MSEC);

//printf("YEAR: %d month: %d day: %d HOUR: %d MINUTE: %d SEC: %d MSEC: %d\n",YEAR,month,day,HOUR,MINU,SEC,MSEC);
//    MonthDay(YEAR, DOY, &month, &day);
//printf("month: %d day: %d\n",month,day);

//    double sec = (double) sachdr.B + (double) SEC + (double) MSEC / 1000.0;
    double sec = (double) SEC + (double) MSEC / 1000.0;
    // id fields
    char onset[] = "?";
    char* kstnm;
    kstnm = calloc(1, 16 * sizeof (char));
    sprintf(kstnm, "STAT");
    char* kinst;
    kinst = calloc(1, 16 * sizeof (char));
    sprintf(kinst, "DIG");
    if (strstr(kinst, "(count") != NULL)
        strcpy(kinst, "(counts)");
    char* kcmpnm;
    kcmpnm = calloc(1, 16 * sizeof (char));
    sprintf(kcmpnm, "CMP");
    char phase[16];
    // create NLL picks
    char* pick_str;
    pick_str = calloc(1, 1024 * sizeof (char));
    for (n = 0; n < num_picks_definative; n++) {
        sprintf(phase, "P%d_", n);
        pick_str = printNlloc(pick_str,
                *(pick_list_definative + n), dtt, kstnm, kinst, kcmpnm, onset, phase,
                YEAR, month, day, HOUR, MINU, sec);
        // write pick to <pick_file> in NLLOC_OBS format
//        fprintf(fp, "%s\n", pick_str);
//        printf("%s\n",pick_str);
        sprintf(PICKLINES[n],"%s",pick_str);
    }
if(prt >=0)
{
    for(n=0;n<num_picks_definative; n++)
    {
      printf("%2d %s PICKLINES: %s\n",n,station,PICKLINES[n]);
//      printf("%c\n",PICKLINES[n][26]);
    }
}
    for(n=0;n<200;n++)
      PICKTIMES[cmp][n] = '\0';
    sprintf(PICKTIMES[cmp],"%3d",num_picks_definative);
    if(num_picks_definative > 0)
    {
      for(n=0;n<num_picks_definative; n++)
      {
        space=n * 21 +1 + 3 + n;
        PICKTIMES[cmp][space-1] = ' ';

        if(PICKLINES[n][34]==' ')PICKLINES[n][34]='0';
        PICKTIMES[cmp][space] = PICKLINES[n][34];       // DAY
        if(PICKLINES[n][35]==' ')PICKLINES[n][35]='0';        
        PICKTIMES[cmp][space +  1] = PICKLINES[n][35];  // DAY
        PICKTIMES[cmp][space +  2] = '/';
        if(PICKLINES[n][32]==' ')PICKLINES[n][32]='0';
        PICKTIMES[cmp][space +  3] = PICKLINES[n][32];  // MONTH
        PICKTIMES[cmp][space +  4] = PICKLINES[n][33];  // MONTH
        PICKTIMES[cmp][space +  5] = '/';
        PICKTIMES[cmp][space +  6] = PICKLINES[n][30];  // YEAR
        PICKTIMES[cmp][space +  7] = PICKLINES[n][31];  // YEAR
        PICKTIMES[cmp][space +  8] = ' ';
        if(PICKLINES[n][37]==' ')PICKLINES[n][37]='0';
        PICKTIMES[cmp][space +  9] = PICKLINES[n][37];  // HOUR
        PICKTIMES[cmp][space + 10] = PICKLINES[n][38];  // HOUR
        PICKTIMES[cmp][space + 11] = ':';
        if(PICKLINES[n][39]==' ')PICKLINES[n][39]='0';
        PICKTIMES[cmp][space + 12] = PICKLINES[n][39];  // MINUTE
        PICKTIMES[cmp][space + 13] = PICKLINES[n][40];  // MINUTE
        PICKTIMES[cmp][space + 14] = ':';
        if(PICKLINES[n][44]==' ')PICKLINES[n][44]='0';
        PICKTIMES[cmp][space + 15] = PICKLINES[n][44];  // SEC.MSEC
        PICKTIMES[cmp][space + 16] = PICKLINES[n][45];  // SEC.MSEC
        PICKTIMES[cmp][space + 17] = PICKLINES[n][46];  // SEC.MSEC
        if(PICKLINES[n][47]==' ')PICKLINES[n][47]='0';
        PICKTIMES[cmp][space + 18] = PICKLINES[n][47];  // SEC.MSEC
        if(PICKLINES[n][48]==' ')PICKLINES[n][48]='0';
        PICKTIMES[cmp][space + 19] = PICKLINES[n][48];  // SEC.MSEC
        PICKTIMES[cmp][space + 20] = PICKLINES[n][49];  // SEC.MSEC
      }
//      printf("PICKTIMES: %s\n",PICKTIMES[cmp]);
    }

    // clean up
//    fclose(fp);
    free(pick_str);
    free(kcmpnm);
    free(kinst);
    free(kstnm);
    free_PickList(pick_list_definative, num_picks_definative); // PickData objects freed here
    free_FilterPicker5_Memory(&mem);	
/*************************************************************************************/
/*     End of PICKER                                                                 */
/*************************************************************************************/

//printf("%3d %s %6d %4d %s\n",cmp,station,nsamp,sr,timestamp);
  sprintf(STATION_NAMES[cmp],"%s",station);
//printf("process_component STATION_NAMES: %s\n",STATION_NAMES[cmp]);  
  sprintf(TIME_START[cmp],"%s",timestamp);
  TIME_START[cmp][10] = ' ';
  for(i=0;i<MAXSMP;i++)
  {
    fseismo[i]=0;
    env[i]=0;
    h[i]=0;
    xtap[i]=0;
    ratios[i]=0;
    smratios[i]=0;
    pk[0][i]=0;
    pk[1][i]=0;
  }  
  srates[cmp]=sr;
//  sprintf(TRIGGER_STATIONS[cmp],"%s",station);
  dc=DC(nsamp,value);                                      // find dc for raw data
  for(i=0;i<nsamp;i++)                                     // remove dc
    value[i] = (value[i] - dc)*1.0;
  for(i=0;i<nsamp;i++)
    vraw[cmp][i]=value[i];
  nusmp[cmp]=nsamp;
  for(i=0;i<nsamp;i++)          
    fseismo[i] = value[i];
  dt = 1.0/(float)(sr);                                    // samplerate in seconds

  eqavst=dt;                                               // samplerate for filter
  orden=4;                                                 // filter order
  design(orden,"BP","BU",30.0,0.3,fl,fh,eqavst);           // design filter coeficients
  for(j=0;j<(nsects*3+1);j++)
  {
    snc[cmp][j]=sn[j];
    sdc[cmp][j]=sd[j];
  }
  cset[cmp]=nsects;

  appl(cmp,nsamp);                                         // filter one channel (fsamp)
  ZERO_START();                                            // adjust start of signal after filtering
//  MaxAmp(nsamp,&MxAmp);                                    // Find maximum amplitude in fsamp
//  Scale_to_1(nsamp,MxAmp);                                 // Scale fsamp to max 1
  for(i=0;i<nsamp;i++)
  {

    vflt[cmp][i]=fseismo[i];
//    if(cmp==0 && (i>300 &&i<400))
//      printf("%2d %10.6f %10.2f %10.6f %10.6f\n",i,fseismo[i],value[i],vraw[cmp][i],vflt[cmp][i]);
  }
  STALTAPHFL(nsamp, ws, wl);                               // compute STA/LTA ratios of fseismo, store in rat
  pindex=FIND_P(nsamp);                                    // find P in ratios
//  TRIGGER_TIME_P(timestamp,ws,pindex,cmp,dt);              // detection buffer-time


// pick phase arrivals

  nep= 120*sr+pindex;
  if( nep>= nsamp)
    nep=nsamp;
  npks=MAKEPICKS(ratios,tstart,dt,nep,cutoff,pindex,cmp);

  pixA = PeakA(npks,dt,cmp);
  pixB = PeakB(pixA,npks,dt,cmp);



  if(pixA < pixB)
  {
    pp=pixA;
    ss=pixB;
  }else{
    pp=pixB;
    ss=pixA;
  }


  npks=MAKEPICKS_N(nsamp,dt,tstart,sta,lta,cutoff,nenvfil,hanfil,npks);   

  if(npks != 0)
  {
    minval=(int)(pk[0][0]/dt);
    if(minval < 100)
      minval = 200000;
    for(i=1;i<2;i++)
    {
      if((int)(pk[0][i]/dt) < minval)
	minval = (int)(pk[0][i]/dt);
    }
    pn=minval;
  }
  if(npks != 0)
  {
    maxval=(int)(pk[0][0]/dt);
    for(i=1;i<3;i++)
    {
      if((int)(pk[0][i]/dt) > maxval)
	maxval = (int)(pk[0][i]/dt);
    }
    sx=maxval;
  }

    


//  TRIGGER_TIME_S(timestamp,ws,ss,cmp,dt);                 // detection buffer-time + time to P


sprintf(stname,"%s",station);
//  printf("%s\n",stname);

  diff= (pn - pindex) * dt;
  if(diff < 0.0)             // P earlier
  {
    if( (fabs)(diff) < 0.5)                                // use pn if earlier than pindex
      pindex = pn;

  }

  diffs= (sx - ss) * dt;
  if(diffs < 0.0)             // S earlier
  {
    if( (fabs)(diff) < 0.5)
      ss = sx;

  }

  cod = CODA(pindex,ss,ltas,nsamp);                        // find coda in ltas after s-phase
  varighet[cmp] = (int)((float)cod * dt);

  TRIGGER_TIME_P(timestamp,ws,pindex,cmp,dt);              // detection buffer-time
  TRIGGER_TIME_S(timestamp,ws,ss,cmp,dt);                 // detection buffer-time + time to P


  idif=abs(pindex-pp);
if(prt >=10)
{
  if(idif <= 12)
  {
    printf("%3d*%s %s %3d %6.3f  %5d  %5d  %5d %5d %5d %5d %5d %8.3f %8.3f *\n",cmp,stname,timestamp,sr,dt,nsamp,pindex,pp,pn,ss,sx,cod,diff,diffs);
  }else{
    printf("%3d %s %s %3d %6.3f  %5d  %5d  %5d %5d %5d %5d \n",cmp,stname,timestamp,sr,dt,nsamp,pindex,pp,pn,ss,sx);
  }
}
  fseismo[pindex] = -1.0;                                  // mark p-phase in checkfiles for GNUPLOT
  fseismo[pixB] = -1.0;                                    // mark s-phase in checkfiles for GNUPLOT

  ltas[cod] = 1.0;                                         // mark coda in checkfiles for GNUPLOT




}
les_linje(FILE *streamfp,char linje[])
{
  int l;
  int br;
  br=0;
  for(l=0;l<85;l++)
  {
    linje[l]=fgetc(streamfp);
    if(linje[l] == 0xa && l == 0)
    {
      br=2;
      break;
    }
    if(linje[l] == 0xa && l > 0)
    {
      linje[l]='\0';
      break;        
    }
    if(linje[l] == EOF)
    {
//      printf("End of file. File contains: %2d items\n",i);

      br=1;
      return(br);
    }
  }
  return(br);
}
/*----------------------------------------------------------------*/
/* Read wavefile converted from miniseed to askii: trigger.ask    */
/*----------------------------------------------------------------*/
read_trigger(int comp_cnt)
{
  FILE *streamfp;
  char trg_lines[500][100];
  char linje[100];

  char station[256];
  char stacpy[256];
  char current_station[256];
  char stdstation[256];
  char sequence[15];
  char type[5];
  char block[5];
  int  ns;
  char text1[10];
  int  sr;
  int  sr_first;
  char text2[10];
  char timestamp[256];
  char timestampcpy[256];
  float value[MAXSMP];
  int nsamp;
  int tegn;
  int ret;
  int heldiv;
  int antall_linjer;
  int v1,v2,v3,v4,v5,v6;
  int first;
  char filnam1[80];
  int printarg;
  int i,l,n,no,p;
  int no_to_record;
  int br=0;
  int cmp;
  int teller;
  int len;
  char dum[100];



  p=0;
  first = 0;
  cmp=0;

for(teller=0;teller < comp_cnt;teller++)
{
//  printf("%2d %s\n",teller,COMP_HDR[teller]);
  sprintf(linje,"%s",COMP_HDR[teller]);
  sscanf(linje,"%s %s %s %d %s %d %s %s",dum,station,dum,&ns,dum,&sr,dum,timestamp);
//  printf("%2d %s %5d %5d %s\n",teller,station,ns,sr,timestamp);
  for(no=0;no<ns;no++)
    value[no]=verdier[teller][no];
  sr_first = sr;
  sprintf(current_station,"%s",station);
  sprintf(timestampcpy,"%s",timestamp);
//  printf("%2d STATION: %s %2d\n",teller,current_station,strlen(current_station));
  len=strlen(current_station);
  
//  printf("current_station: %s\n",current_station);
  
// Check station name. A standard station name of the form XX_OSL___00_EHZ is generated
// where the XX represent missing networkname
//-------------------------------------------------------------------
// 000000000011111
// 012345678901234
//
// CX_ABCDE_BHE_D
// CX_PB09__BHE_D
// CX_ABC___BHE_D
// GE_LVC_00_BHZ_D
// GE_SFJD_10_BH1_D
// GE_ABCDE_10_HHZ_D  
  if(station[0] != '_')                      // check if station name does NOT contain network name
  {
    if(station[2] == '_' && station[6] == '_' && station[12] == '_')  // CX_PB0___BHE_D 
    {
      stdstation[0]  = current_station[0];     // networkname
      stdstation[1]  = current_station[1];     // networkname
      stdstation[2]  = '_';
      stdstation[3]  = current_station[3];
      stdstation[4]  = current_station[4];
      stdstation[5]  = current_station[5];
      stdstation[6]  = '_';
      stdstation[7]  = '_';
      stdstation[8]  = '_';      
      stdstation[9]  = '0';
      stdstation[10] = '0';
      stdstation[11] = '_';
      stdstation[12] = current_station[9];
      stdstation[13] = current_station[10];
      stdstation[14] = current_station[11];
      stdstation[15] = '\0';
//      printf("%s\n",stdstation);
    }    
    if(station[2] == '_' && station[7] == '_' && station[12] == '_')  // CX_PB09__BHE_D 
    {
      stdstation[0]  = current_station[0];     // networkname
      stdstation[1]  = current_station[1];     // networkname
      stdstation[2]  = '_';
      stdstation[3]  = current_station[3];
      stdstation[4]  = current_station[4];
      stdstation[5]  = current_station[5];
      stdstation[6]  = current_station[6];
      stdstation[7]  = '_';
      stdstation[8]  = '_';      
      stdstation[9]  = '0';
      stdstation[10] = '0';
      stdstation[11] = '_';
      stdstation[12] = current_station[9];
      stdstation[13] = current_station[10];
      stdstation[14] = current_station[11];
      stdstation[15] = '\0';
//      printf("%s\n",stdstation);
    }
    if(station[2] == '_' && station[8] == '_' && station[13] == '_')  // CX_PB099__BHE_D 
    {
      stdstation[0]  = current_station[0];     // networkname
      stdstation[1]  = current_station[1];     // networkname
      stdstation[2]  = '_';
      stdstation[3]  = current_station[3];
      stdstation[4]  = current_station[4];
      stdstation[5]  = current_station[5];
      stdstation[6]  = current_station[6];
      stdstation[7]  = current_station[7];
      stdstation[8]  = '_';
      stdstation[9]  = '0';
      stdstation[10] = '0';
      stdstation[11] = '_';
      stdstation[12] = current_station[10];
      stdstation[13] = current_station[11];
      stdstation[14] = current_station[12];
      stdstation[15] = '\0';
//      printf("stdstation:  %s\n",stdstation);
    }    
    if(station[2] == '_' && station[6] == '_' && station[9] == '_') // GE_LVC_00_BHZ_D
    {
//      printf("LOC: %c%c\n",station[7],station[8]);
      stdstation[0]  = current_station[0];     // networkname
      stdstation[1]  = current_station[1];     // networkname
      stdstation[2]  = '_';
      stdstation[3]  = current_station[3];
      stdstation[4]  = current_station[4];
      stdstation[5]  = current_station[5];
      stdstation[6]  = '_';
      stdstation[7]  = '_';
      stdstation[8]  = '_';      
      stdstation[9]  = current_station[7];
      stdstation[10] = current_station[8];
      stdstation[11] = '_';
      stdstation[12] = current_station[10];
      stdstation[13] = current_station[11];
      stdstation[14] = current_station[12];
      stdstation[15] = '\0';
//      printf("%s\n",stdstation);
    }
    if(station[2] == '_' && station[7] == '_' && station[10] == '_') // GE_SFJD_10_BHZ_D
    {
//      printf("LOC: %c%c\n",station[7],station[8]);
      stdstation[0]  = current_station[0];     // networkname
      stdstation[1]  = current_station[1];     // networkname
      stdstation[2]  = '_';
      stdstation[3]  = current_station[3];
      stdstation[4]  = current_station[4];
      stdstation[5]  = current_station[5];
      stdstation[6]  = current_station[6];
      stdstation[7]  = '_';
      stdstation[8]  = '_';      
      stdstation[9]  = current_station[8];
      stdstation[10] = current_station[9];
      stdstation[11] = '_';
      stdstation[12] = current_station[11];
      stdstation[13] = current_station[12];
      stdstation[14] = current_station[13];
      stdstation[15] = '\0';
//      printf("%s\n",stdstation);
    }

/*    
    stdstation[0]  = current_station[0];     // networkname
    stdstation[1]  = current_station[1];     // networkname
    stdstation[2]  ='_';
    stdstation[3]  = current_station[3];
    stdstation[4]  = current_station[4];
    stdstation[5]  = current_station[5];
    stdstation[6]  = current_station[6];
    stdstation[7]  = current_station[7];
    stdstation[8]  ='_';
    stdstation[9]  ='0';
    stdstation[10] ='0';
    stdstation[11] ='_';
    tegn=0;
    for(i=0;i<5;i++)
    {
      if(current_station[i+3] != '_')
        tegn++;
    }
    switch(tegn)
    {
      case 3:
      stdstation[12]  = current_station[8];
      stdstation[13]  = current_station[9];
      stdstation[14]  = current_station[10];
      stdstation[15]  ='\0';
      break;
      case 4:
      stdstation[12]  = current_station[9];
      stdstation[13]  = current_station[10];
      stdstation[14]  = current_station[11];
      stdstation[15]  ='\0';
      break;
      case 5:
      stdstation[12]  = current_station[10];
      stdstation[13]  = current_station[11];
      stdstation[14]  = current_station[12];
      stdstation[15]  ='\0';
      break;
    }
*/    
    
  }else{
    stdstation[0]  ='X';
    stdstation[1]  ='X';
    stdstation[2]  ='_';
    stdstation[3]  = current_station[1];
    stdstation[4]  = current_station[2];
    stdstation[5]  = current_station[3];
    stdstation[6]  = current_station[4];
    stdstation[7]  = current_station[5];
    stdstation[8]  ='_';
    stdstation[9]  ='0';
    stdstation[10] ='0';
    stdstation[11] ='_';
    tegn=0;
    for(i=0;i<5;i++)
    {
      if(current_station[i+1] != '_')
        tegn++;
    }
    switch(tegn)
    {
      case 3:
      stdstation[12]  = current_station[6];
      stdstation[13]  = current_station[7];
      stdstation[14]  = current_station[8];
      stdstation[15]  ='\0';
      if(stdstation[14] == '_')
      {
        stdstation[14] = stdstation[13];
	stdstation[13] = '_';
      }
      break;
      case 4:
      stdstation[12]  = current_station[7];
      stdstation[13]  = current_station[8];
      stdstation[14]  = current_station[9];
      stdstation[15]  ='\0';
      if(stdstation[14] == '_')
      {
        stdstation[14] = stdstation[13];
	stdstation[13] = '_';
      }
      break;
      case 5:
      stdstation[12]  = current_station[8];
      stdstation[13]  = current_station[9];
      stdstation[14]  = current_station[10];
      stdstation[15]  ='\0';
      if(stdstation[14] == '_')
      {
        stdstation[14] = stdstation[13];
	stdstation[13] = '_';
      }
      break;
    }
  }
//--------------------------------------------------------------
// process this component for phases
//--------------------------------------------------------------
      process_component(cmp,stdstation,ns,sr_first,timestampcpy,value);
      cmp++;

}
//exit(0);





      return(cmp);
}


/*  Subroutine to apply an iir filter to a data sequence.                   */
/*    The filter is assumed to be stored as second order sections.          */
/*    Filtering is in-place.                                                */
/*    Zero-phase (forward and reverse) is an option.                        */

//appl(sw,chan,nsamp,peker)
void appl(int fnd,int nsamp)
//int      sw;
//int      chan;
//int      nsamp;
//float    *peker;
{
//float    *tpek;
    /* System generated locals */
int      i1;
int      i2;

    /* Local variables */
int      jptr;
int      i;
int      j;
float    b0;
float    b1;
float    b2;
float    a1;
float    a2;
float    x1;
float    x2;
float    y1;
float    y2;
float    output;
//tpek=peker;

/* Parameter adjustments                                                    */
/*    --data;*/

/* Function Body                                                            */
jptr = 1;
i1 = cset[fnd];

for (j = 1; j <= i1; ++j)
{
    x1 = save_x1[fnd][j];
    x2 = save_x2[fnd][j];
    y1 = save_y1[fnd][j];
    y2 = save_y2[fnd][j];
    b0 = snc[fnd][jptr];
    b1 = snc[fnd][jptr + 1];
    b2 = snc[fnd][jptr + 2];
    a1 = sdc[fnd][jptr + 1];
    a2 = sdc[fnd][jptr + 2];
    i2 = nsamp;
//    i2 = ns[fnd];
//    peker=tpek;
    for (i = 0; i < i2; ++i)
    {
//      output = b0 * (*peker) + b1 * x1 + b2 * x2;
      output = b0 * fseismo[i] + b1 * x1 + b2 * x2;

      output -= a1 * y1 + a2 * y2;
      y2 = y1;
      y1 = output;
      x2 = x1;
      x1 = fseismo[i];
      fseismo[i] = output;
//      peker++;
    }

    jptr += 3;

    save_y2[fnd][j]=y2;
    save_y1[fnd][j]=y1;
    save_x2[fnd][j]=x2;
    save_x1[fnd][j]=x1;

}

}
/*  Subroutine to design IIR digital filters from analog prototypes. */
/*  Input Arguments: */
/*  ---------------- */
/*    IORD                Filter order (10 MAXIMUM) */
/*    TYPE                Character*2 variable containing filter type */
/*                          LOWPASS (LP) */
/*                          HIGHPASS (HP) */
/*                          BANDPASS (BP) */
/*                          BANDREJECT (BR) */
/*   APROTO              Character*2 variable designating analog prototype*/
/*                          Butterworth (BU) */
/*                          Bessel (BE) */
/*                          Chebyshev Type I (C1) */
/*                          Chebyshev Type II (C2) */
/*    A                   Chebyshev stopband attenuation factor */
/*    TRBNDW              Chebyshev transition bandwidth (fraction of */
/*                          lowpass prototype passband width) */
/*    FL                  Low-frequency cutoff */
/*    FH                  High-frequency cutoff */
/*    TS                  Sampling interval (in seconds) */
/*  Output Arguments: */
/*  ----------------- */
/*    SN                  Array containing numerator coefficients of */
/*                        second-order sections packed head-to-tail. */
/*    SD                  Array containing denominator coefficients */
/*                        of second-order sections packed head-to-tail. */

void design( iord, type,aproto,  a,trbndw,  fl,  fh,  ts)
int iord;
char *type;
char *aproto;
float a;
float trbndw;
float fl;
float fh;
float ts;
{
    /* System generated locals */
    float r__1;

    /* Local variables */
    complex p[10], z[10];
    char stype[3*10];
    float omegar, ripple;
    float fhw, eps, flw, dcvalue;

    /*  Analog prototype selection */

    if (strncmp(aproto, "BU", 2) == 0) {
    buroots(p, stype, &dcvalue, iord);
    } else if (strncmp(aproto, "BE", 2) == 0) {
    beroots(p, stype, &dcvalue, iord);
    } else if (strncmp(aproto, "C1", 2) == 0) {
    chebparm(a, trbndw, iord, &eps, &ripple);
    c1roots(p, stype, &dcvalue, iord, &eps);
    } else if (strncmp(aproto, "C2", 2) == 0) {
    omegar = trbndw + (float)1.;
    c2roots(p, z, stype, &dcvalue, iord, a, omegar);
    }

    /*  Analog mapping selection */

    if (strncmp(type, "BP", 2) == 0) {
    r__1 = fl * ts / (float)2.;
    flw = warp(&r__1, &c_b12);
    r__1 = fh * ts / (float)2.;
    fhw = warp(&r__1, &c_b12);
    lptbp(p, z, stype, &dcvalue, &flw, &fhw);
    } else if (strncmp(type, "BR", 2) == 0) {
    r__1 = fl * ts / (float)2.;
    flw = warp(&r__1, &c_b12);
    r__1 = fh * ts / (float)2.;
    fhw = warp(&r__1, &c_b12);
    lptbr(p, z, stype, &dcvalue, &flw, &fhw, &sn[1], &sd[1]);
    } else if (strncmp(type, "LP", 2) == 0) {
    r__1 = fh * ts / (float)2.;
    fhw = warp(&r__1, &c_b12);
    lp(p, z, stype, &dcvalue, &sn[1], &sd[1]);
    cutoffs(&sn[1], &sd[1], &fhw);
    } else if (strncmp(type, "HP", 2) == 0) {
    r__1 = fl * ts / (float)2.;
    flw = warp(&r__1, &c_b12);
    lpthp(p, z, stype, &dcvalue, &sn[1], &sd[1]);
    cutoffs(&sn[1], &sd[1], &flw);
    }

    /*  Bilinear analog to digital transformation */
    bilin2(&sn[1], &sd[1]);

    return;
}

/* BUROOTS -- SUBROUTINE TO COMPUTE BUTTERWORTH POLES FOR */
/*   NORMALIZED LOWPASS FILTER */
/* LAST MODIFIED:  SEPTEMBER 7, 1990 */
/*  OUTPUT ARGUMENTS: */
/*  ----------------- */
/*      P              COMPLEX ARRAY CONTAINING POLES */
/*                       CONTAINS ONLY ONE FROM EACH */
/*                       COMPLEX CONJUGATE PAIR, AND */
/*                       ALL REAL POLES */
/*      RTYPE          CHARACTER ARRAY INDICATING 2ND ORDER SECTION */
/*                       TYPE: */
/*                         (SP)  SINGLE REAL POLE */
/*                         (CP)  COMPLEX CONJUGATE POLE PAIR */
/*                         (CPZ) COMPLEX CONJUGATE POLE-ZERO PAIRS */
/*      DCVALUE        MAGNITUDE OF FILTER AT ZERO FREQUENCY */
/*  INPUT ARGUMENTS: */
/*  ---------------- */
/*      IORD           DESIRED FILTER ORDER */
int buroots(p,rtype, dcvalue, iord)
complex *p;
char *rtype;
float *dcvalue;
int iord;
{
    /* System generated locals */
    int i__1, i__2;
    double d__1, d__2;
    complex q__1;

    /* Local variables */
    int half, k;
    float angle, pi;

    /* Parameter adjustments */
    rtype -= 3;
    --p;

    /* Function Body */
    pi = (float)3.14159265;

    half = iord / 2;

    /* TEST FOR ODD ORDER, AND ADD POLE AT -1 */

    nsects = 0;
    if (half << 1 < iord) {
    p[1].r = (float)-1., p[1].i = (float)0.;
    strncpy(rtype + 3, "SP", 2);
    nsects = 1;
    }
    i__1 = half;
    for (k = 1; k <= i__1; ++k) {
    angle = pi * ((float) ((k << 1) - 1)/(float) (iord << 1)+(float).5);
    ++(nsects);
    i__2 = nsects;
    d__1 = cos(angle);
    d__2 = sin(angle);
    q__1.r = d__1, q__1.i = d__2;
    p[i__2].r = q__1.r, p[i__2].i = q__1.i;
        strncpy(rtype + nsects * 3, "CP", 2);
    }
    *dcvalue = (float)1.;
    return TRUE;
}

/*  CHEBPARM - Calculates Chebyshev type I and II design parameters */
/*  INPUT ARGUMENTS */
/*  --------------- */
/*       A                Desired stopband attenuation */
/*                          i.e. max stopband amplitude is 1/ATTEN */
/*       TRBNDW           Transition bandwidth between stop and passbands */
/*                          as a fraction of the passband width */
/*       IORD             Filter order (number of poles) */
/*  OUTPUT ARGUMENTS */
/*  ---------------- */
/*       EPS              Chebyshev passband parameter */
/*       RIPPLE           Passband ripple */

static int chebparm(a, trbndw, iord, eps, ripple)
float a, trbndw;
int iord;
float *eps, *ripple;
{
    /* System generated locals */
    float r__1, r__2;

    /* Local variables */
    float g, alpha, omegar;

    omegar = trbndw + 1.;
    /* Computing 2nd power */
    r__2 = omegar;
    r__1 = omegar + sqrt(r__2 * r__2 - 1.);
    alpha = pow((double)r__1, (double)iord);
    /* Computing 2nd power */
    r__1 = alpha;
    g = (r__1 * r__1 + (float)1.) / (alpha * (float)2.);
    /* Computing 2nd power */
    r__1 = a;
    *eps = sqrt(r__1 * r__1 - (float)1.) / g;
    /* Computing 2nd power */
    r__1 = *eps;
    *ripple = (float)1. / sqrt(r__1 * r__1 + (float)1.);

    return TRUE;
}

/* C1ROOTS -- SUBROUTINE TO COMPUTE CHEBYSHEV TYPE I POLES FOR */
/*   NORMALIZED LOWPASS FILTER */
/* LAST MODIFIED:  SEPTEMBER 7, 1990 */
/*  OUTPUT ARGUMENTS: */
/*  ----------------- */
/*      P              COMPLEX ARRAY CONTAINING POLES */
/*                       CONTAINS ONLY ONE FROM EACH */
/*                       COMPLEX CONJUGATE PAIR, AND */
/*                       ALL REAL POLES */
/*      RTYPE          CHARACTER ARRAY INDICATING 2ND ORDER SECTION */
/*                       TYPE: */
/*                         (SP)  SINGLE REAL POLE */
/*                         (CP)  COMPLEX CONJUGATE POLE PAIR */
/*                         (CPZ) COMPLEX CONJUGATE POLE-ZERO PAIRS */
/*      DCVALUE        RESPONSE OF FILTER AT ZERO FREQUENCY */
/*  INPUT ARGUMENTS: */
/*  ---------------- */
/*      IORD           DESIRED FILTER ORDER */
/*      EPS            CHEBYSHEV PAR RELATED TO PASSBAND RIPPLE */

static int c1roots(p, rtype, dcvalue, iord, eps)
complex *p;
char *rtype;
float *dcvalue;
int iord;
float *eps;
{
    /* System generated locals */
    int i__1, i__2;
    float r__1;
    double d__1;
    complex q__1;

    /* Local variables */
    int half;
    float c;
    int i;
    float gamma, s, angle, omega, sigma, pi;

    /* Parameter adjustments */
    rtype -= 3;
    --p;

    /* Function Body */
    pi = (float)3.14159265;
    half = iord / 2;

    /*  INTERMEDIATE DESIGN PARS */

    gamma = (sqrt(*eps * *eps + (float)1.) + (float)1.) / *eps;
    gamma = log(gamma) / (float) (iord);
    gamma = exp(gamma);
    s = (gamma - (float)1. / gamma) * (float).5;
    c = (gamma + (float)1. / gamma) * (float).5;

    /*  CALCULATE POLES */

    nsects = 0;
    i__1 = half;
    for (i = 1; i <= i__1; ++i) {
    strncpy(rtype + i * 3, "CP", 2);
    angle = (float) ((i << 1) - 1) * pi / (float) (iord << 1);
    sigma = -(double)s * sin(angle);
    omega = c * cos(angle);
    i__2 = i;
    q__1.r = sigma, q__1.i = omega;
    p[i__2].r = q__1.r, p[i__2].i = q__1.i;
    ++(nsects);
    }
    if (half << 1 < iord) {
    strncpy(rtype + (half + 1) * 3, "SP", 2);
    i__1 = half + 1;
    d__1 = -(double)s;
    q__1.r = d__1, q__1.i = (float)0.;
    p[i__1].r = q__1.r, p[i__1].i = q__1.i;
    ++(nsects);
    *dcvalue = (float)1.;
    } else {
    /* Computing 2nd power */
    r__1 = *eps;
    *dcvalue = (float)1. / sqrt(r__1 * r__1 + 1);
    }
    return TRUE;
}

/* C2ROOTS -- SUBROUTINE TO COMPUTE ROOTS FOR NORMALIZED LOWPASS */
/*   CHEBYSHEV TYPE 2 FILTER */
/* LAST MODIFIED:  SEPTEMBER 7, 1990 */
/*  OUTPUT ARGUMENTS: */
/*  ----------------- */
/*      P              COMPLEX ARRAY CONTAINING POLES */
/*                       CONTAINS ONLY ONE FROM EACH */
/*                       COMPLEX CONJUGATE PAIR, AND */
/*                       ALL REAL POLES */
/*      Z              COMPLEX ARRAY CONTAINING ZEROS */
/*                       CONTAINS ONLY ONE FROM EACH */
/*                       COMPLEX CONJUGATE PAIR, AND */
/*                       ALL REAL ZEROS */
/*      RTYPE          CHARACTER ARRAY INDICATING 2ND ORDER SECTION */
/*                       TYPE: */
/*                         (SP)  SINGLE REAL POLE */
/*                         (CP)  COMPLEX CONJUGATE POLE PAIR */
/*                         (CPZ) COMPLEX CONJUGATE POLE-ZERO PAIRS */
/*      DCVALUE        MAGNITUDE OF FILTER AT ZERO FREQUENCY */
/*  INPUT ARGUMENTS: */
/*  ---------------- */
/*      IORD           DESIRED FILTER ORDER */
/*      A              STOPBAND ATTENUATION FACTOR */
/*      OMEGAR         CUTOFF FREQUENCY OF STOPBAND */
/*                     PASSBAND CUTOFF IS AT 1.0 HERTZ */

static int c2roots(p, z, rtype, dcvalue, iord, a, omegar)
complex *p, *z;
char *rtype;
float *dcvalue;
int iord;
float a, omegar;
{
    /* System generated locals */
    int i__1, i__2;
    double d__1;
    complex q__1;

    /* Local variables */
    int half;
    float beta, c;
    int i;
    float gamma, s, alpha, angle, omega, sigma, denom, pi;

    /* Parameter adjustments */
    rtype -= 3;
    --z;
    --p;

    /* Function Body */
    pi = (float)3.14159265;
    half = iord / 2;

/*  INTERMEDIATE DESIGN PARS */

    gamma = a + sqrt(a * a - (float)1.);
    gamma = log(gamma) / (float) (iord);
    gamma = exp(gamma);
    s = (gamma - (float)1. / gamma) * (float).5;
    c = (gamma + (float)1. / gamma) * (float).5;
    nsects = 0;
    i__1 = half;
    for (i = 1; i <= i__1; ++i) {

/*  CALCULATE POLES */

    strncpy(rtype + i * 3, "CPZ", 3);
    angle = (float) ((i << 1) - 1) * pi / (float) (iord << 1);
    alpha = -(double)s * sin(angle);
    beta = c * cos(angle);
    denom = alpha * alpha + beta * beta;
    sigma = omegar * alpha / denom;
    omega = -(double)(omegar) * beta / denom;
    i__2 = i;
    q__1.r = sigma, q__1.i = omega;
    p[i__2].r = q__1.r, p[i__2].i = q__1.i;

/*  CALCULATE ZEROS */

    omega = omegar / cos(angle);
    i__2 = i;
    q__1.r = (float)0., q__1.i = omega;
    z[i__2].r = q__1.r, z[i__2].i = q__1.i;
    ++(nsects);
    }

/*  ODD-ORDER FILTERS */

    if (half << 1 < iord) {
    strncpy(rtype + (half + 1) * 3, "SP", 2);
    i__1 = half + 1;
    d__1 = -(double)(omegar) / s;
    q__1.r = d__1, q__1.i = (float)0.;
    p[i__1].r = q__1.r, p[i__1].i = q__1.i;
    ++(nsects);
    }

/*  DC VALUE */

    *dcvalue = 1.;
    return TRUE;
}

/* Subroutine to convert an prototype lowpass filter to a bandpass filter via
*/
/*    the analog polynomial transformation.  The lowpass filter is */
/*   described in terms of its poles and zeros (as input to this routine).*/
/*    The output consists of the parameters for second order sections. */
/*  Input Arguments: */
/*  ---------------- */
/*    P                       Array containing poles */
/*    Z                       Array containing zeros */
/*    RTYPE                   Character array containing type information */
/*                              (SP) single real pole  or */
/*                              (CP) complex conjugate pole pair  or */
/*                              (CPZ) complex conjugate pole/zero pairs */
/*    DCVALUE                 Zero frequency value of filter */
/*    FL                      Low-frequency cutoff */
/*    FH                      High-frequency cutoff */
/*  Output Arguments: */
/*  ----------------- */
/*    SN                      Numerator polynomials for second order */
/*                              sections. */
/*    SD                      Denominator polynomials for second order */
/*                              sections. */
/*                              This subroutine doubles the number of */
/*                              sections. */

static int lptbp(p, z, rtype, dcvalue, fl, fh)
complex *p, *z;
char *rtype;
float *dcvalue;
float *fl, *fh;
{
    /* System generated locals */
    int i__1, i__2, i__3, i__4, i__5, i__6, i__7;
    double d__1;
    complex q__1, q__2, q__3, q__4, q__5, q__6, q__7, q__8, q__9, q__10;

    /* Local variables */
    int iptr;
    float a, b;
    complex h;
    int i, n;
    complex s;
    float scale;
    complex ctemp, p1, p2;
    float twopi;
    complex z1, z2;
    float pi;

    /* Parameter adjustments */
    rtype -= 3;
    --z;
    --p;

    /* Function Body */
    pi = (float)3.14159265;
    twopi = pi * (float)2.;
    a = twopi * twopi * *fl * *fh;
    b = twopi * (*fh - *fl);

    n = nsects;
    nsects = 0;
    iptr = 1;
    i__1 = n;
    for (i = 1; i <= i__1; ++i) {
    if (strncmp(rtype + i * 3, "CPZ", 3) == 0) {
        i__2 = i;
        q__3.r = b * z[i__2].r, q__3.i = b * z[i__2].i;
            q__2 = cpowi(q__3, 2);
        d__1 = a * (float)4.;
        q__1.r = q__2.r - d__1, q__1.i = q__2.i;
        ctemp.r = q__1.r, ctemp.i = q__1.i;
            q__1 = csqrt1(ctemp);
        ctemp.r = q__1.r, ctemp.i = q__1.i;
        i__2 = i;
        q__3.r = b * z[i__2].r, q__3.i = b * z[i__2].i;
        q__2.r = q__3.r + ctemp.r, q__2.i = q__3.i + ctemp.i;
        q__1.r = q__2.r * (float).5, q__1.i = q__2.i * (float).5;
        z1.r = q__1.r, z1.i = q__1.i;
        i__2 = i;
        q__3.r = b * z[i__2].r, q__3.i = b * z[i__2].i;
        q__2.r = q__3.r - ctemp.r, q__2.i = q__3.i - ctemp.i;
        q__1.r = q__2.r * (float).5, q__1.i = q__2.i * (float).5;
        z2.r = q__1.r, z2.i = q__1.i;
        i__2 = i;
        q__3.r = b * p[i__2].r, q__3.i = b * p[i__2].i;
            q__2 = cpowi(q__3, 2);
        d__1 = a * (float)4.;
        q__1.r = q__2.r - d__1, q__1.i = q__2.i;
        ctemp.r = q__1.r, ctemp.i = q__1.i;
            q__1 = csqrt1(ctemp);
        ctemp.r = q__1.r, ctemp.i = q__1.i;
        i__2 = i;
        q__3.r = b * p[i__2].r, q__3.i = b * p[i__2].i;
        q__2.r = q__3.r + ctemp.r, q__2.i = q__3.i + ctemp.i;
        q__1.r = q__2.r * (float).5, q__1.i = q__2.i * (float).5;
        p1.r = q__1.r, p1.i = q__1.i;
        i__2 = i;
        q__3.r = b * p[i__2].r, q__3.i = b * p[i__2].i;
        q__2.r = q__3.r - ctemp.r, q__2.i = q__3.i - ctemp.i;
        q__1.r = q__2.r * (float).5, q__1.i = q__2.i * (float).5;
        p2.r = q__1.r, p2.i = q__1.i;
            q__2 = conjg(z1);
        q__1.r = z1.r * q__2.r - z1.i * q__2.i, q__1.i = z1.r * q__2.i +
            z1.i * q__2.r;
        sn[iptr] = q__1.r;
        sn[iptr + 1] = z1.r * (float)-2.;
        sn[iptr + 2] = (float)1.;
            q__2 = conjg(p1);
        q__1.r = p1.r * q__2.r - p1.i * q__2.i, q__1.i = p1.r * q__2.i +
            p1.i * q__2.r;
        sd[iptr] = q__1.r;
        sd[iptr + 1] = p1.r * (float)-2.;
        sd[iptr + 2] = (float)1.;
        iptr += 3;
            q__2 = conjg(z2);
        q__1.r = z2.r * q__2.r - z2.i * q__2.i, q__1.i = z2.r * q__2.i +
            z2.i * q__2.r;
        sn[iptr] = q__1.r;
        sn[iptr + 1] = z2.r * (float)-2.;
        sn[iptr + 2] = (float)1.;
            q__2 = conjg(p2);
        q__1.r = p2.r * q__2.r - p2.i * q__2.i, q__1.i = p2.r * q__2.i +
            p2.i * q__2.r;
        sd[iptr] = q__1.r;
        sd[iptr + 1] = p2.r * (float)-2.;
        sd[iptr + 2] = (float)1.;
        iptr += 3;
        nsects += 2;
    } else if (strncmp(rtype + i * 3, "CP", 2) == 0) {
        i__2 = i;
        q__3.r = b * p[i__2].r, q__3.i = b * p[i__2].i;
            q__2 = cpowi(q__3, 2);
        d__1 = a * (float)4.;
        q__1.r = q__2.r - d__1, q__1.i = q__2.i;
        ctemp.r = q__1.r, ctemp.i = q__1.i;
            q__1 = csqrt1(ctemp);
        ctemp.r = q__1.r, ctemp.i = q__1.i;
        i__2 = i;
        q__3.r = b * p[i__2].r, q__3.i = b * p[i__2].i;
        q__2.r = q__3.r + ctemp.r, q__2.i = q__3.i + ctemp.i;
        q__1.r = q__2.r * (float).5, q__1.i = q__2.i * (float).5;
        p1.r = q__1.r, p1.i = q__1.i;
        i__2 = i;
        q__3.r = b * p[i__2].r, q__3.i = b * p[i__2].i;
        q__2.r = q__3.r - ctemp.r, q__2.i = q__3.i - ctemp.i;
        q__1.r = q__2.r * (float).5, q__1.i = q__2.i * (float).5;
        p2.r = q__1.r, p2.i = q__1.i;
        sn[iptr] = (float)0.;
        sn[iptr + 1] = b;
        sn[iptr + 2] = (float)0.;
            q__2 = conjg(p1);
        q__1.r = p1.r * q__2.r - p1.i * q__2.i, q__1.i = p1.r * q__2.i +
            p1.i * q__2.r;
        sd[iptr] = q__1.r;
        sd[iptr + 1] = p1.r * (float)-2.;
        sd[iptr + 2] = (float)1.;
        iptr += 3;
        sn[iptr] = (float)0.;
        sn[iptr + 1] = b;
        sn[iptr + 2] = (float)0.;
            q__2 = conjg(p2);
        q__1.r = p2.r * q__2.r - p2.i * q__2.i, q__1.i = p2.r * q__2.i +
            p2.i * q__2.r;
        sd[iptr] = q__1.r;
        sd[iptr + 1] = p2.r * (float)-2.;
        sd[iptr + 2] = (float)1.;
        iptr += 3;
        nsects += 2;
    } else if (strncmp(rtype + i * 3, "SP", 2) == 0) {
        sn[iptr] = (float)0.;
        sn[iptr + 1] = b;
        sn[iptr + 2] = (float)0.;
        sd[iptr] = a;
        i__2 = i;
        sd[iptr + 1] = -(double)b * p[i__2].r;
        sd[iptr + 2] = (float)1.;
        iptr += 3;
        ++(nsects);
    }
    }
/* Scaling - use fact that the bandpass filter amplitude at sqrt( omega_l **/
/*            equals the amplitude of the lowpass prototype at d.c. */
    d__1 = sqrt(a);
    q__1.r = (float)0., q__1.i = d__1;
    s.r = q__1.r, s.i = q__1.i;
    h.r = (float)1., h.i = (float)0.;

    iptr = 1;
    i__1 = nsects;
    for (i = 1; i <= i__1; ++i) {
    i__2 = iptr + 2;
    q__6.r = sn[i__2] * s.r, q__6.i = sn[i__2] * s.i;
    i__3 = iptr + 1;
    q__5.r = q__6.r + sn[i__3], q__5.i = q__6.i;
    q__4.r = q__5.r * s.r - q__5.i * s.i, q__4.i = q__5.r * s.i + q__5.i *
         s.r;
    i__4 = iptr;
    q__3.r = q__4.r + sn[i__4], q__3.i = q__4.i;
    q__2.r = h.r * q__3.r - h.i * q__3.i, q__2.i = h.r * q__3.i + h.i *
        q__3.r;
    i__5 = iptr + 2;
    q__10.r = sd[i__5] * s.r, q__10.i = sd[i__5] * s.i;
    i__6 = iptr + 1;
    q__9.r = q__10.r + sd[i__6], q__9.i = q__10.i;
    q__8.r = q__9.r * s.r - q__9.i * s.i, q__8.i = q__9.r * s.i + q__9.i *
         s.r;
    i__7 = iptr;
    q__7.r = q__8.r + sd[i__7], q__7.i = q__8.i;
        q__1 = cdiv(q__2, q__7);
    h.r = q__1.r, h.i = q__1.i;
    iptr += 3;
    }
    q__2.r = *dcvalue, q__2.i = (float)0.;
    d__1 = h.r;
    q__5 = conjg(h);
    q__4.r = d__1 * q__5.r, q__4.i = d__1 * q__5.i;
    q__3 = csqrt1(q__4);
    q__1 = cdiv(q__2, q__3);
    scale = q__1.r;
    sn[1] *= scale;
    sn[2] *= scale;
    sn[3] *= scale;
    return TRUE;
}

/*  Subroutine to generate second order section parameterization */
/*    from an pole-zero description for lowpass filters. */
/*  Input Arguments: */
/*  ---------------- */
/*    P                       Array containing poles */
/*    Z                       Array containing zeros */
/*   RTYPE                   Character array containing root type information
*/
/*                              (SP)  Single pole or */
/*                              (CP)  Complex conjugate pole pair */
/*                             (CPZ) Complex conjugate pole and zero pairs*/
/*    DCVALUE                 Zero-frequency value of prototype filter */
/*  Output Arguments: */
/*  ----------------- */
/*    SN                      Numerator polynomials for second order */
/*                              sections. */
/*    SD                      Denominator polynomials for second order */
/*                              sections. */

static int lp(p, z, rtype, dcvalue, sn, sd)
complex *p, *z;
char *rtype;
float *dcvalue;
float *sn, *sd;
{
    /* System generated locals */
    int i__1, i__2, i__3;
    complex q__1, q__2, q__3, q__4;

    /* Local variables */
    int iptr, i;
    float scale;

    /* Parameter adjustments */
    --sd;
    --sn;
    rtype -= 3;
    --z;
    --p;

    /* Function Body */
    iptr = 1;
    i__1 = nsects;
    for (i = 1; i <= i__1; ++i) {
    if (strncmp(rtype + i * 3, "CPZ", 3) == 0) {
        i__2 = i;
        q__2 = conjg(p[i]);
        q__1.r = p[i__2].r * q__2.r - p[i__2].i * q__2.i, q__1.i = p[i__2]
            .r * q__2.i + p[i__2].i * q__2.r;
        i__3 = i;
        q__4 = conjg(z[i]);
        q__3.r = z[i__3].r * q__4.r - z[i__3].i * q__4.i, q__3.i = z[i__3]
            .r * q__4.i + z[i__3].i * q__4.r;
        scale = q__1.r / q__3.r;
        i__2 = i;
        q__2 = conjg(z[i]);
        q__1.r = z[i__2].r * q__2.r - z[i__2].i * q__2.i, q__1.i = z[i__2]
            .r * q__2.i + z[i__2].i * q__2.r;
        sn[iptr] = q__1.r * scale;
        i__2 = i;
        sn[iptr + 1] = z[i__2].r * (float)-2. * scale;
        sn[iptr + 2] = scale * (float)1.;
        i__2 = i;
        q__2 = conjg(p[i]);
        q__1.r = p[i__2].r * q__2.r - p[i__2].i * q__2.i, q__1.i = p[i__2]
            .r * q__2.i + p[i__2].i * q__2.r;
        sd[iptr] = q__1.r;
        i__2 = i;
        sd[iptr + 1] = p[i__2].r * (float)-2.;
        sd[iptr + 2] = (float)1.;
        iptr += 3;
    } else if (strncmp(rtype + i * 3, "CP", 2) == 0) {
        i__2 = i;
        q__2 = conjg(p[i]);
        q__1.r = p[i__2].r * q__2.r - p[i__2].i * q__2.i, q__1.i = p[i__2]
            .r * q__2.i + p[i__2].i * q__2.r;
        scale = q__1.r;
        sn[iptr] = scale;
        sn[iptr + 1] = (float)0.;
        sn[iptr + 2] = (float)0.;
        i__2 = i;
        q__2 = conjg(p[i]);
        q__1.r = p[i__2].r * q__2.r - p[i__2].i * q__2.i, q__1.i = p[i__2]
            .r * q__2.i + p[i__2].i * q__2.r;
        sd[iptr] = q__1.r;
        i__2 = i;
        sd[iptr + 1] = p[i__2].r * (float)-2.;
        sd[iptr + 2] = (float)1.;
        iptr += 3;
    } else if (strncmp(rtype + i * 3, "SP", 2) == 0) {
        i__2 = i;
        scale = -(double)p[i__2].r;
        sn[iptr] = scale;
        sn[iptr + 1] = (float)0.;
        sn[iptr + 2] = (float)0.;
        i__2 = i;
        sd[iptr] = -(double)p[i__2].r;
        sd[iptr + 1] = (float)1.;
        sd[iptr + 2] = (float)0.;
        iptr += 3;
    }
    }
    sn[1] = *dcvalue * sn[1];
    sn[2] = *dcvalue * sn[2];
    sn[3] = *dcvalue * sn[3];
    return TRUE;
}

/*  Subroutine to convert a lowpass filter to a band reject filter */
/*    via an analog polynomial transformation.  The lowpass filter is */
/*   described in terms of its poles and zeros (as input to this routine).*/
/*    The output consists of the parameters for second order sections. */
/*  Input Arguments: */
/*  ---------------- */
/*    P                       Array containing poles */
/*    Z                       Array containing zeros */
/*    RTYPE                   Character array containing type information */
/*                              (SP)  single real pole or */
/*                              (CP)  complex conjugate pole pair */
/*                              (CPZ) complex conjugate pole/zero pairs */
/*    DCVALUE                 Zero-frequency value of prototype filter */
/*                              prior to transformation */
/*    FL                      Low-frequency cutoff */
/*    FH                      High-frequency cutoff */
/*  Output Arguments: */
/*  ----------------- */
/*    SN                      Numerator polynomials for second order */
/*                              sections. */
/*    SD                      Denominator polynomials for second order */
/*                              sections. */
/*    NSECTS                  Number of second order sections following */
/*                              transformation.  The number is doubled. */

static int lptbr(p, z, rtype, dcvalue, fl, fh, sn, sd)
complex *p, *z;
char *rtype;
float *dcvalue;
float *fl, *fh, *sn, *sd;
{
    /* System generated locals */
    int i__1, i__2;
    double d__1;
    complex q__1, q__2, q__3;

    /* Local variables */
    complex cinv;
    int iptr;
    float a, b, h;
    int i, n;
    float scale;
    complex ctemp, p1, p2;
    float twopi;
    complex z1, z2;
    float pi;

    /* Parameter adjustments */
    --sd;
    --sn;
    rtype -= 3;
    --z;
    --p;

    /* Function Body */
    pi = (float)3.14159265;
    twopi = pi * (float)2.;
    a = twopi * twopi * *fl * *fh;
    b = twopi * (*fh - *fl);
    n = nsects;
    nsects = 0;
    iptr = 1;
    i__1 = n;
    for (i = 1; i <= i__1; ++i) {
    if (strncmp(rtype + i * 3, "CPZ", 3) == 0) {
        q__1 = cdiv(c_b43, z[i]);
        cinv.r = q__1.r, cinv.i = q__1.i;
        q__3.r = b * cinv.r, q__3.i = b * cinv.i;
            q__2 = cpowi(q__3, 2);
        d__1 = a * (float)4.;
        q__1.r = q__2.r - d__1, q__1.i = q__2.i;
        ctemp.r = q__1.r, ctemp.i = q__1.i;
        q__1 = csqrt1(ctemp);
        ctemp.r = q__1.r, ctemp.i = q__1.i;
        q__3.r = b * cinv.r, q__3.i = b * cinv.i;
        q__2.r = q__3.r + ctemp.r, q__2.i = q__3.i + ctemp.i;
        q__1.r = q__2.r * (float).5, q__1.i = q__2.i * (float).5;
        z1.r = q__1.r, z1.i = q__1.i;
        q__3.r = b * cinv.r, q__3.i = b * cinv.i;
        q__2.r = q__3.r - ctemp.r, q__2.i = q__3.i - ctemp.i;
        q__1.r = q__2.r * (float).5, q__1.i = q__2.i * (float).5;
        z2.r = q__1.r, z2.i = q__1.i;
        q__1 = cdiv(c_b43, p[i]);
        cinv.r = q__1.r, cinv.i = q__1.i;
        q__3.r = b * cinv.r, q__3.i = b * cinv.i;
            q__2 = cpowi(q__3, 2);
        d__1 = a * (float)4.;
        q__1.r = q__2.r - d__1, q__1.i = q__2.i;
        ctemp.r = q__1.r, ctemp.i = q__1.i;
        q__1 = csqrt1(ctemp);
        ctemp.r = q__1.r, ctemp.i = q__1.i;
        q__3.r = b * cinv.r, q__3.i = b * cinv.i;
        q__2.r = q__3.r + ctemp.r, q__2.i = q__3.i + ctemp.i;
        q__1.r = q__2.r * (float).5, q__1.i = q__2.i * (float).5;
        p1.r = q__1.r, p1.i = q__1.i;
        q__3.r = b * cinv.r, q__3.i = b * cinv.i;
        q__2.r = q__3.r - ctemp.r, q__2.i = q__3.i - ctemp.i;
        q__1.r = q__2.r * (float).5, q__1.i = q__2.i * (float).5;
        p2.r = q__1.r, p2.i = q__1.i;
        q__2 = conjg(z1);
        q__1.r = z1.r * q__2.r - z1.i * q__2.i, q__1.i = z1.r * q__2.i +
            z1.i * q__2.r;
        sn[iptr] = q__1.r;
        sn[iptr + 1] = z1.r * (float)-2.;
        sn[iptr + 2] = (float)1.;
        q__2 = conjg(p1);
        q__1.r = p1.r * q__2.r - p1.i * q__2.i, q__1.i = p1.r * q__2.i +
            p1.i * q__2.r;
        sd[iptr] = q__1.r;
        sd[iptr + 1] = p1.r * (float)-2.;
        sd[iptr + 2] = (float)1.;
        iptr += 3;
        q__2 = conjg(z2);
        q__1.r = z2.r * q__2.r - z2.i * q__2.i, q__1.i = z2.r * q__2.i +
            z2.i * q__2.r;
        sn[iptr] = q__1.r;
        sn[iptr + 1] = z2.r * (float)-2.;
        sn[iptr + 2] = (float)1.;
        q__2 = conjg(p2);
        q__1.r = p2.r * q__2.r - p2.i * q__2.i, q__1.i = p2.r * q__2.i +
            p2.i * q__2.r;
        sd[iptr] = q__1.r;
        sd[iptr + 1] = p2.r * (float)-2.;
        sd[iptr + 2] = (float)1.;
        iptr += 3;
        nsects += 2;
    } else if (strncmp(rtype + i * 3, "CP", 2) == 0) {
        q__1 = cdiv(c_b43, p[i]);
        cinv.r = q__1.r, cinv.i = q__1.i;
        q__3.r = b * cinv.r, q__3.i = b * cinv.i;
            q__2 = cpowi(q__3, 2);
        d__1 = a * (float)4.;
        q__1.r = q__2.r - d__1, q__1.i = q__2.i;
        ctemp.r = q__1.r, ctemp.i = q__1.i;
        q__1 = csqrt1(ctemp);
        ctemp.r = q__1.r, ctemp.i = q__1.i;
        q__3.r = b * cinv.r, q__3.i = b * cinv.i;
        q__2.r = q__3.r + ctemp.r, q__2.i = q__3.i + ctemp.i;
        q__1.r = q__2.r * (float).5, q__1.i = q__2.i * (float).5;
        p1.r = q__1.r, p1.i = q__1.i;
        q__3.r = b * cinv.r, q__3.i = b * cinv.i;
        q__2.r = q__3.r - ctemp.r, q__2.i = q__3.i - ctemp.i;
        q__1.r = q__2.r * (float).5, q__1.i = q__2.i * (float).5;
        p2.r = q__1.r, p2.i = q__1.i;
        sn[iptr] = a;
        sn[iptr + 1] = (float)0.;
        sn[iptr + 2] = (float)1.;
        q__2 = conjg(p1);
        q__1.r = p1.r * q__2.r - p1.i * q__2.i, q__1.i = p1.r * q__2.i +
            p1.i * q__2.r;
        sd[iptr] = q__1.r;
        sd[iptr + 1] = p1.r * (float)-2.;
        sd[iptr + 2] = (float)1.;
        iptr += 3;
        sn[iptr] = a;
        sn[iptr + 1] = (float)0.;
        sn[iptr + 2] = (float)1.;
        q__2 = conjg(p2);
        q__1.r = p2.r * q__2.r - p2.i * q__2.i, q__1.i = p2.r * q__2.i +
            p2.i * q__2.r;
        sd[iptr] = q__1.r;
        sd[iptr + 1] = p2.r * (float)-2.;
        sd[iptr + 2] = (float)1.;
        iptr += 3;
        nsects += 2;
    } else if (strncmp(rtype + i * 3, "SP", 2) == 0) {
        sn[iptr] = a;
        sn[iptr + 1] = (float)0.;
        sn[iptr + 2] = (float)1.;
        i__2 = i;
        sd[iptr] = -(double)a * p[i__2].r;
        sd[iptr + 1] = b;
        i__2 = i;
        sd[iptr + 2] = -(double)p[i__2].r;
        iptr += 3;
        ++(nsects);
    }
    }
/*  Scaling - use the fact that the bandreject filter amplitude  at d.c. */
/*            equals the lowpass prototype amplitude at d.c. */
    h = (float)1.;
    iptr = 1;
    i__1 = nsects;
    for (i = 1; i <= i__1; ++i) {
    h = h * sn[iptr] / sd[iptr];
    iptr += 3;
    }
    scale = *dcvalue / ABS(h);
    sn[1] *= scale;
    sn[2] *= scale;
    sn[3] *= scale;
    return TRUE;
}

/*  Subroutine to alter the cutoff of a filter.  Assumes that the */
/*    filter is structured as second order sections.  Changes */
/*    the cutoffs of normalized lowpass and highpass filters through */
/*    a simple polynomial transformation. */
/*  Input Arguments: */
/*  ---------------- */
/*    F                       New cutoff frequency */
/*  Input/Output Arguments: */
/*  ----------------------- */
/*    SN                      Numerator polynomials for second order */
/*                              sections. */
/*    SD                      Denominator polynomials for second order */
/*                              sections. */

static int cutoffs(sn, sd, f)
float *sn, *sd;
float *f;
{
    /* System generated locals */
    int i__1;

    /* Local variables */
    int iptr, i;
    float scale;

    /* Parameter adjustments */
    --sd;
    --sn;

    /* Function Body */
    scale = *f * (float)6.2831853000000004;

    iptr = 1;
    i__1 = nsects;
    for (i = 1; i <= i__1; ++i) {

    sn[iptr + 1] /= scale;
    sn[iptr + 2] /= scale * scale;
    sd[iptr + 1] /= scale;
    sd[iptr + 2] /= scale * scale;
    iptr += 3;
    }
    return TRUE;
}

/*  Subroutine to convert a lowpass filter to a highpass filter via */
/*    an analog polynomial transformation.  The lowpass filter is */
/*   described in terms of its poles and zeroes (as input to this routine).*/
/*    The output consists of the parameters for second order sections. */
/*  Input Arguments: */
/*  ---------------- */
/*    P                       Array containing poles */
/*    Z                       Array containing zeroes */
/*   RTYPE                   Character array containing root type information
*/
/*                              (SP) single real pole or */
/*                              (CP)  complex conjugate pair */
/*                              (CPZ) complex pole/zero pairs */
/*    DCVALUE                 Zero-frequency value of prototype filter */
/*  Output Arguments: */
/*  ----------------- */
/*    SN                      Numerator polynomials for second order */
/*                              sections. */
/*    SD                      Denominator polynomials for second order */
/*                              sections. */

static int lpthp(p, z, rtype, dcvalue, sn, sd)
complex *p, *z;
char *rtype;
float *dcvalue;
float *sn, *sd;
{
    /* System generated locals */
    int i__1, i__2, i__3;
    complex q__1, q__2, q__3, q__4;

   /* Local variables */
    int iptr, i;
    float scale;

    /* Parameter adjustments */
    --sd;
    --sn;
    rtype -= 3;
    --z;
    --p;

    /* Function Body */
    iptr = 1;
    i__1 = nsects;
    for (i = 1; i <= i__1; ++i) {
    if (strncmp(rtype + i * 3, "CPZ", 3) == 0) {
        i__2 = i;
        q__2 = conjg(p[i]);
        q__1.r = p[i__2].r * q__2.r - p[i__2].i * q__2.i, q__1.i = p[i__2]
            .r * q__2.i + p[i__2].i * q__2.r;
        i__3 = i;
        q__4 = conjg(z[i]);
        q__3.r = z[i__3].r * q__4.r - z[i__3].i * q__4.i, q__3.i = z[i__3]
            .r * q__4.i + z[i__3].i * q__4.r;
        scale = q__1.r / q__3.r;
        sn[iptr] = scale * (float)1.;
        i__2 = i;
        sn[iptr + 1] = z[i__2].r * (float)-2. * scale;
        i__2 = i;
        q__2 = conjg(z[i]);
        q__1.r = z[i__2].r * q__2.r - z[i__2].i * q__2.i, q__1.i = z[i__2]
            .r * q__2.i + z[i__2].i * q__2.r;
        sn[iptr + 2] = q__1.r * scale;
        sd[iptr] = (float)1.;
        i__2 = i;
        sd[iptr + 1] = p[i__2].r * (float)-2.;
        i__2 = i;
        q__2 = conjg(p[i]);
        q__1.r = p[i__2].r * q__2.r - p[i__2].i * q__2.i, q__1.i = p[i__2]
            .r * q__2.i + p[i__2].i * q__2.r;
        sd[iptr + 2] = q__1.r;
        iptr += 3;
    } else if (strncmp(rtype + i * 3, "CP", 2) == 0) {
        i__2 = i;
        q__2 = conjg(p[i]);
        q__1.r = p[i__2].r * q__2.r - p[i__2].i * q__2.i, q__1.i = p[i__2]
            .r * q__2.i + p[i__2].i * q__2.r;
        scale = q__1.r;
        sn[iptr] = (float)0.;
        sn[iptr + 1] = (float)0.;
        sn[iptr + 2] = scale;
        sd[iptr] = (float)1.;
        i__2 = i;
        sd[iptr + 1] = p[i__2].r * (float)-2.;
        i__2 = i;
        q__2 = conjg(p[i]);
        q__1.r = p[i__2].r * q__2.r - p[i__2].i * q__2.i, q__1.i = p[i__2]
            .r * q__2.i + p[i__2].i * q__2.r;
        sd[iptr + 2] = q__1.r;
        iptr += 3;
    } else if (strncmp(rtype + i * 3, "SP", 2) == 0) {
        i__2 = i;
        scale = -(double)p[i__2].r;
        sn[iptr] = (float)0.;
        sn[iptr + 1] = scale;
        sn[iptr + 2] = (float)0.;
        sd[iptr] = (float)1.;
        i__2 = i;
        sd[iptr + 1] = -(double)p[i__2].r;
        sd[iptr + 2] = (float)0.;
        iptr += 3;
    }
    }
    sn[1] *= *dcvalue;
    sn[2] *= *dcvalue;
    sn[3] *= *dcvalue;
    return TRUE;
}

/* Transforms an analog filter to a digital filter via the bilinear transforma
ti*/
/*   Assumes both are stored as second order sections.  The transformation is
*/
/*    done in-place. */
/*  Input Arguments: */
/*  ---------------- */
/*   SN                   Array containing numerator polynomial coefficients f
or*/
/*                           second order sections.  Packed head-to-tail. */
/*   SD                   Array containing denominator polynomial coefficients
 f*/
/*                           second order sections.  Packed head-to-tail. */

static int bilin2(sn, sd)
float *sn, *sd;
{
    /* System generated locals */
    int i__1;

    /* Local variables */
    int iptr, i;
    float scale, a0, a1, a2;

    /* Parameter adjustments */
    --sd;
    --sn;

    /* Function Body */
    iptr = 1;
    i__1 = nsects;
    for (i = 1; i <= i__1; ++i) {
    a0 = sd[iptr];
    a1 = sd[iptr + 1];
    a2 = sd[iptr + 2];
    scale = a2 + a1 + a0;
    sd[iptr] = (float)1.;
    sd[iptr + 1] = (a0 - a2) * (float)2. / scale;
    sd[iptr + 2] = (a2 - a1 + a0) / scale;
    a0 = sn[iptr];
    a1 = sn[iptr + 1];
    a2 = sn[iptr + 2];
    sn[iptr] = (a2 + a1 + a0) / scale;
    sn[iptr + 1] = (a0 - a2) * (float)2. / scale;
    sn[iptr + 2] = (a2 - a1 + a0) / scale;
    iptr += 3;
    }
    return TRUE;
}

/* WARP -- FUNCTION, APPLIES TANGENT FREQUENCY WARPING TO COMPENSATE */
/*         FOR BILINEAR ANALOG -> DIGITAL TRANSFORMATION */
/* ARGUMENTS: */
/* ---------- */
/*      F       ORIGINAL DESIGN FREQUENCY SPECIFICATION (HERTZ) */
/*      TS      SAMPLING INTERVAL (SECONDS) */
/*  LAST MODIFIED:  SEPTEMBER 20, 1990 */

static double warp(f, ts)
float *f, *ts;
{
    /* System generated locals */
    float ret_val;

    /* Local variables */
    float angle, twopi;

    twopi = (float)6.2831853;
    angle = twopi * *f * *ts / (float)2.;
    ret_val = tan(angle) * (float)2. / *ts;
    ret_val /= twopi;
    return ret_val;
}

/* BEROOTS -- SUBROUTINE TO RETURN BESSEL POLES FOR */
/*   NORMALIZED LOWPASS FILTER */
/* LAST MODIFIED:  SEPTEMBER 7, 1990 */
/*  OUTPUT ARGUMENTS: */
/*  ----------------- */
/*      P              COMPLEX ARRAY CONTAINING POLES */
/*                       CONTAINS ONLY ONE FROM EACH */
/*                       COMPLEX CONJUGATE PAIR, AND */
/*                       ALL REAL POLES */
/*      RTYPE          CHARACTER ARRAY INDICATING 2ND ORDER SECTION */
/*                       TYPE: */
/*                         (SP)  SINGLE REAL POLE */
/*                         (CP)  COMPLEX CONJUGATE POLE PAIR */
/*                         (CPZ) COMPLEX CONJUGATE POLE-ZERO PAIRS */
/*      DCVALUE        MAGNITUDE OF FILTER AT ZERO FREQUENCY */
/*  INPUT ARGUMENTS: */
/*  ---------------- */
/*      IORD           DESIRED FILTER ORDER */

static int beroots(p, rtype, dcvalue, iord)
complex *p;
char *rtype;
float *dcvalue;
int iord;
{
    /* Parameter adjustments */
    rtype -= 3;
    --p;

    /* Function Body */
    if (iord == 1) {
    p[1].r = (float)-1., p[1].i = (float)0.;
    strncpy(rtype + 3, "SP", 2);
    } else if (iord == 2) {
    p[1].r = (float)-1.1016013, p[1].i = (float).6360098;
    strncpy(rtype + 3, "CP", 2);
    } else if (iord == 3) {
    p[1].r = (float)-1.0474091, p[1].i = (float).9992645;
    strncpy(rtype + 3, "CP", 2);
    p[2].r = (float)-1.3226758, p[2].i = (float)0.;
    strncpy(rtype + 6, "SP", 2);
    } else if (iord == 4) {
    p[1].r = (float)-.9952088, p[1].i = (float)1.2571058;
    strncpy(rtype + 3, "CP", 2);
    p[2].r = (float)-1.3700679, p[2].i = (float).4102497;
    strncpy(rtype + 6, "CP", 2);
    } else if (iord == 5) {
    p[1].r = (float)-.9576766, p[1].i = (float)1.4711244;
    strncpy(rtype + 3, "CP", 2);
    p[2].r = (float)-1.3808774, p[2].i = (float).7179096;
    strncpy(rtype + 6, "CP", 2);
    p[3].r = (float)-1.502316, p[3].i = (float)0.;
    strncpy(rtype + 9, "SP", 2);
    } else if (iord == 6) {
    p[1].r = (float)-.9306565, p[1].i = (float)1.6618633;
    strncpy(rtype + 3, "CP", 2);
    p[2].r = (float)-1.3818581, p[2].i = (float).9714719;
    strncpy(rtype + 6, "CP", 2);
    p[3].r = (float)-1.5714904, p[3].i = (float).3208964;
    strncpy(rtype + 9, "CP", 2);
    } else if (iord == 7) {
    p[1].r = (float)-.9098678, p[1].i = (float)1.8364514;
    strncpy(rtype + 3, "CP", 2);
    p[2].r = (float)-1.3789032, p[2].i = (float)1.1915667;
    strncpy(rtype + 6, "CP", 2);
    p[3].r = (float)-1.6120388, p[3].i = (float).5892445;
    strncpy(rtype + 9, "CP", 2);
    p[4].r = (float)-1.6843682, p[4].i = (float)0.;
    strncpy(rtype + 12, "SP", 2);
    } else if (iord == 8) {
    p[1].r = (float)-.892871, p[1].i = (float)1.9983286;
    strncpy(rtype + 3, "CP", 2);
    p[2].r = (float)-1.3738431, p[2].i = (float)1.3883585;
    strncpy(rtype + 6, "CP", 2);
    p[3].r = (float)-1.6369417, p[3].i = (float).8227968;
    strncpy(rtype + 9, "CP", 2);
    p[4].r = (float)-1.7574108, p[4].i = (float).2728679;
    strncpy(rtype + 12, "CP", 2);
    }
    nsects = iord - iord / 2;
    *dcvalue = (float)1.;
    return TRUE;
}

static complex cmul(a,b)
complex a,b;
{
    complex c;
    c.r = a.r * b.r - a.i * b.i;
    c.i = a.i * b.r + a.r * b.i;
    return c;
}

static complex cpowi(a,n)
complex a;
int n;
{
    int i;
    complex c;
    c.i = 0.;
    c.r = 1.;
    for (i = 1; i <= n; ++i) {
    c = cmul(c, a);
    }
    return c;

}
static complex csqrt1(z)
complex z;
{   complex c;
    float x,y,w,r;
    if ((z.r == 0.0) && (z.i == 0.0)) {
        c.r=0.0;
        c.i=0.0;
        return c;
    } else {
        x=fabs(z.r);
        y=fabs(z.i);
        if (x >= y) {
            r=y/x;
            w=sqrt(x)*sqrt(0.5*(1.0+sqrt(1.0+r*r)));
        } else {
            r=x/y;
            w=sqrt(y)*sqrt(0.5*(r+sqrt(1.0+r*r)));
        }
        if (z.r >= 0.0) {
            c.r=w;
            c.i=z.i/(2.0*w);
        } else {
            c.i=(z.i >= 0) ? w : -w;
            c.r=z.i/(2.0*c.i);
        }
        return c;
    }
}

static complex conjg(z)
complex z;
{   complex c;
    c.r=z.r;
    c.i = -z.i;
    return c;
}

static complex cdiv(a,b)
complex a,b;
{   complex c;
    float r,den;
    if (fabs(b.r) >= fabs(b.i)) {
        r=b.i/b.r;
        den=b.r+r*b.i;
        c.r=(a.r+r*a.i)/den;
        c.i=(a.i-r*a.r)/den;
    } else {
        r=b.r/b.i;
        den=b.i+r*b.r;
        c.r=(a.r*r+a.i)/den;
        c.i=(a.i*r-a.r)/den;
    }
    return c;
}
//*=========================================================================      
//*                                                                               
//* Routine to convert from day-of-year to date                                   
//*                                                         Leif Kvamme 12-1-85   
void DTE (DOY,DAY,MON,YR)                                          
int DOY;
int *DAY;
int *MON;
int YR;
{
int MTH[14],J,M,N;
for(J=1;J<8;J=J+2)
    MTH[J] = 31;
for(J=8;J<13;J=J+2)
    MTH[J] = 31;
MTH[2] = 28;
MTH[4] = 30;
MTH[6] = 30;
MTH[9] = 30;
MTH[11]= 30;
if ((YR%4) == 0) MTH[2] = 29;
M = 0;                                                                    
for(J=1;J<13;J++)
{
    M = M + MTH[J];
    N = DOY - M;                                                              
    if (N <= 0)
    {                                                       
        *MON = J;
        *DAY = N + MTH[J];                                                     
        goto four;                                                               
    }
}
four:                                                                   
J=0;
}
void TIMSEC (YR,MTH,DAY,HR,MIN,SECS,MSECS)                         
int YR;
int MTH;
int DAY;
int HR;
int MIN;
float SECS;
double *MSECS;                         /* total seconds to be returned      */
{
int DYR;                               /* flag for leap year                */
int IYR;                               /* number of leap-years since 1900   */
int YDY;                               /* number of days in current year    */

if(YR <= 50) YR=YR+100;                /* check for 2000 or more in 2 digits*/

if(YR >= 1900) YR=YR - 1900;
DYR = 0;
if((YR % 4) == 0) DYR = 1;
IYR = YR/4 - DYR;
//                -- Seconds to beginning of     
*MSECS = (float)(IYR*366) + (float)((YR-IYR)*365);
//                -- current year                
*MSECS = *MSECS*86400.0;
//                -- January                     
if (MTH == 1) YDY = DAY;
//                -- February                    
if (MTH == 2) YDY = DAY + 31;
//                -- ....                        
if (MTH == 3) YDY = DAY + DYR + 59;       
if (MTH == 4) YDY = DAY + DYR + 90;                                     
if (MTH == 5) YDY = DAY + DYR + 120;                                    
if (MTH == 6) YDY = DAY + DYR + 151;                                    
if (MTH == 7) YDY = DAY + DYR + 181;                                    
if (MTH == 8) YDY = DAY + DYR + 212;                                    
if (MTH == 9) YDY = DAY + DYR + 243;                                    
if (MTH ==10) YDY = DAY + DYR + 273;                                    
if (MTH ==11) YDY = DAY + DYR + 304;                                    
if (MTH ==12) YDY = DAY + DYR + 334;                                    
*MSECS = *MSECS + (float)(YDY*86400 + HR*3600 + MIN*60) + SECS;
}
SECTIM (MSECS,YR,DOY,MTH,DAY,HR,MIN,SEC)                      
double MSECS;
int *YR;
int *DOY;
int *MTH;
int *DAY;
int *HR;
int *MIN;
float *SEC;
{
int DDAY,MMTH;
//             -- Temporary seconds         
double SSEC;
//             -- Seconds per day           
double SECDAY;
//             -- Seconds per year          
double SECYR;
//              -- Leap-year indicators      
float  IND, SIN;
//              -- Counter                   
int    I;
SECDAY=86400.0;
SECYR =31536000.0;
//*-------------------------------------------------------------------------      
//*  Find year:                                                                   
//*                                                                               
SSEC = 0.0;
for(I=1;I<201;I++)
{
    IND = 0.0;
//               -- Leap year                 
    if ((I%4) == 0) IND = 1.0;
//               -- Add years                 
    SSEC = SSEC + SECYR + IND*SECDAY;
//               -- Year found                
    if (SSEC > MSECS) goto found;              
        SIN = IND;
}
found:                                                              
*YR = I - 1;
//*                                                                               
//*  Find day-of-year and date:                                                   
//*                                                                               
//c-- Reset remaining seconds   
SSEC = MSECS - (SSEC - SECYR - IND*SECDAY);

for(I=1;I<367;I++)
{
    SSEC = SSEC - SECDAY;
    if (SSEC < 0) goto incr;                                              
}
incr:

I = I - 1;
//           -- 366 days of year          
if (SIN > 0.0) I = I + 1;
//           --                           
if (I == 0)
{                           
//           -- Justifying if change      
    *YR = *YR - 1;
//           -- of year                   
    I  = 365;                                  
    if ((*YR%4) == 0) I = 366;
}


*DOY = I;
//-- convert from doy to date  
DTE (I,&DDAY,&MMTH,*YR);
*DAY=DDAY;
*MTH=MMTH;


if(*YR >= 100) *YR=*YR-100;            /* test for year 2000 or more        */





//*                                                                               
//*  Find hour:                                                                   
//*                                                                               
SSEC = SSEC + SECDAY;

for(I=1;I<25;I++)
{
    SSEC = SSEC - 60.0*60.0;
    if (SSEC < 0) goto seks;                                              
}
seks:
*HR = I - 1;
//*                                                                               
//*  Find minutes and remaining seconds:                                          
//*                                                                               
SSEC = SSEC + 60.0*60.0;
for(I=1;I<61;I++)
{
    SSEC = SSEC - 60.0;
    if (SSEC < 0) goto otte;
}
otte:                                              
   *MIN = I - 1;
   *SEC = SSEC + 60.0;
if(*YR < 50)
{
    *YR=*YR+2000;
}else{
    *YR=*YR+1900;
}
}

/***************************************************************************
 * writeascii:
 * 
 * Write data buffer to output file as ASCII.
 *
 * Returns the number of samples written or -1 on error.
 ***************************************************************************/
static int
writeascii (MSTrace *mst,int chan)
{
  char outfile[1024];
  char *outname;
  char timestr[50];
  char srcname[50];
  char *samptype;
  
  int line, col, cnt, samplesize;
  int lines;
  void *sptr;


  
  if ( ! mst )
    return -1;
  
  if ( mst->numsamples == 0 || mst->samprate == 0.0 )
    return 0;
  
  if ( verbose )
    fprintf (stderr, "Writing ASCII for %.8s.%.8s.%.8s.%.8s\n",
	     mst->network, mst->station, mst->location, mst->channel);
  
  /* Generate source name and ISO time string */
  mst_srcname (mst, srcname, 1);
  ms_hptime2isotimestr (mst->starttime, timestr, 1);
    
  /* Set sample type description */

  if ( mst->sampletype == 'i' )
    {
      samptype = "INTEGER";
    }

if(prt == -10)
{  
  /* Generate and open output file name if single file not being used */
  if ( ! ofp )
    {
      /* Create output file name: Net.Sta.Loc.Chan.Qual.Year-Month-DayTHour:Min:Sec.Subsec.txt */
      snprintf (outfile, sizeof(outfile), "%s.%s.%s.%s.%c.%s.txt",
		mst->network, mst->station, mst->location, mst->channel,
		mst->dataquality, timestr);
      
         
      
      /* Open output file */
      if ( (ofp = fopen (outfile, "wb")) == NULL )
	{
	  fprintf (stderr, "Cannot open output file: %s (%s)\n",
		   outfile, strerror(errno));
	  return -1;
	}
      
      outname = outfile;
    }
}  
  /* Header format:
   * "TIMESERIES Net_Sta_Loc_Chan_Qual, ## samples, ## sps, isotime, SLIST|TSPAIR, INTEGER|FLOAT|ASCII, Units" */
  
       /* Print header line */
if(prt == -10)
{
      fprintf (ofp, "TIMESERIES %s, %lld samples, %g sps, %s, SLIST, %s\n",
	       srcname, (long long int)mst->numsamples, mst->samprate, timestr, samptype);
}
      sprintf (COMP_HDR[chan], "TIMESERIES %s , %lld samples, %g sps, %s , SLIST, %s",
	       srcname, (long long int)mst->numsamples, mst->samprate, timestr, samptype);

//TU printf("COMP_HDR: %s\n",COMP_HDR[chan]);
      
//      lines = (mst->numsamples / slistcols) + ((slistcols == 1) ? 0 : 1);
      lines = mst->numsamples;
//      printf("lines: %d numsamples: %d\n",lines,mst->numsamples);
      
      if ( (samplesize = ms_samplesize(mst->sampletype)) == 0 )
	{
	  fprintf (stderr, "Unrecognized sample type: %c\n", mst->sampletype);
	}
      

      else
	for ( cnt = 0, line = 0; line < lines; line++ )
	  {
		if ( cnt < mst->numsamples )
		  {
		    sptr = (char*)mst->datasamples + (cnt * samplesize);
		    
		    if ( mst->sampletype == 'i' )
		      {
if(prt == -10)
			  fprintf (ofp, "%d", *(int32_t *)sptr);
//			printf("%d\n",*(int32_t *)sptr);
        verdier[chan][cnt]=*(int32_t *)sptr;
		      }

		      cnt++;
		  }
if(prt == -10)
	    fprintf (ofp, "\n");
	  }

if(prt == -10)
{  
  if ( outname == outfile )
    {
      fclose (ofp);
      ofp = 0;
    }
}  
//  fprintf (stderr, "Wrote %lld samples from %s to %s\n",
//	   (long long int)mst->numsamples, srcname, outname);
  
  return mst->numsamples;
}  /* End of writeascii() */



/***************************************************************************
 * parameter_proc:
 *
 * Process the command line parameters.
 *
 * Returns 0 on success, and -1 on failure
 ***************************************************************************/
static int
parameter_proc (int argcount, char **argvec)
{
  int optind;
  int error = 0;

  if (argcount <= 1)
    error++;

  /* Process all command line arguments */
for (optind = 1; optind < argcount; optind++)
{
  if (strcmp (argvec[optind], "-h") == 0)
  {
    usage();
    exit (0);
  }
  else if (strcmp (argvec[optind], "-sfile") == 0)
  {
    strcpy(sfilename, argvec[++optind]);
    if(prt > 0)
      printf("rtpick:sfilename..........  ....: %s\n",sfilename);	  
  }
  else if (strcmp (argvec[optind], "-wavefile") == 0)
  {
    wavefiledef = 1;
    strcpy(wfilename, argvec[++optind]);
    if(prt > 0)    
      printf("rtpick:wfilename..........  ....: %s\n",wfilename);	  
  }
  else if (strcmp (argvec[optind], "-locate") == 0)
  {
    locate = atoi(argvec[++optind]);
    if(prt > 0)   
      printf("rtpick:autolocate...............: %d\n",locate);	  
  }
  else if (strcmp (argvec[optind], "-iter") == 0)
  {
    iterations = atoi(argvec[++optind]);
    if(prt > 0)    
      printf("rtpick:iterations...............: %d\n",iterations);    
  }
  else if (strcmp (argvec[optind], "-maxres") == 0)
  {
    maxres = atof(argvec[++optind]);
    if(prt > 0)    
      printf("rtpick:maxres...................: %6.1f\n",maxres);    
  }
  else if (strcmp (argvec[optind], "-prt") == 0)
  {
    prt = atoi (argvec[++optind]);
    if(prt > 0)    
      printf("rtpick:prt......................: %d\n",prt);
  }
  else if (strcmp (argvec[optind], "-keep") == 0)
  {
    keep = atoi (argvec[++optind]);
    if(prt > 0)    
      printf("rtpick:keep.....................: %d\n",keep);
  }
  else if (strcmp (argvec[optind], "-mail") == 0)
  {
    mail1 = 1;
    strcpy(mailaddress1, argvec[++optind]);
    if(prt > 0)    
      printf("rtpick:mailaddress1........  ....: %s\n",mailaddress1);	  
  }
  else if (strncmp (argvec[optind], "-", 1 ) == 0)
  {
    fprintf(stderr, "Unknown option: %s\n", argvec[optind]);
    exit (1);
  }


//  else
//  {
      /* Add the file name to the intput file list */
//    if ( ! addnode (&filelist, NULL, 0, argvec[optind], strlen(argvec[optind])+1) )
//    {
//      fprintf (stderr, "Error adding file name to list\n");
//    }
//  }
//-------------------------------------------------
}    // for
//---------------------------------------------
// mseed2askii
//---------------------------------------------
  /* Make sure an input files were specified */
//  if ( filelist == 0 )
//    {
//      fprintf (stderr, "No input files were specified\n\n");
//      fprintf (stderr, "%s version %s\n\n", PACKAGE, VERSION);
//      fprintf (stderr, "Try %s -h for usage\n", PACKAGE);
//      exit (1);
//    }
//------------------------------------------------
  /* If errors then report the usage message and quit */
  if (error)
    {
      usage ();
      exit (1);
    }
//---------------------------------------------
// mseed2askii
//---------------------------------------------
  /* Check the input files for any list files, if any are found
   * remove them from the list and add the contained list */
  if ( filelist )
    {
      struct listnode *prevln, *ln;
      char *lfname;

      prevln = ln = filelist;
      while ( ln != 0 )
        {
          lfname = ln->data;

          if ( *lfname == '@' )
            {
              /* Remove this node from the list */
              if ( ln == filelist )
                filelist = ln->next;
              else
                prevln->next = ln->next;

              /* Skip the '@' first character */
              if ( *lfname == '@' )
                lfname++;

              /* Read list file */
//              readlistfile (lfname);

              /* Free memory for this node */
              if ( ln->key )
                free (ln->key);
              free (ln->data);
              free (ln);
            }
          else
            {
              prevln = ln;
            }
	  
          ln = ln->next;
        }
    }
//----------------------------------------------------  

  return 0;
}				/* End of parameter_proc() */

/***************************************************************************
 * getoptval:
 * Return the value to a command line option; checking that the value is 
 * itself not an option (starting with '-') and is not past the end of
 * the argument list.
 *
 * argcount: total arguments in argvec
 * argvec: argument list
 * argopt: index of option to process, value is expected to be at argopt+1
 *
 * Returns value on success and exits with error message on failure
 ***************************************************************************/
static char *
getoptval (int argcount, char **argvec, int argopt)
{
  if ( argvec == NULL || argvec[argopt] == NULL ) {
    fprintf (stderr, "getoptval(): NULL option requested\n");
    exit (1);
  }
  
  /* Special case of '-o -' usage */
  if ( (argopt+1) < argcount && strcmp (argvec[argopt], "-o") == 0 )
    if ( strcmp (argvec[argopt+1], "-") == 0 )
      return argvec[argopt+1];
  
  if ( (argopt+1) < argcount && *argvec[argopt+1] != '-' )
    return argvec[argopt+1];

  fprintf (stderr, "Option %s requires a value\n", argvec[argopt]);
  exit (1);
}  /* End of getoptval() */
/***************************************************************************
 * addnode:
 *
 * Add node to the specified list.
 *
 * Return a pointer to the added node on success and NULL on error.
 ***************************************************************************/
static struct listnode *
addnode (struct listnode **listroot, void *key, int keylen,
	 void *data, int datalen)
{
  struct listnode *lastlp, *newlp;
  
  if ( data == NULL )
    {
      fprintf (stderr, "addnode(): No data specified\n");
      return NULL;
    }
  
  lastlp = *listroot;
  while ( lastlp != 0 )
    {
      if ( lastlp->next == 0 )
        break;

      lastlp = lastlp->next;
    }

  /* Create new listnode */
  newlp = (struct listnode *) malloc (sizeof (struct listnode));
  memset (newlp, 0, sizeof (struct listnode));
  
  if ( key )
    {
      newlp->key = malloc (keylen);
      memcpy (newlp->key, key, keylen);
    }
  
  if ( data)
    {
      newlp->data = malloc (datalen);
      memcpy (newlp->data, data, datalen);
    }
  
  newlp->next = 0;
  
  if ( lastlp == 0 )
    *listroot = newlp;
  else
    lastlp->next = newlp;
  
  return newlp;
}  /* End of addnode() */
/***************************************************************************
 * usage:
 * Print the usage message and exit.
 ***************************************************************************/
static void
usage (void)
{


  fprintf (stderr, "\nUsage: %s [options]\n\n", PACKAGE);
  fprintf (stderr,
	   "General program options:\n"
	   " -h             show this usage message\n"
           "\n"
           "Options:\n"
           " -sfile    txt  name of s-file, full path\n"
           " -wavefile txt  name of waveform file, full path\n"
           " -locate   n    0 - no autolocation, 1 - autolocation (default = 0)\n"
           " -iter     n    number of iterations (default = 200)\n"
           " -maxres   n    maximum residual (default = 3.0)\n"
           " -prt      n    1-10, the higher number, the more printout (debug)\n"
           " -keep     n    1 - keep first part of s-file, 0 - write new sfile\n"
           " -mail     txt  valid email address to send event info\n"
           "\n");


}				/* End of usage() */









