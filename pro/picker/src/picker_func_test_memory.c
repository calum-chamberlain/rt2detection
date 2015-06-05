#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <math.h>

//#include "sachdr.h"

#include "ew_bridge.h"
#include "PickData.h"
#include "FilterPicker5_Memory.h"
#include "FilterPicker5.h"





#define DEBUG 1

//#define PACKET_SIZE 10 //simulated packet size (seconds)
  float value[90000];
  char timestampcpy[25];
  int  sr_first;
  char current_station[25];
  int nsamp;













void MonthDay(int year, int yearday, int* pmonth, int* pday);
int les_linje(FILE *streamfp,char linje[])
{
  int l;
  int br;
  br=0;
  for(l=0;l<80;l++)
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
/* Read recording file                                            */
/*----------------------------------------------------------------*/
int read_trigger(  char streamfile[])
{
  FILE *streamfp;
//  char trg_lines[500][100];
  char linje[100];

  char station[25];
//  char stacpy[25];
//  char current_station[25];
  char sequence[15];
  char type[5];
  char block[5];
  int  ns;
  char text1[10];
  int  sr;
//  int  sr_first;
  char text2[10];
  char timestamp[25];

//  float value[90000];
//  int nsamp;
  int tegn=0;
  int ret;
  int heldiv;
  int antall_linjer;
  int v1,v2,v3,v4,v5,v6;
  int first;
  char filnam1[80];
  int printarg=1;
  int i,n,p;
//  int no_to_record;
  int br=0;
  int cmp;
//  char dum[10];

/* Open trigger-file */
  sprintf(filnam1,"%s\0",streamfile);
  
  printf("Recording: %s\n",filnam1);
  if ((streamfp = fopen (filnam1, "rb")) == NULL)
  {
    if(printarg >= 1)
    {
      printf("Can't open recording_stream file: %s\n",filnam1);
      printf("Check if file exists, or is wrong.\n");
    }
    exit(0);
  }

  p=0;
  first = 0;
  cmp=0;


  br = les_linje(streamfp,linje);  // read first line in new station, and save samplerate + timestamp
  printf("%s\n",linje);
  sscanf(linje,"%s %s %s %s %d %s %d %s %s",station,sequence,type,block,&ns,text1,&sr,text2,timestamp);
  sr_first = sr;
  sprintf(current_station,"%s\0",station);
//printf("current_station: %s\n",current_station);
  sprintf(timestampcpy,"%s\0",timestamp);
  for(i=0;i<20;i++)                // count characters in station name
  {
    if(station[i]=='\0')
    {
      tegn=i;
      break;
    }
  }
  for(n=0;n<2000;n++)   // continue until EOF
  {
    sscanf(linje,"%s %s %s %s %d %s %d %s %s",station,sequence,type,block,&ns,text1,&sr,text2,timestamp);
    heldiv=ns - (ns/6)*6;  // rest values
    antall_linjer = ns/6;  // number of lines to read
//    printf("antall linjer: %d  rest values: %d\n",antall_linjer,heldiv);
    for(i=0;i<antall_linjer;i++)
    {
      br = les_linje(streamfp,linje);  // read line of samples
      sscanf(linje,"%d %d %d %d %d %d",&v1,&v2,&v3,&v4,&v5,&v6);
      value[p] = v1;
      p++;
      value[p] = v2;
      p++;
      value[p] = v3;
      p++;
      value[p] = v4;
      p++;
      value[p] = v5;
      p++;
      value[p] = v6;
      p++;
    }
    if(heldiv != 0)
    {
      br = les_linje(streamfp,linje);
      switch(heldiv)
      {
        case 1:
        sscanf(linje,"%d",&v1);
        value[p] = v1;
        p++;
        break;
        case 2:
        sscanf(linje,"%d %d",&v1,&v2);
        value[p] = v1;
        p++;
        value[p] = v2;
        p++;
        break;
        case 3:
        sscanf(linje,"%d %d %d",&v1,&v2,&v3);
        value[p] = v1;
        p++;
        value[p] = v2;
        p++;
        value[p] = v3;
        p++;
        break;
        case 4:
        sscanf(linje,"%d %d %d %d",&v1,&v2,&v3,&v4);
        value[p] = v1;
        p++;
        value[p] = v2;
        p++;
        value[p] = v3;
        p++;
        value[p] = v4;
        p++;
        break;
        case 5:
        sscanf(linje,"%d %d %d %d %d",&v1,&v2,&v3,&v4,&v5);
        value[p] = v1;
        p++;
        value[p] = v2;
        p++;
        value[p] = v3;
        p++;
        value[p] = v4;
        p++;
        value[p] = v5;
        p++;
        break;
      } // switch
    }   // if heldiv
    br = les_linje(streamfp,linje);
    if(br == 1)  // EOF
      break;     // end read loop
    if(br == 2)  // CR first column
    {
  //    printf("dobbel LF\n");
      br = les_linje(streamfp,linje);
      if(br == 1)  // EOF
        break;     // end read loop
    }
//    printf("%s\n",linje);
    sscanf(linje,"%s %s %s %s %d %s %d %s %s",station,sequence,type,block,&ns,text1,&sr,text2,timestamp);
    ret=strncmp(current_station,station,tegn);
    if(ret != 0)  // NEW STATION ?
    {
      nsamp=p;
//      printf("NSAMP: %d\n",nsamp);
printf("cmp: %2d current_station: %s nsamp: %5d sr_first: %4d %s\n",cmp,current_station,nsamp,sr_first,timestampcpy);
return(cmp);
//      process_component(cmp,current_station,nsamp,sr_first,timestampcpy,value);
      cmp++;
      p=0;
//      printf("new station, %s\n",station);
//      printf("%s  ",linje);
      sscanf(linje,"%s %s %s %s %d %s %d %s %s",station,sequence,type,block,&ns,text1,&sr,text2,timestamp);
      sr_first = sr;
      sprintf(current_station,"%s\n",station);
      sprintf(timestampcpy,"%s\n",timestamp);
      for(i=0;i<20;i++)                // count characters in station name
      {
        if(station[i]=='\0')
        {
          tegn=i;
          break;
        }
      }
    }
//    br = les_linje(streamfp,linje);
//    if(br == 1)  // EOF
//      break;     // end read loop
}  // -> 2000
      nsamp=p;
//      process_component(cmp,current_station,nsamp,sr,timestampcpy,value);
      cmp++;
      nsamp=p;
      return(cmp);
}



int main(int argc, char *argv[]) {

    int n,i;

    if (argc < 3) {
        printf("Usage: %s <SAC_file> <pick_file>\n", argv[0]);
        printf("  Picks are appended to end of <pick_file> in NLLOC_OBS format. \n");
        return 0;
    }
    FILE *fp;
    char filnam[200];
    char *triggerfile    = 0;
    int nchannels; 
    int length;
    int YEAR;
    int DOY;
    int HOUR;
    int MINU;
    int SEC;
    int MSEC;
    float sr;
    char dummy[100];


    BOOLEAN_INT useMemory = TRUE_INT; // set to TRUE_INT (=1) if function is called for packets of data in sequence, FALSE_INT (=0) otherwise


    sprintf(filnam,"/home/seismo/req/trigger0.ask\0");
    printf("name: %s\n",filnam);
    triggerfile = filnam;
    nchannels=read_trigger(filnam);


    double longTermWindow = 10.0; // NOTE: auto set below
    double threshold1 = 10.0;
    double threshold2 = 10.0;
    double tUpEvent = 0.5; // NOTE: auto set below
    double filterWindow = 4.0; // NOTE: auto set below
    //
    // auto set values
    // get dt
    double dt = 1.0/(float)sr_first;
    sr=1.0/(float)sr_first;
//    if (DEBUG) printf("dt: %f\n", dt);
    //dt = dt < 0.02 ? 0.02 : dt;     // aviod too-small values for high sample rate data
    //
    filterWindow = 300.0 * dt;
    long iFilterWindow = (long) (0.5 + filterWindow * 1000.0);
    if (iFilterWindow > 1)
        filterWindow = (double) iFilterWindow / 1000.0;
    //
    longTermWindow = 500.0 * dt; // seconds
    long ilongTermWindow = (long) (0.5 + longTermWindow * 1000.0);
    if (ilongTermWindow > 1)
        longTermWindow = (double) ilongTermWindow / 1000.0;
    //
    //tUpEvent = 10.0 * dt;   // time window to take integral of charFunct version
    tUpEvent = 20.0 * dt; // AJL20090522
    long itUpEvent = (long) (0.5 + tUpEvent * 1000.0);
    if (itUpEvent > 1)
        tUpEvent = (double) itUpEvent / 1000.0;

    printf("dt            : %f\n",dt);
    printf("filterWindow  : %f\n",filterWindow);
    printf("longTermWindow: %f\n",longTermWindow);
    printf("threshold1    : %f\n",threshold1);
    printf("threshold2    : %f\n",threshold2);
    printf("tUpEvent      : %f\n",tUpEvent);

    // do picker function test
    // definitive pick data
    PickData** pick_list_definative = NULL;
    int num_picks_definative = 0;
    // persistent memory
    FilterPicker5_Memory* mem = NULL;

    int proc_samples = 0;
    int read_samples = nsamp;

    printf("read_samples: %4d\n",read_samples);

        // temporary data
    PickData** pick_list = NULL; // array of num_picks ptrs to PickData structures/objects containing returned picks
        int num_picks = 0;

        Pick(
                dt,
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
        printf("num_picks: %d\n",num_picks);

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




    // create NLLOC_OBS picks
    // open pick file
    if ((fp = fopen(argv[2], "a")) == 0) {
        perror(argv[2]);
        return -1;
    }
    // date
printf("timestampcpy: %s\n",timestampcpy);

    int month, day;


sscanf(timestampcpy,"%4d%1c%3d%1c%2d%1c%2d%1c%2d%1c%3d",&YEAR,dummy,&DOY,dummy,&HOUR,dummy,&MINU,dummy,&SEC,dummy,&MSEC);

printf("YEAR: %d DOY: %d HOUR: %d MINUTE: %d SEC: %d MSEC: %d\n",YEAR,DOY,HOUR,MINU,SEC,MSEC);
    MonthDay(YEAR, DOY, &month, &day);
printf("month: %d day: %d\n",month,day);

//    double sec = (double) sachdr.B + (double) SEC + (double) MSEC / 1000.0;
    double sec = (double) SEC + (double) MSEC / 1000.0;
    // id fields
    char onset[] = "?";
    char* kstnm;
    kstnm = calloc(1, 16 * sizeof (char));
    sprintf(kstnm, "PB01\0");
    char* kinst;
    kinst = calloc(1, 16 * sizeof (char));
    sprintf(kinst, "EAM\0");
    if (strstr(kinst, "(count") != NULL)
        strcpy(kinst, "(counts)");
    char* kcmpnm;
    kcmpnm = calloc(1, 16 * sizeof (char));
    sprintf(kcmpnm, "BHZ\0");
    char phase[16];
    // create NLL picks
    char* pick_str;
    pick_str = calloc(1, 1024 * sizeof (char));
    for (n = 0; n < num_picks_definative; n++) {
        sprintf(phase, "P%d_", n);
        pick_str = printNlloc(pick_str,
                *(pick_list_definative + n), dt, kstnm, kinst, kcmpnm, onset, phase,
                YEAR, month, day, HOUR, MINU, sec);
        // write pick to <pick_file> in NLLOC_OBS format
        fprintf(fp, "%s\n", pick_str);
    }


    // clean up
    fclose(fp);
    free(pick_str);
    free(kcmpnm);
    free(kinst);
    free(kstnm);
    free_PickList(pick_list_definative, num_picks_definative); // PickData objects freed here
    free_FilterPicker5_Memory(&mem);
//    free(sample);

    return (0);

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




