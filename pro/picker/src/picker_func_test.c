#include <stdlib.h>
#include <stdio.h>
#include <string.h>

#include "sachdr.h"

#include "ew_bridge.h"
#include "PickData.h"
#include "FilterPicker5_Memory.h"
#include "FilterPicker5.h"

#define DEBUG 0


void MonthDay(int year, int yearday, int* pmonth, int* pday);

int main(int argc, char *argv[]) {


    if (argc < 3) {
        printf("Usage: %s <SAC_file> <pick_file>\n", argv[0]);
        printf("  Picks are appended to end of <pick_file> in NLLOC_OBS format. \n");
        return 0;
    }



    BOOLEAN_INT useMemory = FALSE_INT; // set to TRUE_INT (=1) if function is called for packets of data in sequence, FALSE_INT (=0) otherwise



    // open and read SAC file
    FILE *fp;
    if ((fp = fopen(argv[1], "r")) == 0) {
        perror(argv[1]);
        return -1;
    }
    // read header
    struct HDR sachdr;
///    fread(&sachdr, sizeof (sachdr), 1, fp);
    // allocate array for data
    sachdr.NPTS = 329;
    printf("sachdr.NPTS: %d\n", sachdr.NPTS);
    float* sample = calloc(sachdr.NPTS, sizeof (float));
    // read data
    int k1;
    int seq;
    float value;
    for(k1=0;k1<sachdr.NPTS;k1++)
    {
      fscanf(fp,"%d %f",&seq,&value);
      printf("%3d %8.3f\n",seq,value);
      sample[k1]=value;
//      fscanf(fp,"%d %f",seq,&sample[k1]);
    }
//    fread(sample, sizeof (float), sachdr.NPTS, fp);
    fclose(fp);
printf("%8.3f\n",sample[325]);


    // set picker paramters (TODO: make program arguments?)
    // SEE: _DOC_ in FilterPicker5.c for more details on the picker parameters
    // defaults
    // filp_test filtw 4.0 ltw 10.0 thres1 8.0 thres2 8.0 tupevt 0.2 res PICKS...
    double filterWindow = 4.0; // NOTE: auto set below
    double longTermWindow = 10.0; // NOTE: auto set below
    double threshold1 = 10.0;
    double threshold2 = 10.0;
    double tUpEvent = 0.5; // NOTE: auto set below
    //
    // auto set values
    // get dt
    sachdr.DELTA = 0.020;
    double dt = sachdr.DELTA;
    printf("sachdr.DELTA: %f\n", sachdr.DELTA);
    //dt = dt < 0.02 ? 0.02 : dt;     // aviod too-small values for high sample rate data
    //
    filterWindow = 300.0 * dt;
    printf("filterWindow: %8.3f\n",filterWindow);
    long iFilterWindow = (long) (0.5 + filterWindow * 1000.0);
    printf("iFilterWindow: %8.3f\n",iFilterWindow);
    if (iFilterWindow > 1)
        filterWindow = (double) iFilterWindow / 1000.0;
    printf("iFilterWindow: %8.3f\n",iFilterWindow);
    //
    longTermWindow = 500.0 * dt; // seconds
    long ilongTermWindow = (long) (0.5 + longTermWindow * 1000.0);
    if (ilongTermWindow > 1)
        longTermWindow = (double) ilongTermWindow / 1000.0;
    //
    tUpEvent = 20.0 * dt; // time window to take integral of charFunct version
    long itUpEvent = (long) (0.5 + tUpEvent * 1000.0);
    if (itUpEvent > 1)
        tUpEvent = (double) itUpEvent / 1000.0;
    //
    printf("picker_func_test: filp_test filtw %f ltw %f thres1 %f thres2 %f tupevt %f res PICKS\n",
            filterWindow, longTermWindow, threshold1, threshold2, tUpEvent);



    // do picker function test
    PickData** pick_list = NULL; // array of num_picks ptrs to PickData structures/objects containing returned picks
    int num_picks = 0;
    FilterPicker5_Memory* mem = NULL;

    Pick(
            sachdr.DELTA,
            sample,
            sachdr.NPTS,
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
    printf("picker_func_test: num_picks: %d\n", num_picks);

    // create NLLOC_OBS picks
    // open pick file
    if ((fp = fopen(argv[2], "a")) == 0) {
        perror(argv[2]);
        return -1;
    }
    // date
    int month, day;
    MonthDay(sachdr.NZYEAR, sachdr.NZJDAY, &month, &day);
    double sec = (double) sachdr.B + (double) sachdr.NZSEC + (double) sachdr.NZMSEC / 1000.0;
    // id fields
    char onset[] = "?";
    char* kstnm;
    kstnm = calloc(1, 16 * sizeof (char));
    strncpy(kstnm, sachdr.KSTNM, 6);
    char* kinst;
    kinst = calloc(1, 16 * sizeof (char));
    strncpy(kinst, sachdr.KINST, 6);
    if (strstr(kinst, "(count") != NULL)
        strcpy(kinst, "(counts)");
    char* kcmpnm;
    kcmpnm = calloc(1, 16 * sizeof (char));
    strncpy(kcmpnm, sachdr.KCMPNM, 6);
    char phase[16];
    // create NLL picks
    char* pick_str;
    pick_str = calloc(1, 1024 * sizeof (char));
    int n;
    for (n = 0; n < num_picks; n++) {
        sprintf(phase, "P%d_", n);
        pick_str = printNlloc(pick_str,
                *(pick_list + n), sachdr.DELTA, kstnm, kinst, kcmpnm, onset, phase,
                sachdr.NZYEAR, month, day, sachdr.NZHOUR, sachdr.NZMIN, sec);
        // write pick to <pick_file> in NLLOC_OBS format
        fprintf(fp, "%s\n", pick_str);
    }


    // clean up
    fclose(fp);
    free(pick_str);
    free(kcmpnm);
    free(kinst);
    free(kstnm);
    free_FilterPicker5_Memory(&mem);
    free_PickList(pick_list, num_picks);
    free(sample);

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

