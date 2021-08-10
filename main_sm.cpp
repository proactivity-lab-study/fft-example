/**
 * Example of usage of fast Fourier Transform algorithm.
 * 
 * FFT algorithm from http://paulbourke.net/miscellaneous/dft/ (2020 April 22)
 * The algorithm expects input signal length value to be a power of two.
 * 
 * Test signals are in directory signals. There are four signals with length 6.5s
 * (sin01, sin02, sin03, random) and two signals with length 4s (sin800, sin2367).
 * 
 * Sampling rate for all test signals is 10 kHz. All signals are normalized (-1..1)
 * except sin03 which is -1.1..1.1. 
 * sin01 - pure sine wave 800 Hz 4s long from beginning, then 2.5s of silence
 * sin02 - pure sine wave 2367 Hz 4s long at the end, beginning 2.5s is silence
 * sin03 - sum of signals sin01, sin02, random
 * random - 6.5s of white noise (or close to wn)
 * sin800 - 4s sine wave 800 Hz
 * sin2367 - 4s sine wave 2367 Hz
 * 
 * Compile with:        g++ -o fft_calc main_sm.cpp fft.c -DQUIET
 * Usage:               ./fft_calc signals/smenete_1.txt
 * 
 * Some #define values need to be changed when reading 4s signal files (sin800, 
 * sin2367).
 * 
 *
 * 
 * Author Johannes Ehala, Prolab, Taltech
 * Date 22 April 2020
 *
 */
#include <string.h>
#include <cstdio>
#include <cstdlib>
#include <inttypes.h>
#include <math.h>
#include "fft.h"

int main(int argc, char *argv[]) {

    #define SIZE_OF_LINE        50      // maximum expected size of any line in test signal files
   // #define NUM_OF_SAMPLES      25000   // number of samples in test signals
    #define SIGNAL_BUF_SIZE     4096   // must be power of two and greater or equal to NUM_OF_SAMPLES
    #define MOVING_WINDOW_SIZE  2048    // value used here must be a power of two
    #define SAMPLING_RATE       3000   // Hz

    FILE *fp, *fp2;

    float signal[SIGNAL_BUF_SIZE * 2];  // needs to be twice as big because calculations use complex numbers
    float fft[SIGNAL_BUF_SIZE / 2];     // temporary storage of results
    char res_file_name[80];
    uint32_t i, k, m, num_windows;
    char line[SIZE_OF_LINE];
    float fq_step_full, fq_step_win;

    argc--; //don't read our own name
    *argv++; //don't read samples size

    int NUM_OF_SAMPLES = atoi(argv[0]);
    printf("%i\n", NUM_OF_SAMPLES);

    while(argc--)
    {
        *argv++;
        if (*argv == NULL)break;
        
        fp = fopen(*argv, "r");
        if (!fp)
        {
            printf("Failed to open %s!\n", *argv);
        }
        else 
        {
            printf("Reading file %s\n", *argv);

            i = 0;
            
            while (fgets(line, SIZE_OF_LINE, fp) != NULL)
            {
                signal[i++] = atof(line); 
            }
            
            if (ferror(fp)) perror("Error: ");
            
            printf("Read %u lines\n", i); // entire signal read into memory
            fclose(fp);

            // Normalize around 0
            // https://stats.stackexchange.com/questions/178626/how-to-normalize-data-between-1-and-1
            
            // calc min & max

            /*float min = 0;
            float max= 4096;
            float den = 0;
            for(int i = 0; i<NUM_OF_SAMPLES; i++){
                if(signal[i] < min){
                    min = signal[i];
                }else if (signal[i] > max){
                    max = signal[i];
                }
            }
            den = max-min;
            */
            //for(int i = 0; i<NUM_OF_SAMPLES; i++){
            //    signal[i] = (signal[i]-2048) / 2048;
            //}

            // Calculate FFT for entire signal.
            // First pad zeros to the end.
            for (k = NUM_OF_SAMPLES; k < SIGNAL_BUF_SIZE; k++)
            {
                signal[k] = 0.0;
            }
            fft_complex(signal, SIGNAL_BUF_SIZE);

            // Calc frequency step for full signal FFT result.
            fq_step_full = ((float) SAMPLING_RATE / 2) / (SIGNAL_BUF_SIZE / 2);

            /*****************************************************************************
             * Write results to file. File name is the same as the input file name
             * appended with '_fft'. 
             * 
             * Results are put into columns. First column is full signal FFT x-axis 
             * values and second column is the full signal fft result. 
             * 
             * Third column is the x-axis values for all sliding window results. The 
             * rest of the columns are sliding window results in temporal order.
            ******************************************************************************/
            
            strcpy(res_file_name, *argv);
            strcpy(strstr(res_file_name, ".txt"), "_fft.txt\0");
            
            fp2 = fopen(res_file_name, "w");
            
            if (!fp2)
            {
                printf("Failed to open %s!\n", res_file_name);
            }
            else 
            {
                printf("Writing results to file %s\n", res_file_name);

                // Print results to file. Full signal and each window result in separate column.
                for (k = 0; k < SIGNAL_BUF_SIZE / 2; k++)
                {
                    fprintf(fp2, "%lf\n", signal[k]);

                } 

                fclose(fp2);
            }
            

        }
    }
	
    printf("ENDED OK!\n");
    return 0;
}

// TODO: at the moment, skipping normalization bc it's done in matlab
// TODO: perform windowing here too
