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
 * Compile with:        g++ -o fft_calc main.c
 * Usage:               ./fft_calc signals/sin01
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
    #define NUM_OF_SAMPLES      65000   // number of samples in test signals
    #define SIGNAL_BUF_SIZE     65536   // must be power of two and greater or equal to NUM_OF_SAMPLES
    #define MOVING_WINDOW_SIZE  4096    // value used here must be a power of two
    #define SAMPLING_RATE       10000   // Hz

    FILE *fp, *fp2;

    float signal[SIGNAL_BUF_SIZE * 2];  // needs to be twice as big because calculations use complex numbers
    float fft_window[MOVING_WINDOW_SIZE * 2];  // needs to be twice as big because calculations use complex numbers
    float fft[SIGNAL_BUF_SIZE / 2];     // temporary storage of results
    char res_file_name[80];
    uint32_t i, k, m, num_windows;
    char line[SIZE_OF_LINE];
    float fq_step_full, fq_step_win;

    argc--; //don't read our own name

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

            /*****************************************************************************
             * Calculate FFT of signal first using a ~0.4 second sliding window (kind
             * of like a spectogram) and then on the whole signal aswell. 
             * 
             * The FFT algorithm does its calculations in place (directly on the buffer
             * given to it) so the original signal is lost. Thats why we need to 
             * copy it to a temporary buffer 'fft_window'.
             * 
             * We also need to copy intermediat results to a temporary buffer 'fft'
             * before we write the results to file in the very end.
             * 
             * We also need to pad zeros to those buffers that have less samples
             * than a power two value (ie the last sliding window and full signal)
            ******************************************************************************/
            
            num_windows = ceil((float) i / MOVING_WINDOW_SIZE);
            for (k = 0; k < num_windows - 1; k++)
            {
                //copy window to work buffer
                for (m = 0; m < MOVING_WINDOW_SIZE; m++)
                {
                    fft_window[m] = signal[m + MOVING_WINDOW_SIZE * k];
                }
                
                //perform fft on window
                fft_complex(fft_window, MOVING_WINDOW_SIZE);
                
                //store results until we write them to file
                for (m = 0; m < MOVING_WINDOW_SIZE / 2; m++)
                {
                    fft[m + (MOVING_WINDOW_SIZE / 2) * k] = fft_window[m];
                }
            }

            //do last window
            //last window might be smaller than 4096, pad zeros to the end, to get full window
            for (m = 0; m < MOVING_WINDOW_SIZE; m++)
            {
                if (m < NUM_OF_SAMPLES - MOVING_WINDOW_SIZE * (num_windows - 1))
                {
                    fft_window[m] = signal[m + MOVING_WINDOW_SIZE * k];
                }
                else fft_window[m] = 0.0;
            }
            
            //perform fft on window
            fft_complex(fft_window, MOVING_WINDOW_SIZE);
            
            //store results until we write them to file
            for (m = 0; m < MOVING_WINDOW_SIZE / 2; m++)
            {
                fft[m + MOVING_WINDOW_SIZE / 2 * k] = fft_window[m];
            }

            // Calc frequency step for window FFT result.
            fq_step_win = ((float) SAMPLING_RATE / 2) / (MOVING_WINDOW_SIZE / 2);

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

                // Print header row.
                fprintf(fp2, "xaxis_fq\tfull_sig\txaxis_fq\t");
                
                m = calc_power(SIGNAL_BUF_SIZE);
                for(k = 1; k <= m; k++)
                {
                    if (k == m)
                    {
                        fprintf(fp2, "win%u\n", k);
                    }
                    else
                    {
                        fprintf(fp2, "win%u\t", k);
                    }
                }

                // Print results to file. Full signal and each window result in separate column.
                for (k = 0; k < MOVING_WINDOW_SIZE / 2; k++)
                {
                    fprintf(fp2, "%lf\t%lf\t%lf\t", (fq_step_full * k), signal[k], (fq_step_win * k));
                    
                    for ( m = 0; m < num_windows - 1; m++)
                    {
                        fprintf(fp2, "%lf\t", fft[k + (MOVING_WINDOW_SIZE / 2) * m]);
                    }
                    
                    fprintf(fp2, "%lf\n", fft[k + (MOVING_WINDOW_SIZE / 2) * m]);
                }
                
                // Now only remainder of full signal FFT.
                for ( ; k < SIGNAL_BUF_SIZE / 2; k++)
                {
                    fprintf(fp2, "%lf\t%lf\n", fq_step_full*k, signal[k]);
                }

                fclose(fp2);
            }
            
            strcpy(strstr(res_file_name, ".txt"), "_only_windows.txt\0");

            fp2 = fopen(res_file_name, "w");
            
            if (!fp2)
            {
                printf("Failed to open %s!\n", res_file_name);
            }
            else 
            {
                printf("Writing results to file %s\n", res_file_name);

                // Print header row.
                fprintf(fp2, "xaxis_fq\t");
                
                m = calc_power(SIGNAL_BUF_SIZE);
                for(k = 1; k <= m; k++)
                {
                    if (k == m)
                    {
                        fprintf(fp2, "win%u\n", k);
                    }
                    else
                    {
                        fprintf(fp2, "win%u\t", k);
                    }
                }

                // Print results to file. Full signal and each window result in separate column.
                for (k = 0; k < MOVING_WINDOW_SIZE / 2; k++)
                {
                    fprintf(fp2, "%lf\t", (fq_step_win * k));
                    
                    for ( m = 0; m < num_windows - 1; m++)
                    {
                        fprintf(fp2, "%lf\t", fft[k + (MOVING_WINDOW_SIZE / 2) * m]);
                    }
                    
                    fprintf(fp2, "%lf\n", fft[k + (MOVING_WINDOW_SIZE / 2) * m]);
                }
            
                fclose(fp2);
            }          
            

        }
    }
	
    printf("ENDED OK!\n");
    return 0;
}
