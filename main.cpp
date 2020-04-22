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
 * Compile with:		g++ -o fft_calc main.c
 * Usage: 				./fft_calc signals/sin01
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
#include "fft.c"

int main(int argc, char **argv) {

	#define SIZE_OF_LINE 		50
	#define NUM_OF_SAMPLES 		65000
	#define SIGNAL_BUF_SIZE		65536	//must be power of two and more than NUM_OF_SAMPLES
	#define WINDOW_SIZE			4096	//value used here must be a power of two
	#define SAMPLING_RATE 		10000

	FILE *fp, *fp2;

	float signal[SIGNAL_BUF_SIZE * 2];	//needs to be twice as big because complex numbers
	float fft_window[WINDOW_SIZE * 2];	//needs to be twice as big because complex numbers
	float fft[SIGNAL_BUF_SIZE / 2];		//temporary storage of results
	char fnimi2[80];
	int i,k,m,num_windows;
	char line[SIZE_OF_LINE];
	float fq_step_full, fq_step_win;

	*argv++;//don't read our own name
	argc--;//don't read our own name

	while(argc--)
	{
		if(*argv == NULL)break;
		fp=fopen(*argv, "r");
		if(!fp){printf("Failed to open %s!\n", *argv);continue;}
		else printf("Reading file %s\n", *argv);

		i=0;
		while(fgets(line, SIZE_OF_LINE, fp) != NULL)
		{
			signal[i++] = atof(line);
		}
		if(ferror(fp))perror("Error while reading file");
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
		num_windows = ceil((float)i/WINDOW_SIZE);
		for(k=0;k<num_windows-1;k++)
		{
			//need to copy each window to work buffer, because fft performs operation in place and orignal signal is lost
			for(m=0;m<WINDOW_SIZE;m++)fft_window[m] = signal[ m + WINDOW_SIZE * k ];
			//perform fft in window
			fft_complex(fft_window, WINDOW_SIZE);
			//store results until we write them to file
			for(m=0;m<WINDOW_SIZE/2;m++)fft[m + (WINDOW_SIZE / 2) * k] = fft_window[m];
		}

		//do last window
		//the last window will be smaller than 4096, so we pad zeros to the end, to get full window
		for(m=0;m<WINDOW_SIZE;m++)
		{
			if(m < NUM_OF_SAMPLES - WINDOW_SIZE * (num_windows - 1))fft_window[m] = signal[ m + WINDOW_SIZE * k ];
			else fft_window[m] = 0.0;
		}
		//perform fft in window
		fft_complex(fft_window, WINDOW_SIZE);
		//store results until we write them to file
		for(m=0;m<WINDOW_SIZE/2;m++)fft[m + WINDOW_SIZE / 2 * k] = fft_window[m];

		//calc frequency step for window fft result
		fq_step_win = ((float)SAMPLING_RATE / 2) / (WINDOW_SIZE / 2);

		//calc fft for entire signal
		//first pad zeros to the end
		for(k=65000;k<SIGNAL_BUF_SIZE;k++)signal[k] = 0.0;
		fft_complex(signal, SIGNAL_BUF_SIZE);

		//calc frequency step for full signal fft result
		fq_step_full = ((float)SAMPLING_RATE / 2) / (SIGNAL_BUF_SIZE / 2);
		
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
		strcpy(fnimi2, *argv);
		strcpy(strstr(fnimi2,".txt"), "_fft.txt\0");
		fp2=fopen(fnimi2, "w");
		if(!fp2){printf("Failed to open %s!\n", fnimi2);continue;}
		else printf("Writing to results file %s\n", fnimi2);

		fprintf(fp2,"xaxis_fq\tfull_sig\txaxis_fq\t");
		fprintf(fp2,"win1\twin2\twin3\twin4\twin5\t");
		fprintf(fp2,"win6\twin7\twin8\twin9\twin10\t");
		fprintf(fp2,"win11\twin12\twin13\twin14\twin15\twin16\n");

		for(k=0;k<WINDOW_SIZE/2;k++)
		{
			fprintf(fp2, "%lf\t%lf\t%lf\t", (fq_step_full*k), signal[k], (fq_step_win*k));
			for(m=0;m<num_windows-1;m++)fprintf(fp2, "%lf\t", fft[k + (WINDOW_SIZE/2) * m]);
			fprintf(fp2, "%lf\n", fft[k + (WINDOW_SIZE/2) * m]);
		}
		for(;k<SIGNAL_BUF_SIZE/2;k++)fprintf(fp2, "%lf\t%lf\n", fq_step_full*k, signal[k]);
		fclose(fp2);

		/*****************************************************************************
		 * Next file to read.
		 ******************************************************************************/
		*argv++;
	}
	
	printf("ENDED OK!\n");
	return 0;
}
