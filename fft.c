/**
 * Algorithm from http://paulbourke.net/miscellaneous/dft/ (2020 April 22)
 * Copyright Paul Bourke
 * Modifications by Johannes Ehala, Prolab, Taltech
 */

#include "math.h"

int calc_power(int num_samples)
{
	int p;
	int power_of_two = 0;
	if(num_samples != 0 && !(num_samples & (num_samples -1)))
	{
		p = num_samples;
		while(((p & 1) == 0) && p > 1)
		{
			p >>= 1;
			power_of_two++;
		}
	}
	return power_of_two;
}

void bit_alignment(float *data_startpoint, float *y_startpoint, int data_points) {
	int i,j,k,i2;
	double tx;

	i2 = data_points >> 1;
	j = 0;
	for (i=0;i<data_points-1;i++)
	{
		if (i < j)
		{
			tx = *(data_startpoint + i);
			*(data_startpoint + i) = *(data_startpoint + j);
			*(data_startpoint + j) = tx;
		}
		k = i2;
		while (k <= j)
		{
			j -= k;
			k >>= 1;
		}
		j += k;
	*(y_startpoint + i) = 0.0;
	}
}

void fft_complex(float *signal, int num_samples) 
{
	int i,i1,j,l,l1,l2;
	float c1,c2,t1,t2,u1,u2,z;
	float* y_startpoint;
	int p2;

	p2 = calc_power(num_samples);
	y_startpoint = signal + num_samples;

	bit_alignment(signal, y_startpoint, num_samples);

	if(p2 > 0)
	{
		c1 = -1.0;
		c2 = 0.0;
		l2 = 1;
		for (l=0;l<p2;l++)
		{
			l1 = l2;
			l2 <<= 1;
			u1 = 1.0;
			u2 = 0.0;
			for (j=0;j<l1;j++)
			{
				for (i=j;i<num_samples;i+=l2)
				{
					i1 = i + l1;
					t1 = u1 * *(signal + i1) - u2 * y_startpoint[i1];
					t2 = u1 * y_startpoint[i1] + u2 * *(signal + i1);
					*(signal + i1) = *(signal + i) - t1;
					y_startpoint[i1] = y_startpoint[i] - t2;
					*(signal + i) += t1;
					y_startpoint[i] += t2;
				}
				z =  u1 * c1 - u2 * c2;
				u2 = u1 * c2 + u2 * c1;
				u1 = z;
			}
			c2 = -(sqrt((1.0 - c1) / 2.0));
			c1 = sqrt((1.0 + c1) / 2.0);
		}

		for (i=0;i<num_samples;i++)
		{
			 *(signal + i) /= (double)num_samples;
			 *(y_startpoint + i) /= (double)num_samples;
			 *(signal + i) = sqrt(*(signal + i)* *(signal + i) + *(y_startpoint + i)* *(y_startpoint + i));
		}
	}
}
