#include <stdio.h>
#include <dsplib\wavefile.h>
#include <string.h>
#include <stdfix.h>
#include "stdfix_emu.h"
#include "common.h"

#define CHANNEL_NUM 5

#define MINUS_ONE_DB FRACT_NUM(0.891251)
#define MINUS_FOUR_DB FRACT_NUM(0.630957)
#define MINUS_TEN_DB FRACT_NUM(0.316228)
#define MINUS_THREE_POINT_NINE_DB FRACT_NUM(0.638263)
#define MINUS_NINE_POINT_FIVE_DB FRACT_NUM(0.334965)

#define ZERO FRACT_NUM(0.0)

__memX DSPfract sampleBuffer[MAX_NUM_CHANNEL][BLOCK_SIZE];

__memX DSPfract LPF11kHz[6] = {FRACT_NUM(0.2751861110627305),
					    FRACT_NUM(0.2751861110627305),
					    FRACT_NUM(0.2751861110627305),
					    FRACT_NUM(-0.08262237437738445),
					    FRACT_NUM(0.2659891930056907)};

__memX DSPfract HPF3kHz[6] = {FRACT_NUM(0.7860324093253385),
					   FRACT_NUM(-0.7860324093253385),
					   FRACT_NUM(0.7860324093253385),
					   FRACT_NUM(-0.75493214891505035),
					   FRACT_NUM(0.634265339471253)};

__memX DSPfract HPF5kHz[6] = {FRACT_NUM(0.6612655734076553),
					   FRACT_NUM(-0.66126557340765535),
					   FRACT_NUM(0.6612655734076553),
					   FRACT_NUM(-0.5850684744562007),
					   FRACT_NUM(0.47492534471822007)};

DSPint outputMode = 2;

__memX DSPfract x_history0[] = {FRACT_NUM(0), FRACT_NUM(0)};	//left
__memY DSPfract y_history0[] = {FRACT_NUM(0), FRACT_NUM(0)};

__memX DSPfract x_history1[] = {FRACT_NUM(0), FRACT_NUM(0)};	//right
__memY DSPfract y_history1[] = {FRACT_NUM(0), FRACT_NUM(0)};

__memX DSPfract x_history2[] = {FRACT_NUM(0), FRACT_NUM(0)};	//ls
__memY DSPfract y_history2[] = {FRACT_NUM(0), FRACT_NUM(0)};

__memX DSPfract x_history3[] = {FRACT_NUM(0), FRACT_NUM(0)};	//rs
__memY DSPfract y_history3[] = {FRACT_NUM(0), FRACT_NUM(0)};

__memX DSPfract x_history4[] = {FRACT_NUM(0), FRACT_NUM(0)};	//c
__memY DSPfract y_history4[] = {FRACT_NUM(0), FRACT_NUM(0)};

__memX DSPfract XTempHistory[] = {FRACT_NUM(0), FRACT_NUM(0)};
__memY DSPfract YTempHistory[] = {FRACT_NUM(0), FRACT_NUM(0)};


__memX DSPfract *p;


DSPfract second_order_IIR(DSPfract input, __memX const DSPfract* coefficients, __memX DSPfract* x_history, __memY DSPfract* y_history);
void processing();

 
int main(int argc, char *argv[])
 {
    WAVREAD_HANDLE *wav_in;
    WAVWRITE_HANDLE *wav_out;
	
	char WavInputName[256];
	char WavOutputName[256];
	
    int nChannels;
	int bitsPerSample;
    int sampleRate;
    int iNumSamples;
    int i;

	// Init channel buffers
	for(i=0; i<MAX_NUM_CHANNEL; i++)
		memset(&sampleBuffer[i],0,BLOCK_SIZE);

	// Open input wav file
	//-------------------------------------------------
	strcpy(WavInputName,argv[0]);
	wav_in = cl_wavread_open(WavInputName);
	 if(wav_in == NULL)
    {
        printf("Error: Could not open wavefile.\n");
        return -1;
    }

	//-------------------------------------------------
	
	// Read input wav header
	//-------------------------------------------------
	nChannels = cl_wavread_getnchannels(wav_in);
    bitsPerSample = cl_wavread_bits_per_sample(wav_in);
    sampleRate = cl_wavread_frame_rate(wav_in);
    iNumSamples =  cl_wavread_number_of_frames(wav_in);
	//-------------------------------------------------
	
	// Open output wav file
	//-------------------------------------------------
	strcpy(WavOutputName,argv[1]);
	wav_out = cl_wavwrite_open(WavOutputName, bitsPerSample, CHANNEL_NUM, sampleRate);
	if(!wav_out)
    {
        printf("Error: Could not open wavefile.\n");
        return -1;
    }
	//-------------------------------------------------
	
	// Processing loop
	//-------------------------------------------------	
    {
		int i;
		int j;
		int k;
		int sample;
		
		// exact file length should be handled correctly...
		for(i=0; i<iNumSamples/BLOCK_SIZE; i++)
		{	
			for(j=0; j<BLOCK_SIZE; j++)
			{
				for(k=0; k<nChannels; k++)
				{	
					sample = cl_wavread_recvsample(wav_in);
        			sampleBuffer[k][j] = rbits(sample);
				}
			}

			processing();

			for(j=0; j<BLOCK_SIZE; j++)
			{
				for(k=0; k<CHANNEL_NUM; k++)
				{	
					sample = bitsr(sampleBuffer[k][j]);
					cl_wavwrite_sendsample(wav_out, sample);
				}
			}		
		}
	}
	
	// Close files
	//-------------------------------------------------	
    cl_wavread_close(wav_in);
    cl_wavwrite_close(wav_out);
	//-------------------------------------------------	

    return 0;
 }

/**************************************
 * IIR filtar drugog reda
 *
 * input - ulazni odbirak
 * coefficients - koeficijenti [a0 a1 a2 b0 b2 b2]
 * z_x - memorija za ulazne odbirke (niz duzine 2)
 * z_y - memorija za izlazne odbirke (niz duzine 2)
 *
 * povratna vrednost - izlazni odbirak
 *
 *************************************/

DSPfract second_order_IIR(DSPfract input, __memX const DSPfract* coefficients, __memX DSPfract* x_history, __memY DSPfract* y_history)
{
	DSPaccum output = 0;

	DSPfract tmp_x;
	DSPfract tmp_y;

    //output += *(coefficients++) * input;        /* A0 * x(n)     */
	asm("x0 = xmem[%1]; %1+=1 \n\t %0+=x0*%2" : "+a"(output) , "+i"(coefficients) : "y"(input));
    //output += *(coefficients) * *x_history; /* A1 * x(n-1) */
	asm("x0 = xmem[%1] \n\t x1 = xmem[%2] \n\t %0+=x0*x1" : "+a"(output) : "i"(coefficients) , "i"(x_history));
	//output += *(coefficients++) * *(x_history++); /* A1 * x(n-1) */
	asm("x0 = xmem[%1]; %1+=1 \n\t x1 = xmem[%2]; %2+=1 \n\t %0+=x0*x1" : "+a"(output) , "+i"(coefficients) , "+i"(x_history));
    //output += *(coefficients++) * *(x_history--); /* A2 * x(n-2)   */
	asm("x0 = xmem[%1]; %1+=1 \n\t x1 = xmem[%2]; %2-=1 \n\t %0+=x0*x1" : "+a"(output) , "+i"(coefficients) , "+i"(x_history));
    //output -= *(coefficients) * *y_history; /* B1 * y(n-1) */
	asm("x0 = xmem[%1] \n\t y0 = ymem[%2] \n\t %0-=x0*y0" : "+a"(output) : "i"(coefficients), "i"(y_history));
	//output -= *(coefficients++)* *(y_history++); /* B1 * y(n-1) */
	asm("x0 = xmem[%1]; %1+=1 \n\t y0 = ymem[%2]; %2+=1 \n\t %0-=x0*y0" : "+a"(output), "+i"(coefficients), "+i"(y_history));
    output -= *(coefficients) * *(y_history--); /* B2 * y(n-2)   */
	//asm("x0 = xmem[%2] \n\t y0 = ymem[%1]; %1-=1 \n\t %0-=x0*y0" : "+a"(output), "+i"(y_history) : "i"(coefficients));

	tmp_x = *x_history;
	tmp_y = *y_history;

    *(++y_history) = tmp_y;
	*(--y_history) = output;
	*(++x_history) = tmp_x;
	*(--x_history) = input;

    return output;
}

void processing()
{
	DSPint i;
	p = sampleBuffer[0];

	for(i = 0; i < BLOCK_SIZE; p++, i++)
	{
		//***********************************************************
		//L
		//***********************************************************
		*p = *p * MINUS_ONE_DB;

		//Ls
		*(p + 32) = second_order_IIR(*p, LPF11kHz, x_history2, y_history2);

		//C
		*(p + 64) = second_order_IIR(*p, HPF5kHz, x_history4, y_history4);
		*(p + 64) = *(p + 64) * MINUS_TEN_DB;

		//
		*p = second_order_IIR(*(p + 32), HPF3kHz, x_history0, y_history0);
		*p = *p * MINUS_FOUR_DB;

		*p = *p + *(p + 64);

		//delete Ls and C
		if(outputMode - 2 < 0)
		{
			for(p = sampleBuffer[0]; p < sampleBuffer[0] + BLOCK_SIZE; p++)
			{
				*(p + 64) = ZERO;	//C
			}
			if(outputMode - 1 < 0)
			{
				for(p = sampleBuffer[0]; p < sampleBuffer[0] + BLOCK_SIZE; p++)
				{
					*(p + 32) = ZERO;	//Ls
				}
			}
		}

		//***********************************************************
		//R
		//***********************************************************

		*(p + 16) = *(p + 16) * MINUS_ONE_DB;

		*(p + 48) = second_order_IIR(*(p + 16), LPF11kHz, x_history3, y_history3);

		*(p + 80) = second_order_IIR(*(p + 16), HPF5kHz, XTempHistory, YTempHistory);
		*(p + 80) = *(p + 80) * MINUS_NINE_POINT_FIVE_DB;

		*(p + 16) = second_order_IIR(*(p + 48), HPF3kHz, x_history1, y_history1);
		*(p + 16) = *(p + 16) * MINUS_THREE_POINT_NINE_DB;

		*(p + 16) = *(p + 80) + *(p + 16);

		//delete Rs
		if(outputMode - 1 < 0)
		{
			for(p = sampleBuffer[0]; p < sampleBuffer[0] + BLOCK_SIZE; p++)
			{
				*(p + 48) = ZERO;
			}
		}
	}
}
