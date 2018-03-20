#include <stdlib.h>
#include <string.h>
#include "WAVheader.h"

#define BLOCK_SIZE 16
#define MAX_NUM_CHANNEL 8
double sampleBuffer[MAX_NUM_CHANNEL][BLOCK_SIZE];


#define MINUS_ONE_DB 0.891251
#define MINUS_FOUR_DB 0.630957
#define MINUS_TEN_DB 0.316228
#define MINUS_THREE_POINT_NINE_DB 0.638263
#define MINUS_NINE_POINT_FIVE_DB 0.334965


double LPF11kHz[5] = {0.2751861110627305,
					  0.550372222125461,
					  0.2751861110627305,
					  -0.1652447487547689,
					  0.2659891930056907};

double HPF3kHz[5] = {0.7860324093253385,
					 -1.572064818650677,
					 0.7860324093253385,
					 -1.5098642978301007,
					 0.634265339471253};

double HPF5kHz[5] = {0.6612655734076553,
					 -1.3225311468153107,
					 0.6612655734076553, 
					 -1.1701369489124014,
					 0.47492534471822007};

unsigned int outputMode = 2;

double x_history0[] = {0, 0};	//left 
double y_history0[] = {0, 0};

double x_history1[] = {0, 0};	//right
double y_history1[] = {0, 0};

double x_history2[] = {0, 0};	//ls
double y_history2[] = {0, 0};

double x_history3[] = {0, 0};	//rs
double y_history3[] = {0, 0};

double x_history4[] = {0, 0};	//c
double y_history4[] = {0, 0};

double XTempHistory[] = {0, 0};	
double YTempHistory[] = {0, 0};

//double *tempBufferR;

double *p;


double second_order_IIR(double input, double* coefficients, double* x_history, double* y_history);
void processing();


int main(int argc, char* argv[])
{
	FILE *wav_in=NULL;
	FILE *wav_out=NULL;
	
	char WavInputName[256];
	char WavOutputName[256];
	WAV_HEADER inputWAVhdr,outputWAVhdr;	

	// Init channel buffers
	for(int i=0; i<MAX_NUM_CHANNEL; i++)
		memset(&sampleBuffer[i],0,BLOCK_SIZE);

	// Open input and output wav files
	//-------------------------------------------------
	strcpy(WavInputName, argv[1]);
	wav_in = OpenWavFileForRead (WavInputName,"rb");
	strcpy(WavOutputName, argv[2]);
	wav_out = OpenWavFileForRead (WavOutputName,"wb");
	//-------------------------------------------------
	
	// Read input wav header
	//-------------------------------------------------
	ReadWavHeader(wav_in,inputWAVhdr);
	//-------------------------------------------------
	
	// Set up output WAV header
	//-------------------------------------------------	
	outputWAVhdr = inputWAVhdr;
	outputWAVhdr.fmt.NumChannels = 5; // change number of channels

	int oneChannelSubChunk2Size = inputWAVhdr.data.SubChunk2Size/inputWAVhdr.fmt.NumChannels;
	int oneChannelByteRate = inputWAVhdr.fmt.ByteRate/inputWAVhdr.fmt.NumChannels;
	int oneChannelBlockAlign = inputWAVhdr.fmt.BlockAlign/inputWAVhdr.fmt.NumChannels;
	
	outputWAVhdr.data.SubChunk2Size = oneChannelSubChunk2Size*outputWAVhdr.fmt.NumChannels;
	outputWAVhdr.fmt.ByteRate = oneChannelByteRate*outputWAVhdr.fmt.NumChannels;
	outputWAVhdr.fmt.BlockAlign = oneChannelBlockAlign*outputWAVhdr.fmt.NumChannels;


	// Write output WAV header to file
	//-------------------------------------------------
	WriteWavHeader(wav_out,outputWAVhdr);


	// Processing loop
	//-------------------------------------------------	
	{
		int sample;
		int BytesPerSample = inputWAVhdr.fmt.BitsPerSample/8;
		const double SAMPLE_SCALE = -(double)(1 << 31);		//2^31
		int iNumSamples = inputWAVhdr.data.SubChunk2Size/(inputWAVhdr.fmt.NumChannels*inputWAVhdr.fmt.BitsPerSample/8);
		
		// exact file length should be handled correctly...
		for(int i=0; i<iNumSamples/BLOCK_SIZE; i++)
		{	
			for(int j=0; j<BLOCK_SIZE; j++)
			{
				for(int k=0; k<inputWAVhdr.fmt.NumChannels; k++)
				{	
					sample = 0; //debug
					fread(&sample, BytesPerSample, 1, wav_in);
					sample = sample << (32 - inputWAVhdr.fmt.BitsPerSample); // force signextend
					sampleBuffer[k][j] = sample / SAMPLE_SCALE;				// scale sample to 1.0/-1.0 range		
				}
			}
			//***********************************************************
			//***********************************************************
			processing();
			

			for(int j=0; j<BLOCK_SIZE; j++)
			{
				for(int k=0; k<outputWAVhdr.fmt.NumChannels; k++)
				{	
					sample = sampleBuffer[k][j] * SAMPLE_SCALE ;	// crude, non-rounding 			
					sample = sample >> (32 - inputWAVhdr.fmt.BitsPerSample);
					fwrite(&sample, outputWAVhdr.fmt.BitsPerSample/8, 1, wav_out);		
				}
			}		
		}
	}
	
	// Close files
	//-------------------------------------------------	
	//compare(wav_out, model0);

	fclose(wav_in);
	fclose(wav_out);
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

double second_order_IIR(double input, double* coefficients, double* x_history, double* y_history) 
{
    double output = 0;

    output += *(coefficients++) * input;        /* A0 * x(n)     */
    output += *(coefficients++) * *(x_history++); /* A1 * x(n-1) */
    output += *(coefficients++) * *(x_history--); /* A2 * x(n-2)   */
    output -= *(coefficients++) * *(y_history++); /* B1 * y(n-1) */
    output -= *coefficients * *(y_history--); /* B2 * y(n-2)   */

	double tmp_x = *x_history;
	double tmp_y = *y_history;

	*(++y_history) = tmp_y;	/* y(n-2) = y(n-1) */
	*(--y_history) = output;	/* y(n-1) = y(n)   */
	*(++x_history) = tmp_x;	/* x(n-2) = x(n-1) */
	*(--x_history) = input;	/* x(n-1) = x(n)   */

    return output;
}

void processing()
{
	for(p = sampleBuffer[0]; p < sampleBuffer[0] + BLOCK_SIZE; p++)
	{
		//***********************************************************
		//L
		//***********************************************************
		*p = *p * MINUS_ONE_DB;
		
		*(p + 32) = second_order_IIR(*p, LPF11kHz, x_history2, y_history2);
	
		*(p + 64) = second_order_IIR(*p, HPF5kHz, x_history4, y_history4);
		*(p + 64) = *(p + 64) * MINUS_TEN_DB;
	
		*p = second_order_IIR(*(p + 32), HPF3kHz, x_history0, y_history0);
		*p = *p * MINUS_FOUR_DB;
	
		*p = *p + *(p + 64);
	
		//delete Ls and C
		if(outputMode != 2)
		{
			for(p = sampleBuffer[0]; p < sampleBuffer[0] + BLOCK_SIZE; p++)
			{
				*(p + 64) = 0;
			}
			if(outputMode != 1)
			{
				for(p = sampleBuffer[0]; p < sampleBuffer[0] + BLOCK_SIZE; p++)
				{
					*(p + 32) = 0;
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
		if(outputMode == 0)
		{
			for(p = sampleBuffer[0]; p < sampleBuffer[0] + BLOCK_SIZE; p++)
			{
				*(p + 48) = 0;
			}
		}
	}
}