#include <stdlib.h>
#include <string.h>
#include "WAVheader.h"
#include "stdfix_emu.h"
#include "fixed_point_math.h"
#include "common.h"

#define BLOCK_SIZE 16
#define MAX_NUM_CHANNEL 8


DSPfract sampleBuffer[MAX_NUM_CHANNEL][BLOCK_SIZE];


#define MINUS_ONE_DB FRACT_NUM(0.891251)
#define MINUS_FOUR_DB FRACT_NUM(0.630957)
#define MINUS_TEN_DB FRACT_NUM(0.316228)
#define MINUS_THREE_POINT_NINE_DB FRACT_NUM(0.638263)
#define MINUS_NINE_POINT_FIVE_DB FRACT_NUM(0.334965)

#define ZERO FRACT_NUM(0.0)


DSPfract LPF11kHz[6] = {FRACT_NUM(0.2751861110627305),
					    FRACT_NUM(0.550372222125461/2),
					    FRACT_NUM(0.2751861110627305),
					    FRACT_NUM(-0.1652447487547689/2),
					    FRACT_NUM(0.2659891930056907)};

DSPfract HPF3kHz[6] = {FRACT_NUM(0.7860324093253385),
					   FRACT_NUM(-1.572064818650677/2),
					   FRACT_NUM(0.7860324093253385),
					   FRACT_NUM(-1.5098642978301007/2),
					   FRACT_NUM(0.634265339471253)};

DSPfract HPF5kHz[6] = {FRACT_NUM(0.6612655734076553),
					   FRACT_NUM(-1.3225311468153107/2),
					   FRACT_NUM(0.6612655734076553), 
					   FRACT_NUM(-1.1701369489124014/2),
					   FRACT_NUM(0.47492534471822007)};

DSPint outputMode = 2;

DSPfract x_history0[] = {FRACT_NUM(0), FRACT_NUM(0)};	//left 
DSPfract y_history0[] = {FRACT_NUM(0), FRACT_NUM(0)};

DSPfract x_history1[] = {FRACT_NUM(0), FRACT_NUM(0)};	//right
DSPfract y_history1[] = {FRACT_NUM(0), FRACT_NUM(0)};

DSPfract x_history2[] = {FRACT_NUM(0), FRACT_NUM(0)};	//ls
DSPfract y_history2[] = {FRACT_NUM(0), FRACT_NUM(0)};

DSPfract x_history3[] = {FRACT_NUM(0), FRACT_NUM(0)};	//rs
DSPfract y_history3[] = {FRACT_NUM(0), FRACT_NUM(0)};

DSPfract x_history4[] = {FRACT_NUM(0), FRACT_NUM(0)};	//c
DSPfract y_history4[] = {FRACT_NUM(0), FRACT_NUM(0)};

DSPfract XTempHistory[] = {FRACT_NUM(0), FRACT_NUM(0)};	
DSPfract YTempHistory[] = {FRACT_NUM(0), FRACT_NUM(0)};

//double *tempBufferR;

DSPfract *p;


DSPfract second_order_IIR(DSPfract input, const DSPfract* coefficients, DSPfract* x_history, DSPfract* y_history);
void processing();


int main(int argc, char* argv[])
{
	FILE *wav_in=NULL;
	FILE *wav_out=NULL;
	char WavInputName[256];
	char WavOutputName[256];
	WAV_HEADER inputWAVhdr,outputWAVhdr;	

	// Init channel buffers
	for(DSPint i=0; i<MAX_NUM_CHANNEL; i++)
		for(DSPint j=0; j<BLOCK_SIZE; j++)
			sampleBuffer[i][j] = FRACT_NUM(0.0);

	// Open input and output wav files
	//-------------------------------------------------
	strcpy(WavInputName, argv[1]);
	wav_in = OpenWavFileForRead (WavInputName, "rb");
	strcpy(WavOutputName, argv[2]);
	wav_out = OpenWavFileForRead (WavOutputName, "wb");
	//-------------------------------------------------
	
	// Read input wav header
	//-------------------------------------------------
	ReadWavHeader(wav_in,inputWAVhdr);
	//-------------------------------------------------
	
	// Set up output WAV header
	//-------------------------------------------------	
	outputWAVhdr = inputWAVhdr;
	outputWAVhdr.fmt.NumChannels = 5; // change number of channels

	DSPint oneChannelSubChunk2Size = inputWAVhdr.data.SubChunk2Size/inputWAVhdr.fmt.NumChannels;
	DSPint oneChannelByteRate = inputWAVhdr.fmt.ByteRate/inputWAVhdr.fmt.NumChannels;
	DSPint oneChannelBlockAlign = inputWAVhdr.fmt.BlockAlign/inputWAVhdr.fmt.NumChannels;
	
	outputWAVhdr.data.SubChunk2Size = oneChannelSubChunk2Size*outputWAVhdr.fmt.NumChannels;
	outputWAVhdr.fmt.ByteRate = oneChannelByteRate*outputWAVhdr.fmt.NumChannels;
	outputWAVhdr.fmt.BlockAlign = oneChannelBlockAlign*outputWAVhdr.fmt.NumChannels;


	// Write output WAV header to file
	//-------------------------------------------------
	WriteWavHeader(wav_out,outputWAVhdr);


	// Processing loop
	//-------------------------------------------------	
	{
		DSPint sample;
		DSPint BytesPerSample = inputWAVhdr.fmt.BitsPerSample/8;
		const double SAMPLE_SCALE = -(double)(1 << 31);		//2^31
		DSPint iNumSamples = inputWAVhdr.data.SubChunk2Size/(inputWAVhdr.fmt.NumChannels*inputWAVhdr.fmt.BitsPerSample/8);
		
		// exact file length should be handled correctly...
		for(DSPint i=0; i<iNumSamples/BLOCK_SIZE; i++)
		{	
			for(DSPint j=0; j<BLOCK_SIZE; j++)
			{
				for(DSPint k=0; k<inputWAVhdr.fmt.NumChannels; k++)
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
				for(DSPint k=0; k<outputWAVhdr.fmt.NumChannels; k++)
				{	
					sample = sampleBuffer[k][j].toLong();	// crude, non-rounding 			
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

DSPfract second_order_IIR(DSPfract input, DSPfract* coefficients, DSPfract* x_history, DSPfract* y_history) 
{
	DSPaccum output = 0;

    output += *(coefficients++) * input;        /* A0 * x(n)     */
    output += *(coefficients) * *x_history; /* A1 * x(n-1) */
	output += *(coefficients++) * *(x_history++); /* A1 * x(n-1) */
    output += *(coefficients++) * *(x_history--); /* A2 * x(n-2)   */
    output -= *(coefficients) * *y_history; /* B1 * y(n-1) */
	output -= *(coefficients++)* *(y_history++); /* B1 * y(n-1) */
    output -= *(coefficients) * *(y_history--); /* B2 * y(n-2)   */

	DSPfract tmp_x = *x_history;
	DSPfract tmp_y = *y_history;
    
    *(++y_history) = tmp_y;	
	*(--y_history) = output;	
	*(++x_history) = tmp_x;	
	*(--x_history) = input;	

    return FRACT_NUM(output);
}

void processing()
{
	for(p = sampleBuffer[0]; p < sampleBuffer[0] + BLOCK_SIZE; p++)
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
		if(outputMode != 2)
		{
			for(p = sampleBuffer[0]; p < sampleBuffer[0] + BLOCK_SIZE; p++)
			{
				*(p + 64) = ZERO;
			}
			if(outputMode != 1)
			{
				for(p = sampleBuffer[0]; p < sampleBuffer[0] + BLOCK_SIZE; p++)
				{
					*(p + 32) = ZERO;
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
				*(p + 48) = ZERO;
			}
		}
	}
}