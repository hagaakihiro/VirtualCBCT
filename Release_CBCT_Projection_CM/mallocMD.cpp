/************************************************************************************************

	VOLUME-ONE Library (VOL):
	
		array.c:	memory allocation of multi-dimensional(1D-4D) array


	All Rights Reserved (C) 2002-	VOLUME-ONE developpers group

		2002.12.14:		first version for open source
		2004.03.05:		2nd release of VOL-1.1

************************************************************************************************/
#include <stdio.h>
#include <stdlib.h>
#if !defined(__APPLE__)
#include <malloc.h>
#endif
#include "mallocMD.h"

int
getByteSizeOfUnit(int unit)
{
	switch(unit)
	{
		case UNIT_UINT8:	return (int)sizeof(unsigned char);
		case UNIT_SINT8:	return (int)sizeof(char);
		case UNIT_UINT16:	return (int)sizeof(unsigned short);
		case UNIT_SINT16:	return (int)sizeof(short);
		case UNIT_UINT32:	return (int)sizeof(unsigned long);
		case UNIT_SINT32:	return (int)sizeof(long);
		case UNIT_FLOAT32:	return (int)sizeof(float);
		case UNIT_FLOAT64:	return (int)sizeof(double);
		default:		return 0;
	}
}

void*
new1DArray(int length,int unit)
{
	char*	ret=(char *)calloc(length,getByteSizeOfUnit(unit));

	return (void*)ret;
}

void**
new2DArray(int width, int height, int unit)
{
	char*	data1D;
	int	bytesize = getByteSizeOfUnit(unit);
	char**	data;

	if( (data1D=(char *)new1DArray(width*height,unit)) == NULL )
	{
		return(NULL);
	}

	/* make pointers for 2D */
	if( (data=(char **)malloc(height*sizeof(char*))) == NULL )
	{
		free(data1D);
		return(NULL);
	}

	for(int j=0; j<height ; j++)	data[j] = data1D + bytesize*j*width;

	return( (void**)data );
}

void***
new3DArray(int width, int height, int depth, int unit)
{
	char*		data1D;
	int		bytesize = getByteSizeOfUnit(unit);
	char***	data;
	int		wh = width*height;


	/* allocate as 1D */
	if( (data1D=(char *)new1DArray(wh*depth,unit)) == NULL )
	{	
		return(NULL);
	}

	/* make pointers for 3D */
	if( (data=(char ***)malloc(depth*sizeof(char**))) == NULL )
	{
		free(data1D);
		return(NULL);
	}

	for(int k=0; k<depth ; k++)
	{
		if( (data[k]=(char **)malloc(height*sizeof(char*))) == NULL )
		{
			free(data1D);
			for(int i=0;i<k;i++)	free(data[i]);
			free(data);
			return(NULL);
		}
	}

	for(int k=0; k<depth ; k++)
	{
		for(int j=0; j<height ; j++)
		{
			data[k][j] = data1D + bytesize*(k*wh + j*width);
		}
	}

	return((void***)data);
}

void****
new4DArray(int width, int height, int depth, int channel, int unit)
{
	char*		data1D;
	int		bytesize = getByteSizeOfUnit(unit);
	char****	data;
	int		wh = width*height;
	int		whd = wh*depth;


	/* allocate as 1D */
	if( (data1D=(char *)new1DArray(whd*channel,unit)) == NULL )
	{
		return(NULL);
	}

	/* make pointers for 4D */
	if( (data=(char ****)malloc(channel*sizeof(char**))) == NULL )
	{
		free(data1D);
		return(NULL);
	}

	for(int l=0; l<channel ; l++)
	{
		if( (data[l]=(char ***)malloc(depth*sizeof(char**))) == NULL )
		{
			free(data1D);
			for(int i=0;i<l;i++)	free(data[i]);
			free(data);
			return(NULL);
		}
	}

	for(int l=0; l<channel ; l++)
	{
		for(int k=0; k<depth ; k++)
		{
			if( (data[l][k]=(char **)malloc(height*sizeof(char*))) == NULL )
			{
				free(data1D);
				for(int j=0;j<channel;j++) for(int i=0;i<k;i++) free(data[j][i]);
        			for(int i=0;i<channel;i++) free(data[i]);
				free(data);
				return(NULL);
			}
		}
	}

	for(int l=0; l<channel ; l++)
		for(int k=0; k<depth ; k++)
			for(int j=0; j<height ; j++)
				data[l][k][j] = data1D + bytesize*(l*whd + k*wh + j*width);

	return((void****)data);
}

///////////////////////



void*
get1DArrayOf2DArray(void** data,int unit)
{
	switch(unit)
	{
		case UNIT_UINT8:	return((void*)&((unsigned char**)data)[0][0]);
		case UNIT_SINT8:	return((void*)&((char**)data)[0][0]);
		case UNIT_UINT16:	return((void*)&((unsigned short**)data)[0][0]);
		case UNIT_SINT16:	return((void*)&((short**)data)[0][0]);
		case UNIT_UINT32:	return((void*)&((unsigned long**)data)[0][0]);
		case UNIT_SINT32:	return((void*)&((long**)data)[0][0]);
		case UNIT_FLOAT32:	return((void*)&((float**)data)[0][0]);
		case UNIT_FLOAT64:	return((void*)&((double**)data)[0][0]);

		default:					return(NULL);
	}
}

void*
get1DArrayOf3DArray(void*** data,int unit)
{
	switch(unit)
	{
		case UNIT_UINT8:	return((void*)&((unsigned char***)data)[0][0][0]);
		case UNIT_SINT8:	return((void*)&((char***)data)[0][0][0]);
		case UNIT_UINT16:	return((void*)&((unsigned short***)data)[0][0][0]);
		case UNIT_SINT16:	return((void*)&((short***)data)[0][0][0]);
		case UNIT_UINT32:	return((void*)&((unsigned long***)data)[0][0][0]);
		case UNIT_SINT32:	return((void*)&((long***)data)[0][0][0]);
		case UNIT_FLOAT32:	return((void*)&((float***)data)[0][0][0]);
		case UNIT_FLOAT64:	return((void*)&((double***)data)[0][0][0]);

		default:					return(NULL);
	}
}


void*
get1DArrayOf4DArray(void**** data,int unit)
{
	switch(unit)
	{
		case UNIT_UINT8:	return((void*)&((unsigned char****)data)[0][0][0][0]);
		case UNIT_SINT8:	return((void*)&((char****)data)[0][0][0][0]);
		case UNIT_UINT16:	return((void*)&((unsigned short****)data)[0][0][0][0]);
		case UNIT_SINT16:	return((void*)&((short****)data)[0][0][0][0]);
		case UNIT_UINT32:	return((void*)&((unsigned long****)data)[0][0][0][0]);
		case UNIT_SINT32:	return((void*)&((long****)data)[0][0][0][0]);
		case UNIT_FLOAT32:	return((void*)&((float****)data)[0][0][0][0]);
		case UNIT_FLOAT64:	return((void*)&((double****)data)[0][0][0][0]);

		default:					return(NULL);
	}
}


///////////////////////

void
delete1DArray(void* data)
{
	free(data);
}

void
delete2DArray(void** data,int width,int height,int unit)
{
	if( data == NULL )	return;

	free((char*)get1DArrayOf2DArray(data,unit));

	free( data );

	return;
}

void
delete3DArray(void*** data,int width,int height,int depth,int unit)
{
	if( data == NULL )	return;
	
	free((char*)get1DArrayOf3DArray(data,unit));

	for(int i=0;i<depth;i++)	free( ((char***)data)[i] );

	free( data );

	return;
}

void
delete4DArray(void**** data,int width,int height,int depth,int channel,int unit)
{
	if( data == NULL )	return;

	free((char*)get1DArrayOf4DArray(data,unit));

	for(int j=0;j<channel;j++)
		for(int i=0;i<depth;i++) free( ((char****)data)[j][i] );

	for(int i=0;i<channel;i++)	free( ((char****)data)[i] );

	free( data );

	return;
}





