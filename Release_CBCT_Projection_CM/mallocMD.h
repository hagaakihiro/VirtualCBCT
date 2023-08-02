/************************************************************************************************

	VOLUME-ONE Library (VOL):
	
		array.h:	memory allocation of multi-dimensional(1D-4D) array


	All Rights Reserved (C) 2002-	VOLUME-ONE developpers group

		2002.12.14:		first version for open source
		2004.03.05:		2nd release of VOL-1.1

************************************************************************************************/

#ifndef ARRAY_H


#define UNIT_UINT8		0		// unsigned char
#define UNIT_SINT8		1		// (signed) char
#define UNIT_UINT16		2		// unsigned short
#define UNIT_SINT16		3		// (signed) short
#define UNIT_UINT32		4		// unsigned long
#define UNIT_SINT32		5		// (signed) long
#define UNIT_FLOAT32		6		// float
#define UNIT_FLOAT64		7		// double



int		getByteSizeOfUnit(int unit);

void*		new1DArray(int length, int unit);
void** 	new2DArray(int width, int height, int unit);
void***	new3DArray(int width, int height, int depth, int unit);
void****	new4DArray(int width, int height, int depth, int channel, int unit);

void*		get1DArrayOf2DArray(void** data,int unit);
void*		get1DArrayOf3DArray(void*** data,int unit);
void*		get1DArrayOf4DArray(void**** data,int unit);

void		delete1DArray(void* data);
void		delete2DArray(void** data,int width,int height,int unit);
void		delete3DArray(void*** data,int width,int height,int depth,int unit);
void		delete4DArray(void**** data,int width,int height,int depth,int channel,int unit);

#endif

#define ARRAY_H
