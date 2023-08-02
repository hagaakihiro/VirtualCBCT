
void
prior_weight_production( int pri_size, int reconSize, int reconSlices, int num_material, 
			 float* ProjImagefloat, float* npf_float, float* npf_float_re, float* projVolume,
			 int projImgWidth, int projImgHeight, int usedProjNumber, int startProjNumber, int endProjNumber,
			 double* anglesRad, double* cosTable, double* sinTable, double* XcenterShift, double* YcenterShift, 
			 short* req_reconImageshort, short* Pri_reconImageshort, short* Pri_reconImageshort_0, float* Attenuation, float* pdf_pe, int NE, double xfactor);

void
attenuation_coefficient( char* spectrum_data_name, int num_material, float* Attenuation, float* pdf_pe, int *NE, int nthin, char* mat_type);



