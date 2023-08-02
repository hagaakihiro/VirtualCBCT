#include <stdio.h>
#include <string.h>
#include <math.h>
#include <malloc.h> // for Ubuntu

#include "physParams.h"
#include "mallocMD.h"
#include "IR_ImageCBCT.h"

double double_TV(short a, short ax, short ay, short az, short ax1, short ay1, short az1, short ayx1, short azx1, short axy1, short azy1, short axz1, short ayz1)
{
  double dist1 = sqrt(double((a-ax)*(a-ax)+(a-ay)*(a-ay)+(a-az)*(a-az)))+0.00000001;
  double dist2 = sqrt(double((a-ax1)*(a-ax1)+(ax1-ayx1)*(ax1-ayx1)+(ax1-azx1)*(ax1-azx1)))+0.00000001;
  double dist3 = sqrt(double((ay1-axy1)*(ay1-axy1)+(a-ay1)*(a-ay1)+(ay1-azy1)*(ay1-azy1)))+0.00000001;
  double dist4 = sqrt(double((az1-axz1)*(az1-axz1)+(az1-ayz1)*(az1-ayz1)+(a-az1)*(a-az1)))+0.00000001;
  return ((a-ax)+(a-ay)+(a-az))/dist1+(a-ax1)/dist2+(a-ay1)/dist3+(a-az1)/dist4;
}

void
IR_ImageCBCT( const int ite, int reconSize, int reconSlices, double reconScale, short* reconImageshort_before, float* reconImagefloat,
              float* npf_float, float* npf_float_re, double wbefore,
              double weight_tv, double weight_prior, const int thinout)
{
  double rmax = reconSize*0.5 * 1.0;
#ifdef _OPENMP
#pragma omp parallel for
#endif
  for(int kkz = 0; kkz < reconSlices; kkz++)
    for(int kky = 0; kky < reconSize; kky++)
      for(int kkx = 0; kkx < reconSize; kkx++)
	{
	  double cost_tv = 0.0, cost_prior = 0.0;
	  int kk = kkz*reconSize*reconSize + kky*reconSize + kkx;
	  //TotalVariation
	  if(kkz < reconSlices-1 && kky < reconSize-1 && kkx < reconSize-1 && kkx > 0 && kky >  0 && kkz > 0)
	    {
	      int kx = kkz*reconSize*reconSize + kky*reconSize + kkx+1;
	      int ky = kkz*reconSize*reconSize + (kky+1)*reconSize + kkx;
	      int kz = (kkz+1)*reconSize*reconSize + kky*reconSize + kkx;
	      int kx1 = kkz*reconSize*reconSize + kky*reconSize + kkx-1;
	      int ky1 = kkz*reconSize*reconSize + (kky-1)*reconSize + kkx;
	      int kz1 = (kkz-1)*reconSize*reconSize + kky*reconSize + kkx;

	      int kyx1 = kkz*reconSize*reconSize + (kky+1)*reconSize + kkx-1;
	      int kzx1 = (kkz+1)*reconSize*reconSize + kky*reconSize + kkx-1;

	      int kxy1 = kkz*reconSize*reconSize + (kky-1)*reconSize + kkx+1;
	      int kzy1 = (kkz+1)*reconSize*reconSize + (kky-1)*reconSize + kkx;

	      int kxz1 = (kkz-1)*reconSize*reconSize + kky*reconSize + kkx+1;
	      int kyz1 = (kkz-1)*reconSize*reconSize + (kky+1)*reconSize + kkx;


	      cost_tv =  double_TV(reconImageshort_before[kk],
				    reconImageshort_before[kx],
				    reconImageshort_before[ky],
				    reconImageshort_before[kz],
				    reconImageshort_before[kx1],
				    reconImageshort_before[ky1],
				    reconImageshort_before[kz1],
				    reconImageshort_before[kyx1],
				    reconImageshort_before[kzx1],
				    reconImageshort_before[kxy1],
				    reconImageshort_before[kzy1],
				    reconImageshort_before[kxz1],
				    reconImageshort_before[kyz1]);
	      /*
	      cost_prior =  double_TV(reconImagefloat_before[kk]-Prior_image_2D[kk],
				       reconImagefloat_before[k1]-Prior_image_2D[k1],
				       reconImagefloat_before[k2]-Prior_image_2D[k2],
				       reconImagefloat_before[k3]-Prior_image_2D[k3],
				       reconImagefloat_before[k4]-Prior_image_2D[k4],
				       reconImagefloat_before[k7]-Prior_image_2D[k7],
				       reconImagefloat_before[k8]-Prior_image_2D[k8]);
	      */  //prior image PICCS
	    }
	  double lambda = weight_tv * cost_tv + weight_prior * cost_prior + 0.000001;//for penalty term
	  reconImagefloat[kk] = reconImageshort_before[kk]*(wbefore+(1-wbefore)*npf_float_re[kk]/(npf_float[kk]+lambda));//poisson likelyhood(MLEM)  2020/11/17 madeno yatu(gakkaidemotiiteita AAPMtoka)

	  //reconImagefloat[kk] = reconImageshort_before[kk]*(wbefore+(1-wbefore)*npf_float_re[kk]/(npf_float[kk]+lambda));//weighted square  2020/11/17 kara new made

	  //Gradient Decent
	  //double del_F = npf_float_re[kk] - npf_float[kk] + lambda;//20171023
	  //double del_F = npf_float_re[kk] - npf_float[kk]/thinout + lambda;//20171023 add
	  //reconImagefloat[kk] = (float)(reconImagefloat_before[kk] + wbefore * del_F);

	  if(reconImagefloat[kk] < 0.0 ) reconImagefloat[kk] = 0.0;
	  else if(reconImagefloat[kk] > 65535/2) reconImagefloat[kk] = 5000.0;
	  double xx = -reconSize*reconScale*0.5 + (kkx+0.5)*reconScale;
	  double yy = reconSize*reconScale*0.5 - (kky+0.5)*reconScale;
	  double r = sqrt((xx)*(xx)+(yy)*(yy));
	  if(r >= rmax) reconImagefloat[kk] = 0.0;
	}
}

