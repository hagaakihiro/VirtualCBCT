#include <thrust/host_vector.h>
#include <thrust/device_vector.h>
#include <iostream>
#include <stdio.h>
#include <string.h>
#include <math.h>
//#include <malloc.h> // for Ubuntu

#include "physParams.h"
#include "mallocMD.h"
#include "misc.h"






// kernel code

__constant__ double	devcosTable[360];
__constant__ double	devsinTable[360];
__constant__ double     devXcenterShift[360];
__constant__ double     devYcenterShift[360];

__constant__ int	devreconSize;
__constant__ double   	devreconScale;
__constant__ int	devreconSlices;
__constant__ double   	devreconStep;
__constant__ int	devprojImgHeight;
__constant__ int   	devprojImgWidth;
__constant__ double     devx_reg;
__constant__ double     devy_reg;
__constant__ double     devz_reg;

__constant__ double     devDIST_BTWN_SRC_AND_ISOCENT_CM;
__constant__ double     devDETECTOR_PITCH_CM_AT_ISO_CENT;
__constant__ double     devDETECTOR_PIXEL_NUM;
__constant__ double     devDETECTOR_PITCH_CM;

texture<float, cudaTextureType3D, cudaReadModeElementType> texX;
texture<float, cudaTextureType3D, cudaReadModeElementType> texM;
texture<float, cudaTextureType2D, cudaReadModeElementType> texVGPU;

#define CHECK(call, err)                                                    \
{                                                                               \
  const cudaError_t error = call;					\
  err = error;                                    \
  if(error != cudaSuccess)                                            \
  {                                                                            \
    fprintf(stderr, "Error: %s:%d, ",__FILE__, __LINE__);           \
    fprintf(stderr, "code: %d, reason: %s\n",error,              \
      cudaGetErrorString(error));                                   \
  }                                                                           \
}

__global__ void
reprojectionKernel(float* devprojVolume,
		   unsigned short* devreconImageshort, double xfactor,
		   float* devAttenuation, float* devpdf_pe, int num_mat, int NE, int ke)
{
  int nxm = devreconSize*0.5;
  int nym = devreconSize*0.5;
  int nzm = devreconSlices*0.5;
  int pri_size = devreconSize*devreconSize;

  double rmax = devreconSize*0.5 * devreconScale;
  //double rmax = devreconSlices*0.5 * 1.0;

  /*
//Haga
  ////////////////////////
  int	jjj = blockDim.x;   // gantry angle
  int	jj = blockIdx.y;    // detector height
  int   j  = blockIdx.x;    // detector width

  /////////////////////////////
  */

  /*
 //slow
  ////////////////////////
  int	jjj =threadIdx.x;   // gantry angle
  int	jj = blockIdx.y;    // detector height
  int   j  = blockIdx.x;    // detector width

  /////////////////////////////
  */

  
  ///////////
   int	jjj = blockIdx.z;   // gantry angle
   int   jj = blockIdx.y;    // detector width
   //int   thre  = threadIdx.x;    // detector height
   int j=blockIdx.x*blockDim.x+threadIdx.x;
   /////////////
   

  /*
  //kore
   int	jjj = blockIdx.y;   // gantry angle
   int   jj = blockIdx.x;    // detector width
   int   j  = threadIdx.x;    // detector height
  */


  int jjjj = jjj*devprojImgHeight*devprojImgWidth + (devprojImgHeight-1-j)*devprojImgWidth +  jj; // 18(spot)/36(all) is match to MonteCarlo 0 angle position

  
   //int jjjj = jjj*devprojImgHeight*devprojImgWidth + j*devprojImgWidth +  jj;//down TOP OF HEAD //this is match to CPU

   
   //int jjjj = jjj*devprojImgHeight*devprojImgWidth + j*devprojImgWidth + (devprojImgWidth-1- jj);//down TOP OF HEAD
   
   //int	jjj = blockIdx.y;   // gantry angle
   //int   jj = blockIdx.x;    // detector height
   //int   j  = threadIdx.x;    // detector width  
   //int jjjj = jjj*devprojImgHeight*devprojImgWidth + jj*devprojImgWidth + (devprojImgWidth-1- j);//left TOP OF HEAD
   
   // int jjjj = jjj*devprojImgHeight*devprojImgWidth + jj*devprojImgWidth + j;//RIGHT TOP OF HEAD
  // int jjjj = jjj*512*512 + jj*512 + j;


  devprojVolume[jjjj]=0.0;

  
  
  // Source location in cartesian coodinate
  double sx = -devDIST_BTWN_SRC_AND_ISOCENT_CM*devsinTable[jjj] + devx_reg;
  double sy = devDIST_BTWN_SRC_AND_ISOCENT_CM*devcosTable[jjj] + devy_reg;
  double sz = devz_reg;
  //double sz = 0.0;


  double pz =  -(devYcenterShift[jjj] + j) * devDETECTOR_PITCH_CM_AT_ISO_CENT * (devDETECTOR_PIXEL_NUM / devprojImgHeight) + devz_reg*devDETECTOR_PITCH_CM_AT_ISO_CENT/devDETECTOR_PITCH_CM;
  //meaning:pz =DETECTOR_PITCH_CM_AT_ISO_CENT*-(YcenterShift[jjj]+j)*(800/400)
  double pr = devDETECTOR_PITCH_CM_AT_ISO_CENT * (devXcenterShift[jjj] + jj) * (devDETECTOR_PIXEL_NUM / devprojImgHeight);
  //meaning:pr =DETECTOR_PITCH_CM_AT_ISO_CENT*(XcenterShift+jj)*(800/400)
  
  double px = pr*devcosTable[jjj] + devx_reg*devDETECTOR_PITCH_CM_AT_ISO_CENT/devDETECTOR_PITCH_CM;
  double py = pr*devsinTable[jjj] + devy_reg*devDETECTOR_PITCH_CM_AT_ISO_CENT/devDETECTOR_PITCH_CM;
  double leng = devDIST_BTWN_SRC_AND_ISOCENT_CM*devDIST_BTWN_SRC_AND_ISOCENT_CM/(devDIST_BTWN_SRC_AND_ISOCENT_CM*devDIST_BTWN_SRC_AND_ISOCENT_CM + (pr*pr+pz*pz));
  //meaning:leng =1000*1000/(1000*1000+(pr*pr+pz*pz))
  
  //double dex=devprojImgWidth/2-jj;
  //double dey=-(devprojImgHeight-1-j)+devprojImgHeight/2;


  
  //double leng = 1538*1538/(1538*1538 + (dex*dex+dey*dey));
  leng = leng*sqrt(leng);
  //leng = leng*leng;

    
  double xx[4], yy[4], zz[4];
  int m = 0;
  yy[m] = double(nym)*devreconScale;
  if(fabs((yy[m]-sy)/(py-sy)*(px-sx)+sx) < double(nxm)*devreconScale)
    { 
      xx[m] = (yy[m]-sy)/(py-sy)*(px-sx)+sx;
      zz[m] = (yy[m]-sy)/(py-sy)*(pz-sz)+sz;
      m++;
    }
  yy[m] = -double(nym)*devreconScale;
  if(fabs((yy[m]-sy)/(py-sy)*(px-sx)+sx) < double(nxm)*devreconScale)
    {
      xx[m] = (yy[m]-sy)/(py-sy)*(px-sx)+sx;
      zz[m] = (yy[m]-sy)/(py-sy)*(pz-sz)+sz;
      m++;
    }
  xx[m] = double(nxm)*devreconScale;
  if(fabs((xx[m]-sx)/(px-sx)*(py-sy)+sy) < double(nym)*devreconScale)
    { 
      yy[m] = (xx[m]-sx)/(px-sx)*(py-sy)+sy;
      zz[m] = (xx[m]-sx)/(px-sx)*(pz-sz)+sz;
      m++;
    }
  xx[m] = -double(nxm)*devreconScale;
  if(fabs((xx[m]-sx)/(px-sx)*(py-sy)+sy) < double(nym)*devreconScale)
    {
      yy[m] = (xx[m]-sx)/(px-sx)*(py-sy)+sy;
      zz[m] = (xx[m]-sx)/(px-sx)*(pz-sz)+sz;
      m++;
    }
  double xi, yi, zi, xf, yf, zf;
  if((xx[0]-sx)*(xx[0]-sx)+(yy[0]-sy)*(yy[0]-sy)+(zz[0]-sz)*(zz[0]-sz) < (xx[1]-sx)*(xx[1]-sx)+(yy[1]-sy)*(yy[1]-sy)+(zz[1]-sz)*(zz[1]-sz))
    {
      xi=xx[0]; yi=yy[0]; zi=zz[0]; xf=xx[1]; yf=yy[1]; zf=zz[1];
    }
  else
    {
      xi=xx[1]; yi=yy[1]; zi=zz[1]; xf=xx[0]; yf=yy[0]; zf=zz[0];
    }
  //fprintf(stderr,"%lf %lf %lf %lf %lf %lf\n", xi, yi, zi, xf, yf, zf);
  double xb0 = xi, yb0 = yi, zb0 = zi, xb1, yb1, zb1;
  double xii, yii, zii, xiii, yiii, ziii, xiiii, yiiii, ziiii, yj, xj, zj, rr, transp, total_length;
  int xma, yma, zma, xmb, ymb, zmb, xmc, ymc, zmc, icheck, vnum;

  if(fabs(xi-xf) > fabs(yi-yf)) // dx is larger than dy
    //if(fabs(xi-xf) >= fabs(yi-yf)) // dx is larger than dy
    { 
      int iks = 0;
      double pha = (xf-xi)/fabs(xi-xf);
      int xin, ik=0;
      double chy, chz;
      if((xi < 0 && (xf-xi) < 0) || (xi > 0 && (xf-xi) > 0)) xin = int(xi) + (xf-xi)/fabs(xf-xi);
      else xin = int(xi);
      xii = xin;
      yii = (xii-sx)/(px-sx)*(py-sy)+sy;
      zii = (xii-sx)/(px-sx)*(pz-sz)+sz;
      transp = 0.0;
      total_length = 0.0;
      while( xii <= double(nxm)*devreconScale && xii >= -double(nxm)*devreconScale &&  yii <= double(nym)*devreconScale && yii >= -double(nym)*devreconScale)
	{
	  xii = xin + pha * ik*devreconScale;
	  ymb=yma;
	  yj = yii;
	  yii = (xii-sx)/(px-sx)*(py-sy)+sy;
	  yma=(int)(yii/devreconScale);
	  zmb=zma;
	  zj = zii;
	  zii = (xii-sx)/(px-sx)*(pz-sz)+sz;
	  zma=(int)(zii/devreconStep);
	  icheck = 0;
	  chy = 0;
	  chz = 0;
	  if((ik != 0 && yma != ymb) || (ik != 0 && (yii/fabs(yii) != yj/fabs(yj))))
	    {
	      //if((yii < 0 && (yf-yi) < 0) || (yii > 0 && (yf-yi) > 0)) yiii = (double)(yma*devreconScale);
	      //else yiii = (double)(yma*devreconScale) - (yf-yi)/fabs(yf-yi)*devreconScale;
	      if((yii <= 0 && (yf-yi) < 0) || (yii >= 0 && (yf-yi) > 0)) yiii = (double)(yma*devreconScale);//shimo 2021/12/02
	      else yiii = (double)(yma*devreconScale) - (yf-yi)/fabs(yf-yi)*devreconScale;//shimo 2021/12/02
	      xiii = (yiii-sy)/(py-sy)*(px-sx)+sx;
	      ziii = (yiii-sy)/(py-sy)*(pz-sz)+sz;
	      chy = (xii-xiii)*(xii-xiii)+(yii-yiii)*(yii-yiii)+(zii-ziii)*(zii-ziii);
	      icheck++;
	    }
	  //if((ik != 0 && zma != zmb) || (ik != 0 && (zii/fabs(zii) != zj/fabs(zj))))
	  if(ik != 0 && zma != zmb)
            {
	      if((zii < 0 && (zf-zi) < 0) || (zii > 0 && (zf-zi) > 0)) ziiii = (double)(zma*devreconStep);
	      else ziiii = (double)(zma*devreconStep) - (zf-zi)/fabs(zf-zi)*devreconStep;
	      //ziiii = (double)(zma*devreconStep);
	      xiiii = (ziiii-sz)/(pz-sz)*(px-sx)+sx;
	      yiiii = (ziiii-sz)/(pz-sz)*(py-sy)+sy;
	      chz = (xii-xiiii)*(xii-xiiii)+(yii-yiiii)*(yii-yiiii)+(zii-ziiii)*(zii-ziiii);
	      icheck++;
            }

	  if(icheck == 2 && chz > chy)
            {
	      iks++;
	      xb1=xiiii;
	      yb1=yiiii;
	      zb1=ziiii;
	      rr = sqrt((xb0-xb1)*(xb0-xb1)+(yb0-yb1)*(yb0-yb1)+(zb0-zb1)*(zb0-zb1));
	      xmc = (int)((xb0+xb1)*0.5/devreconScale+nxm);
	      ymc = (int)((-yb0-yb1)*0.5/devreconScale+nym);
	      zmc = (int)((-zb0-zb1)*0.5/devreconStep+nzm);
	      if(xmc >= 0 && ymc >= 0 && zmc >= 0 && xmc <= devreconSize-1 && ymc <= devreconSize-1 && zmc <= devreconSlices-1)
                {
		  double r = sqrt((xb0+xb1)*(xb0+xb1)+(-yb0-yb1)*(-yb0-yb1))*0.5;
		  double ph = 1.0;
		  for(int mat = 0; mat < num_mat; mat++)
		    {
		      double ph = 1.0;
		      vnum = xmc + ymc*devreconSize + (zmc*num_mat+mat)*devreconSize*devreconSize;
		      if(devreconImageshort[vnum] < 0 || r >= rmax) ph = 0.0;
		      transp += devreconImageshort[vnum] * devAttenuation[mat*NE + ke] * rr * ph;
		    }
 		    total_length += rr;
                }
	      xb0=xb1;
	      yb0=yb1;
	      zb0=zb1;
	      iks++;
	      xb1=xiii;
	      yb1=yiii;
	      zb1=ziii;
	      rr = sqrt((xb0-xb1)*(xb0-xb1)+(yb0-yb1)*(yb0-yb1)+(zb0-zb1)*(zb0-zb1));
	      xmc = (int)((xb0+xb1)*0.5/devreconScale+nxm);
	      ymc = (int)((-yb0-yb1)*0.5/devreconScale+nym);
	      zmc = (int)((-zb0-zb1)*0.5/devreconStep+nzm);
	      if(xmc >= 0 && ymc >= 0 && zmc >= 0 && xmc <= devreconSize-1 && ymc <= devreconSize-1 && zmc <= devreconSlices-1)
		{
		  double r = sqrt((xb0+xb1)*(xb0+xb1)+(-yb0-yb1)*(-yb0-yb1))*0.5;
		  double ph = 1.0;
		  for(int mat = 0; mat < num_mat; mat++)
		    {
		      double ph = 1.0;
		      vnum = xmc + ymc*devreconSize + (zmc*num_mat+mat)*devreconSize*devreconSize;
		      if(devreconImageshort[vnum] < 0 || r >= rmax) ph = 0.0;
		      transp += devreconImageshort[vnum] * devAttenuation[mat*NE + ke] * rr * ph;
		    }
 		  total_length += rr;
		}
	      xb0=xb1;
	      yb0=yb1;
	      zb0=zb1;
	    }
	  else if(icheck == 2 && chz < chy)
	    {
	      iks++;
	      xb1=xiii;
	      yb1=yiii;
	      zb1=ziii;
	      rr = sqrt((xb0-xb1)*(xb0-xb1)+(yb0-yb1)*(yb0-yb1)+(zb0-zb1)*(zb0-zb1));
	      xmc = (int)((xb0+xb1)*0.5/devreconScale+nxm);
	      ymc = (int)((-yb0-yb1)*0.5/devreconScale+nym);
	      zmc = (int)((-zb0-zb1)*0.5/devreconStep+nzm);
	      if(xmc >= 0 && ymc >= 0 && zmc >= 0 && xmc <= devreconSize-1 && ymc <= devreconSize-1 && zmc <= devreconSlices-1)
	      {
		double r = sqrt((xb0+xb1)*(xb0+xb1)+(-yb0-yb1)*(-yb0-yb1))*0.5;
		double ph = 1.0;
		for(int mat = 0; mat < num_mat; mat++)
		  {
		    double ph = 1.0;
		    vnum = xmc + ymc*devreconSize + (zmc*num_mat+mat)*devreconSize*devreconSize;
		    if(devreconImageshort[vnum] < 0 || r >= rmax) ph = 0.0;
		    transp += devreconImageshort[vnum] * devAttenuation[mat*NE + ke] * rr * ph;
		  }
		total_length += rr;
	      }
	      xb0=xb1;
	      yb0=yb1;
	      zb0=zb1;
	      iks++;
	      xb1=xiiii;
              yb1=yiiii;
              zb1=ziiii;
              rr = sqrt((xb0-xb1)*(xb0-xb1)+(yb0-yb1)*(yb0-yb1)+(zb0-zb1)*(zb0-zb1));
	      xmc = (int)((xb0+xb1)*0.5/devreconScale+nxm);
	      ymc = (int)((-yb0-yb1)*0.5/devreconScale+nym);
	      zmc = (int)((-zb0-zb1)*0.5/devreconStep+nzm);
              if(xmc >= 0 && ymc >= 0 && zmc >= 0 && xmc <= devreconSize-1 && ymc <= devreconSize-1 && zmc <= devreconSlices-1)
		{
		  double r = sqrt((xb0+xb1)*(xb0+xb1)+(-yb0-yb1)*(-yb0-yb1))*0.5;
		  double ph = 1.0;
		  for(int mat = 0; mat < num_mat; mat++)
		    {
		      double ph = 1.0;
		      vnum = xmc + ymc*devreconSize + (zmc*num_mat+mat)*devreconSize*devreconSize;
		      if(devreconImageshort[vnum] < 0 || r >= rmax) ph = 0.0;
		      transp += devreconImageshort[vnum] * devAttenuation[mat*NE + ke] * rr * ph;
		    }
		  total_length += rr;
              }
              xb0=xb1;
              yb0=yb1;
              zb0=zb1;
	    }
	  else if( icheck == 1 && chz != 0)
	    {
	      iks++;
	      xb1=xiiii;
	      yb1=yiiii;
	      zb1=ziiii;
	      rr = sqrt((xb0-xb1)*(xb0-xb1)+(yb0-yb1)*(yb0-yb1)+(zb0-zb1)*(zb0-zb1));
	      xmc = (int)((xb0+xb1)*0.5/devreconScale+nxm);
	      ymc = (int)((-yb0-yb1)*0.5/devreconScale+nym);
	      zmc = (int)((-zb0-zb1)*0.5/devreconStep+nzm);
	      if(xmc >= 0 && ymc >= 0 && zmc >= 0 && xmc <= devreconSize-1 && ymc <= devreconSize-1 && zmc <= devreconSlices-1)
		{
		  double r = sqrt((xb0+xb1)*(xb0+xb1)+(-yb0-yb1)*(-yb0-yb1))*0.5;
		  double ph = 1.0;
		  for(int mat = 0; mat < num_mat; mat++)
		    {
		      double ph = 1.0;
		      vnum = xmc + ymc*devreconSize + (zmc*num_mat+mat)*devreconSize*devreconSize;
		      if(devreconImageshort[vnum] < 0 || r >= rmax) ph = 0.0;
		      transp += devreconImageshort[vnum] * devAttenuation[mat*NE + ke] * rr * ph;
		    }
		  total_length += rr;
		}
	      xb0=xb1;
	      yb0=yb1;
	      zb0=zb1;
	    }
	  else if( icheck == 1 && chy != 0)
	    {
	      iks++;
	      xb1=xiii;
	      yb1=yiii;
	      zb1=ziii;
	      rr = sqrt((xb0-xb1)*(xb0-xb1)+(yb0-yb1)*(yb0-yb1)+(zb0-zb1)*(zb0-zb1));
	      xmc = (int)((xb0+xb1)*0.5/devreconScale+nxm);
	      ymc = (int)((-yb0-yb1)*0.5/devreconScale+nym);
	      zmc = (int)((-zb0-zb1)*0.5/devreconStep+nzm);
	      if(xmc >= 0 && ymc >= 0 && zmc >= 0 && xmc <= devreconSize-1 && ymc <= devreconSize-1 && zmc <= devreconSlices-1)
		{
		  double r = sqrt((xb0+xb1)*(xb0+xb1)+(-yb0-yb1)*(-yb0-yb1))*0.5;
		  for(int mat = 0; mat < num_mat; mat++)
		    {
		      double ph = 1.0;
		      vnum = xmc + ymc*devreconSize + (zmc*num_mat+mat)*devreconSize*devreconSize;
		      if(devreconImageshort[vnum] < 0 || r >= rmax) ph = 0.0;
		      transp += devreconImageshort[vnum] * devAttenuation[mat*NE + ke] * rr * ph;
		    }
	  
		  total_length += rr;
		}
	      xb0=xb1;
	      yb0=yb1;
	      zb0=zb1;
	    }
	  iks++;
	  xb1=xii;
	  yb1=yii;
	  zb1=zii;
	  rr = sqrt((xb0-xb1)*(xb0-xb1)+(yb0-yb1)*(yb0-yb1)+(zb0-zb1)*(zb0-zb1));
	  xmc = (int)((xb0+xb1)*0.5/devreconScale+nxm);
	  ymc = (int)((-yb0-yb1)*0.5/devreconScale+nym);
	  zmc = (int)((-zb0-zb1)*0.5/devreconStep+nzm);
	  if(xmc >= 0 && ymc >= 0 && zmc >= 0 && xmc <= devreconSize-1 && ymc <= devreconSize-1 && zmc <= devreconSlices-1)
	    {
	      double r = sqrt((xb0+xb1)*(xb0+xb1)+(-yb0-yb1)*(-yb0-yb1))*0.5;
	      for(int mat = 0; mat < num_mat; mat++)
		{
		  double ph = 1.0;
		  vnum = xmc + ymc*devreconSize + (zmc*num_mat+mat)*devreconSize*devreconSize;
		  if(devreconImageshort[vnum] < 0 || r >= rmax) ph = 0.0;
		  transp += devreconImageshort[vnum] * devAttenuation[mat*NE + ke] * rr * ph;
		}
	      total_length += rr;
	    }
	  xb0=xb1;
	  yb0=yb1;
	  zb0=zb1;
	  ik++;
	} // while
      //double intens = devpdf_pe[ke] * exp(-transp/xfactor); // * photon count n_i/n_0;  *leng;

      double intens = devpdf_pe[ke] * exp(-transp/xfactor)*leng;//2021/6/12  // * photon count n_i/n_0;

      //double intens = MAX_AIR_PIXDIFF  * devpdf_pe[ke] * exp(-transp/xfactor)*leng;//2022/6/20  // * photon count n_i/n_0;
 
      //double intens = MAX_PIXVALUE - (MAX_PIXVALUE * exp(-transp/xfactor)*leng);//2022/1/17  //
      //double intens = MAX_PIXVALUE - (MAX_PIXVALUE* devpdf_pe[ke] * exp(-transp/xfactor)*leng);//2022/1/17  //
	devprojVolume[jjjj] = (float)(intens);
    } // dx dy
  
   else     // dy is larger than dx
  //else if(fabs(xi-xf) < fabs(yi-yf))     // dy is larger than dx
    {
      int iks = 0;
      double pha = (yf-yi)/fabs(yi-yf);
      int yin, ik=0;
      double chy, chz;
      if((yi < 0 && (yf-yi) < 0) || (yi > 0 && (yf-yi) > 0)) yin = int(yi) + (yf-yi)/fabs(yf-yi);
      else yin = int(yi);
      yii = yin;
      xii = (yii-sy)/(py-sy)*(px-sx)+sx;
      zii = (yii-sy)/(py-sy)*(pz-sz)+sz;
      transp = 0.0;
      total_length = 0.0;
      while( xii <= double(nxm)*devreconScale && xii >= -double(nxm)*devreconScale &&  yii <= double(nym)*devreconScale && yii >= -double(nym)*devreconScale)
	{
	  yii = yin + pha * ik*devreconScale;
	  xmb=xma;
	  xj = xii;
	  xii = (yii-sy)/(py-sy)*(px-sx)+sx;
	  xma=(int)(xii/devreconScale);
	  zmb=zma;
	  zj = zii;
	  zii = (yii-sy)/(py-sy)*(pz-sz)+sz;
	  zma=(int)(zii/devreconStep);
	  icheck = 0;
	  chy = 0;
	  chz = 0;
	  if((ik != 0 && xma != xmb) || (ik != 0 && (xii/fabs(xii) != xj/fabs(xj))))
	    {
	      //if((xii < 0 && (xf-xi) < 0) || (xii > 0 && (xf-xi) > 0)) xiii = (double)(xma*devreconScale);
	      //else xiii = (double)(xma*devreconScale) - (xf-xi)/fabs(xf-xi)*devreconScale;
	      if((xii <= 0 && (xf-xi) < 0) || (xii >= 0 && (xf-xi) > 0)) xiii = (double)(xma*devreconScale);//shimo 2021/12/02
	      else xiii = (double)(xma*devreconScale) - (xf-xi)/fabs(xf-xi)*devreconScale;//shimo 2021/12/02
	      yiii = (xiii-sx)/(px-sx)*(py-sy)+sy;
	      ziii = (xiii-sx)/(px-sx)*(pz-sz)+sz;
	      chy = (xii-xiii)*(xii-xiii)+(yii-yiii)*(yii-yiii)+(zii-ziii)*(zii-ziii);
	      icheck++;
	    }
	  //if((ik != 0 && zma != zmb) || (ik != 0 && (zii/fabs(zii) != zj/fabs(zj))))
	  if(ik != 0 && zma != zmb)
	    {
	      if((zii < 0 && (zf-zi) < 0) || (zii > 0 && (zf-zi) > 0)) ziiii = (double)(zma*devreconStep);
	      else ziiii = (double)(zma*devreconStep) - (zf-zi)/fabs(zf-zi)*devreconStep;
	      //ziiii = (double)(zma*devreconStep);
	      xiiii = (ziiii-sz)/(pz-sz)*(px-sx)+sx;
	      yiiii = (ziiii-sz)/(pz-sz)*(py-sy)+sy;
	      chz = (xii-xiiii)*(xii-xiiii)+(yii-yiiii)*(yii-yiiii)+(zii-ziiii)*(zii-ziiii);
	      icheck++;
	    }
	  
	  if(icheck == 2 && chz > chy)
	    {
	      iks++;
	      xb1=xiiii;
	      yb1=yiiii;
	      zb1=ziiii;
	      rr = sqrt((xb0-xb1)*(xb0-xb1)+(yb0-yb1)*(yb0-yb1)+(zb0-zb1)*(zb0-zb1));
	      xmc = (int)((xb0+xb1)*0.5/devreconScale+nxm);
	      ymc = (int)((-yb0-yb1)*0.5/devreconScale+nym);
	      zmc = (int)((-zb0-zb1)*0.5/devreconStep+nzm);
	      if(xmc >= 0 && ymc >= 0 && zmc >= 0 && xmc <= devreconSize-1 && ymc <= devreconSize-1 && zmc <= devreconSlices-1)
	      {
		double r = sqrt((xb0+xb1)*(xb0+xb1)+(-yb0-yb1)*(-yb0-yb1))*0.5;
		for(int mat = 0; mat < num_mat; mat++)
		  {
		    double ph = 1.0;
		    vnum = xmc + ymc*devreconSize + (zmc*num_mat+mat)*devreconSize*devreconSize;
		    if(devreconImageshort[vnum] < 0 || r >= rmax) ph = 0.0;
		    transp += devreconImageshort[vnum] * devAttenuation[mat*NE + ke] * rr * ph;
		    
		  }
		total_length += rr;
	      }
	      xb0=xb1;
	      yb0=yb1;
	      zb0=zb1;
	      iks++;
	      xb1=xiii;
	      yb1=yiii;
	      zb1=ziii;
	      rr = sqrt((xb0-xb1)*(xb0-xb1)+(yb0-yb1)*(yb0-yb1)+(zb0-zb1)*(zb0-zb1));
	      xmc = (int)((xb0+xb1)*0.5/devreconScale+nxm);
	      ymc = (int)((-yb0-yb1)*0.5/devreconScale+nym);
	      zmc = (int)((-zb0-zb1)*0.5/devreconStep+nzm);
	      if(xmc >= 0 && ymc >= 0 && zmc >= 0 && xmc <= devreconSize-1 && ymc <= devreconSize-1 && zmc <= devreconSlices-1)
	      {
		double r = sqrt((xb0+xb1)*(xb0+xb1)+(-yb0-yb1)*(-yb0-yb1))*0.5;
		for(int mat = 0; mat < num_mat; mat++)
		  {
		    double ph = 1.0;
		    vnum = xmc + ymc*devreconSize + (zmc*num_mat+mat)*devreconSize*devreconSize;
		    if(devreconImageshort[vnum] < 0 || r >= rmax) ph = 0.0;
		    transp += devreconImageshort[vnum] * devAttenuation[mat*NE + ke] * rr * ph;
		  }
		total_length += rr;
	      }
	      xb0=xb1;
	      yb0=yb1;
	      zb0=zb1;
	    }
	  else if(icheck == 2 && chz < chy)
            {
	      iks++;
	      xb1=xiii;
	      yb1=yiii;
	      zb1=ziii;
	      rr = sqrt((xb0-xb1)*(xb0-xb1)+(yb0-yb1)*(yb0-yb1)+(zb0-zb1)*(zb0-zb1));
	      xmc = (int)((xb0+xb1)*0.5/devreconScale+nxm);
	      ymc = (int)((-yb0-yb1)*0.5/devreconScale+nym);
	      zmc = (int)((-zb0-zb1)*0.5/devreconStep+nzm);
	      if(xmc >= 0 && ymc >= 0 && zmc >= 0 && xmc <= devreconSize-1 && ymc <= devreconSize-1 && zmc <= devreconSlices-1)
		{
		  double r = sqrt((xb0+xb1)*(xb0+xb1)+(-yb0-yb1)*(-yb0-yb1))*0.5;
		  for(int mat = 0; mat < num_mat; mat++)
		    {
		      double ph = 1.0;
		      vnum = xmc + ymc*devreconSize + (zmc*num_mat+mat)*devreconSize*devreconSize;
		      if(devreconImageshort[vnum] < 0 || r >= rmax) ph = 0.0;
		      transp += devreconImageshort[vnum] * devAttenuation[mat*NE + ke] * rr * ph;
		    }
		  total_length += rr;
		}
	      xb0=xb1;
	      yb0=yb1;
	      zb0=zb1;
	      iks++;
	      xb1=xiiii;
	      yb1=yiiii;
	      zb1=ziiii;
	      rr = sqrt((xb0-xb1)*(xb0-xb1)+(yb0-yb1)*(yb0-yb1)+(zb0-zb1)*(zb0-zb1));
	      xmc = (int)((xb0+xb1)*0.5/devreconScale+nxm);
	      ymc = (int)((-yb0-yb1)*0.5/devreconScale+nym);
	      zmc = (int)((-zb0-zb1)*0.5/devreconStep+nzm);
	      if(xmc >= 0 && ymc >= 0 && zmc >= 0 && xmc <= devreconSize-1 && ymc <= devreconSize-1 && zmc <= devreconSlices-1)
		{
		  double r = sqrt((xb0+xb1)*(xb0+xb1)+(-yb0-yb1)*(-yb0-yb1))*0.5;
		  for(int mat = 0; mat < num_mat; mat++)
		    {
		      double ph = 1.0;
		      vnum = xmc + ymc*devreconSize + (zmc*num_mat+mat)*devreconSize*devreconSize;
		      if(devreconImageshort[vnum] < 0 || r >= rmax) ph = 0.0;
		      transp += devreconImageshort[vnum] * devAttenuation[mat*NE + ke] * rr * ph;
		    }
		  total_length += rr;
		}
	      xb0=xb1;
	      yb0=yb1;
	      zb0=zb1;
	    }
	  else if( icheck == 1 && chz != 0)
	    {
	      iks++;
	      xb1=xiiii;
	      yb1=yiiii;
	      zb1=ziiii;
	      rr = sqrt((xb0-xb1)*(xb0-xb1)+(yb0-yb1)*(yb0-yb1)+(zb0-zb1)*(zb0-zb1));
	      xmc = (int)((xb0+xb1)*0.5/devreconScale+nxm);
	      ymc = (int)((-yb0-yb1)*0.5/devreconScale+nym);
	      zmc = (int)((-zb0-zb1)*0.5/devreconStep+nzm);
	      if(xmc >= 0 && ymc >= 0 && zmc >= 0 && xmc <= devreconSize-1 && ymc <= devreconSize-1 && zmc <= devreconSlices-1)
	      {
		double r = sqrt((xb0+xb1)*(xb0+xb1)+(-yb0-yb1)*(-yb0-yb1))*0.5;
		for(int mat = 0; mat < num_mat; mat++)
		  {
		    double ph = 1.0;
		    vnum = xmc + ymc*devreconSize + (zmc*num_mat+mat)*devreconSize*devreconSize;
		    if(devreconImageshort[vnum] < 0 || r >= rmax) ph = 0.0;
		    transp += devreconImageshort[vnum] * devAttenuation[mat*NE + ke] * rr * ph;
		  }
		total_length += rr;
	      }
	      xb0=xb1;
	      yb0=yb1;
	      zb0=zb1;
	    }
	  else if( icheck == 1 && chy != 0)
	    {
	      iks++;
	      xb1=xiii;
	      yb1=yiii;
	      zb1=ziii;
	      rr = sqrt((xb0-xb1)*(xb0-xb1)+(yb0-yb1)*(yb0-yb1)+(zb0-zb1)*(zb0-zb1));
	      xmc = (int)((xb0+xb1)*0.5/devreconScale+nxm);
	      ymc = (int)((-yb0-yb1)*0.5/devreconScale+nym);
	      zmc = (int)((-zb0-zb1)*0.5/devreconStep+nzm);
	      if(xmc >= 0 && ymc >= 0 && zmc >= 0 && xmc <= devreconSize-1 && ymc <= devreconSize-1 && zmc <= devreconSlices-1)
	      {
		double r = sqrt((xb0+xb1)*(xb0+xb1)+(-yb0-yb1)*(-yb0-yb1))*0.5;
		for(int mat = 0; mat < num_mat; mat++)
		  {
		    double ph = 1.0;
		    vnum = xmc + ymc*devreconSize + (zmc*num_mat+mat)*devreconSize*devreconSize;
		    if(devreconImageshort[vnum] < 0 || r >= rmax) ph = 0.0;
		    transp += devreconImageshort[vnum] * devAttenuation[mat*NE + ke] * rr * ph;
		  }
		total_length += rr;
	      }
	      xb0=xb1;
	      yb0=yb1;
	      zb0=zb1;
	    }
	  iks++;
	  xb1=xii;
	  yb1=yii;
	  zb1=zii;
	  rr = sqrt((xb0-xb1)*(xb0-xb1)+(yb0-yb1)*(yb0-yb1)+(zb0-zb1)*(zb0-zb1));
	  xmc = (int)((xb0+xb1)*0.5/devreconScale+nxm);
	  ymc = (int)((-yb0-yb1)*0.5/devreconScale+nym);
	  zmc = (int)((-zb0-zb1)*0.5/devreconStep+nzm);
	  if(xmc >= 0 && ymc >= 0 && zmc >= 0 && xmc <= devreconSize-1 && ymc <= devreconSize-1 && zmc <= devreconSlices-1)
	    {
	      double r = sqrt((xb0+xb1)*(xb0+xb1)+(-yb0-yb1)*(-yb0-yb1))*0.5;
	      for(int mat = 0; mat < num_mat; mat++)
		{
		  double ph = 1.0;
		  vnum = xmc + ymc*devreconSize + (zmc*num_mat+mat)*devreconSize*devreconSize;
		  if(devreconImageshort[vnum] < 0 || r >= rmax) ph = 0.0;
		  transp += devreconImageshort[vnum] * devAttenuation[mat*NE + ke] * rr * ph;
		}
	    }
	  xb0=xb1;
	  yb0=yb1;
	  zb0=zb1;
	  ik++;
	} //while
      
      //double intens = devpdf_pe[ke] * exp(-transp/xfactor);

      double intens = devpdf_pe[ke] * exp(-transp/xfactor)*leng;//2021/6/12

 

      //double intens =  MAX_AIR_PIXDIFF * devpdf_pe[ke] * exp(-transp/xfactor)*leng;//2022/6/20
 
      //double intens = MAX_PIXVALUE - (MAX_PIXVALUE * exp(-transp/xfactor)*leng);//2022/1/17
      //double intens = MAX_PIXVALUE - (MAX_PIXVALUE *devpdf_pe[ke] * exp(-transp/xfactor)*leng);//2022/1/17
      devprojVolume[jjjj] =  (float)(intens);
    } //dx dy

  
  /*

	    }//j
	}//jj
    }//jjj
  */

  
	    }//reprojectionkernel







	  
/*
cudaChannelFormatDesc descVGPU = cudaCreateChannelDesc(32, 0, 0, 0, cudaChannelFormatKindFloat);
cudaArray *cA_VolumeGPU;

cudaChannelFormatDesc desc1 = cudaCreateChannelDesc(32, 0, 0, 0, cudaChannelFormatKindFloat);
cudaChannelFormatDesc desc2 = cudaCreateChannelDesc(32, 0, 0, 0, cudaChannelFormatKindFloat);
cudaArray *cA_Matflo, *cA_xray;
*/

// hostCode

float*	devprojVolume;
unsigned short* devreconImageshort; 
float* devAttenuation;
float* devpdf_pe;

  
int dev=1;//set GPU number
cudaDeviceProp deviceProp;
cudaError_t err;

void
initializeGPU(int reconSize, double reconScale, int reconSlices, double reconStep,
	      int projImgHeight, int projImgWidth, int totalProjNumber,
	      double* cosTable, double* sinTable, double* XcenterShift, double* YcenterShift,
	      double x_reg, double y_reg, double z_reg,
	      double DETECTOR_PITCH_CM, double DETECTOR_PITCH_CM_AT_ISO_CENT,
	      double DIST_BTWN_SRC_AND_ISOCENT_CM, double DETECTOR_PIXEL_NUM,
	      int num_material)
{
  
  cudaSetDevice(dev);//set GPU number
  cudaMemcpyToSymbol( devcosTable, cosTable, sizeof(double)*totalProjNumber );
  cudaMemcpyToSymbol( devsinTable, sinTable, sizeof(double)*totalProjNumber );
  cudaMemcpyToSymbol( devXcenterShift, XcenterShift, sizeof(double)*totalProjNumber );
  cudaMemcpyToSymbol( devYcenterShift, YcenterShift, sizeof(double)*totalProjNumber );

  cudaMemcpyToSymbol( devreconSize, &reconSize, sizeof(int) );
  cudaMemcpyToSymbol( devreconScale, &reconScale, sizeof(double) );
  cudaMemcpyToSymbol( devreconSlices, &reconSlices, sizeof(int) );
  cudaMemcpyToSymbol( devreconStep, &reconStep, sizeof(double) );

  cudaMemcpyToSymbol( devx_reg, &x_reg, sizeof(double) );
  cudaMemcpyToSymbol( devy_reg, &y_reg, sizeof(double) );
  cudaMemcpyToSymbol( devz_reg, &z_reg, sizeof(double) );


  cudaMemcpyToSymbol( devDIST_BTWN_SRC_AND_ISOCENT_CM, &DIST_BTWN_SRC_AND_ISOCENT_CM, sizeof(double) );
  cudaMemcpyToSymbol( devDETECTOR_PITCH_CM_AT_ISO_CENT, &DETECTOR_PITCH_CM_AT_ISO_CENT, sizeof(double) );
  cudaMemcpyToSymbol( devDETECTOR_PIXEL_NUM, &DETECTOR_PIXEL_NUM, sizeof(double) );
  cudaMemcpyToSymbol( devDETECTOR_PITCH_CM, &DETECTOR_PITCH_CM, sizeof(double) );

  /*
  cudaMemcpyToSymbol( devDIST_BTWN_SRC_AND_ISOCENT_CM, &DIST_BTWN_SRC_AND_ISOCENT_CM, sizeof(int) );0  cudaMemcpyToSymbol( devDETECTOR_PITCH_CM_AT_ISO_CENT, &DETECTOR_PITCH_CM_AT_ISO_CENT, sizeof(int) );
  cudaMemcpyToSymbol( devDETECTOR_PIXEL_NUM, &DETECTOR_PIXEL_NUM, sizeof(int) );
  cudaMemcpyToSymbol( devDETECTOR_PITCH_CM, &DETECTOR_PITCH_CM, sizeof(int) );
  */
  cudaMemcpyToSymbol( devprojImgWidth, &projImgWidth, sizeof(int) );
  cudaMemcpyToSymbol( devprojImgHeight, &projImgHeight, sizeof(int) );

  CHECK(cudaMalloc( (void**)&devprojVolume, sizeof(float)* totalProjNumber*projImgHeight*projImgWidth ), err);
  if (err != cudaSuccess){
    printf("devprojVolume cudaMalloc error\n");
    exit(0);
  }

  printf("reconscale%lf¥n",reconScale);
  printf("reconstep%lf¥n",reconStep);
  printf("reconslices%d¥n",reconSlices);
  
  /*
  CHECK(cudaMalloc( (void**)&devreconImageshort, sizeof(unsigned short)*reconSize*reconSize*num_material ), err);
  if (err != cudaSuccess){
    printf("devreconImageshort cudaMalloc error\n");
    exit(0);
  }
  */
  
  //2021/5/19
   CHECK(cudaMalloc( (void**)&devreconImageshort, sizeof(unsigned short)*reconSize*reconSize*reconSlices*num_material ), err);
  if (err != cudaSuccess){
    printf("devreconImageshort cudaMalloc error\n");
    exit(0);
  }
  CHECK(cudaMalloc( (void**)&devpdf_pe, sizeof(float)*1024 ), err);
  if (err != cudaSuccess){
    printf("devpdf_pe cudaMalloc error\n");
    exit(0);
  }
    CHECK(cudaMalloc( (void**)&devAttenuation, sizeof(float)*1024*num_material ), err);
  if (err != cudaSuccess){
    printf("devAttenuation cudaMalloc error\n");
    exit(0);
  }
  /*
  ////////////////preparations of texture for projVolume////////////////  
  CHECK(cudaMallocArray(&cA_VolumeGPU,&descVGPU,(End_Channel-Start_Channel),projImgHeight), err);
  if (err != cudaSuccess){
    printf("cA_VolumeGPU cudaMalloc error\n");
    exit(0);
  }

  //Allocate a cudaArray with cudaMalloc3Drray()
  CHECK(cudaMalloc3DArray(&cA_Matflo,&desc1,make_cudaExtent(projImgHeight,reconSize,reconSize),0), err);
  if (err != cudaSuccess){
    printf("cA_Matflo cudaMalloc3DArray error\n");
    exit(0);
  }
  CHECK(cudaMalloc3DArray(&cA_xray,&desc2,make_cudaExtent(projImgHeight,reconSize,reconSize),0), err);
  if (err != cudaSuccess){
    printf("cA_xray cudaMalloc3DArray error\n");
    exit(0);
  }
  
  //////////page-lock of projVolume and npf_float///////////////
  cudaHostRegister(projVolume,projImgWidth*projImgHeight*sizeof(float),cudaHostRegisterPortable);
  //cudaHostRegister(projVolumeGPU,(End_Channel-Start_Channel)*projImgHeight*sizeof(float),cudaHostRegisterPortable); // 20171126
  cudaHostRegister(npf_float,reconSize*reconSize*sizeof(float),cudaHostRegisterPortable);
  cudaHostRegister(npf_float_re,reconSize*reconSize*sizeof(float),cudaHostRegisterPortable);


  cudaGetDeviceProperties(&deviceProp, dev);
  std::cout << "Device " << dev << std::endl;
  std::cout << "Maximum sizes of each dimension of a block: " << deviceProp.maxThreadsDim[0] << " " << deviceProp.maxThreadsDim[1] << " " << deviceProp.maxThreadsDim[2] << std::endl;
  std::cout << "Maximum sizes of each dimension of a grid: " << deviceProp.maxGridSize[0] << " " << deviceProp.maxGridSize[1] << " " << deviceProp.maxGridSize[2] << std::endl;
  */
}

void
terminateGPU()
{
  cudaFree(devprojVolume);
  cudaFree(devreconImageshort);
  cudaFree(devpdf_pe);
  cudaFree(devAttenuation);
}



void	
reprojectionOnGPU( int ite, int reconSize, int reconSlices, unsigned short* reconImageshort, float* reprojection_float, float* npf_double, float* projVolume,
		      double reconScale, double reconStep,
		      int projImgWidth, int projImgHeight, int startProjNumber, int endProjNumber,
		      double* anglesRad, double* cosTable, double* sinTable, double* XcenterShift, double* YcenterShift, double xfactor, double* xalpha,
		      double DETECTOR_PITCH_CM, double DETECTOR_PITCH_CM_AT_ISO_CENT,
		      double DIST_BTWN_SRC_AND_ISOCENT_CM, double DETECTOR_PIXEL_NUM,
		   double x_reg, double y_reg, double z_reg, float* Attenuation, float* pdf_pe, int num_mat, int NE, char* outputname)
{

  //cudaSetDevice(1);
  /*
  CHECK(cudaMemcpy( devreconImageshort, (void*)&(reconImageshort[0]), sizeof(unsigned short)*reconSize*reconSize*num_mat, cudaMemcpyHostToDevice ), err);
  if (err != cudaSuccess){
    printf("cudaMemcpyHostToDevice  error\n");
    exit(0);
  }
  */

  /*
  unsigned short*   reconImageshort2 = (unsigned short*)malloc( reconSize*reconSize* reconSlices*num_mat*sizeof(unsigned short) );
  
  //2021/5/19 //reverse array
  for(int i=0;i<reconSize*reconSize* reconSlices*num_mat;i++)
  {
    reconImageshort2[reconSize*reconSize* reconSlices*num_mat-1-i]=reconImageshort[i];
  }

  */
  CHECK(cudaMemcpy( devreconImageshort, (void*)&(reconImageshort[0]), sizeof(unsigned short)*reconSize*reconSize* reconSlices*num_mat, cudaMemcpyHostToDevice ), err);
  
  if (err != cudaSuccess){
    printf("cudaMemcpyHostToDevice  error\n");
    exit(0);
  }
    CHECK(cudaMemcpy( devpdf_pe, (void*)&(pdf_pe[0]), sizeof(float)*1024, cudaMemcpyHostToDevice ), err);
  if (err != cudaSuccess){
    printf("cudaMemcpyHostToDevice  error\n");
    exit(0);
  }
  CHECK(cudaMemcpy( devAttenuation, (void*)&(Attenuation[0]), sizeof(float)*1024*num_mat, cudaMemcpyHostToDevice ), err);
  if (err != cudaSuccess){
    printf("cudaMemcpyHostToDevice  error\n");
    exit(0);
  }
  


    float sum_pdf = 0.0;
    for (int i=0; i<NE; i++)
      {
	sum_pdf += pdf_pe[i];
	printf("%d %f sum:%f \n",i, pdf_pe[i], sum_pdf);
      }
    
   

    /* 
      //Haga
   
    int dimGrid_x = endProjNumber;
    int dimBlk_x=projImgWidth;
    int dimBlk_y=projImgHeight;

    dim3 dimBlk(dimBlk_x,dimBlk_y);
    dim3 dimGrid(dimGrid_x);
    */


    /*
    //this is slow mode 20minutes//36proj
    int dimGrid_y =projImgHeight ;//block
    int dimGrid_x=  projImgWidth;//block
    int dimBlk_x=endProjNumber;   //thread
    dim3 dimGrid(dimGrid_x,dimGrid_y);
    dim3 dimBlk(dimBlk_x);
    */


    /*
    ///  // 6min18sec
    int dimGrid_z = endProjNumber;//block
    int dimGrid_y=projImgWidth;//block
    int dimGrid_x=projImgHeight;//block
    int dimBlk_x=32;   //thread   Geforce thread number is 32
    int dimBlk_y=1;   //thread
    int dimBlk_z=1;   //thread
    dim3 dimGrid(dimGrid_x/dimBlk_x,dimGrid_y,dimGrid_z);//block
    dim3 dimBlk(dimBlk_x,dimBlk_y,dimBlk_z);//thread
    //
    */


    /*
    ///  //512512 detector  6min17sec
    int dimGrid_z = endProjNumber;//block
    int dimGrid_y=projImgWidth;//block
    int dimGrid_x=projImgHeight;//block
    int dimBlk_x=64;   //thread  for 512  Geforce thread number is 32
    int dimBlk_y=1;   //thread
    int dimBlk_z=1;   //thread
    dim3 dimGrid(dimGrid_x/dimBlk_x,dimGrid_y,dimGrid_z);//block
    dim3 dimBlk(dimBlk_x,dimBlk_y,dimBlk_z);//thread
    */

  

    ///512 dete 
    ///////////////////////////////////

    //400 dete
    
    /*
    /// 400400detector  5min 48sec   2021/5/25
    int dimGrid_z = endProjNumber;//block
    int dimGrid_y=projImgWidth;//block
    int dimGrid_x=projImgHeight;//block
    int dimBlk_x=40;   //thread   Geforce thread number is 32
    int dimBlk_y=1;   //thread
    int dimBlk_z=1;   //thread
    dim3 dimGrid(dimGrid_x/dimBlk_x,dimGrid_y,dimGrid_z);//block
    dim3 dimBlk(dimBlk_x,dimBlk_y,dimBlk_z);//thread
*/

    /*
    /// 400400detector  4min 39sec   2021/5/25
    int dimGrid_z = endProjNumber;//block
    int dimGrid_y=projImgWidth;//block
    int dimGrid_x=projImgHeight;//block
    int dimBlk_x=80;   //thread   Geforce thread number is 32
    int dimBlk_y=1;   //thread
    int dimBlk_z=1;   //thread
    dim3 dimGrid(dimGrid_x/dimBlk_x,dimGrid_y,dimGrid_z);//block
    dim3 dimBlk(dimBlk_x,dimBlk_y,dimBlk_z);//thread
    */

    
     /// 400400detector  4min 49sec   2021/5/25
    int dimGrid_z = endProjNumber;//block
    int dimGrid_y=projImgWidth;//block
    int dimGrid_x=projImgHeight;//block
    int dimBlk_x=200;   //thread   Geforce thread number is 32
    int dimBlk_y=1;   //thread
    int dimBlk_z=1;   //thread
    dim3 dimGrid(dimGrid_x/dimBlk_x,dimGrid_y,dimGrid_z);//block
    dim3 dimBlk(dimBlk_x,dimBlk_y,dimBlk_z);//thread
    
   
  


     
    
    //float*    projVolume_0 = (float*)new1DArray( projImgWidth*projImgHeight, UNIT_FLOAT32 );
     
    //float*    projVolume_0 = (float*)new1DArray( endProjNumber*projImgWidth*projImgHeight, UNIT_FLOAT32 );

    //float*    projVolume_1 = (float*)new1DArray( endProjNumber*projImgWidth*projImgHeight, UNIT_FLOAT32 );

    



    float*    projVolume_1 = (float*)malloc( endProjNumber*projImgWidth*projImgHeight*sizeof(float) );
    float*    projVolume_2 = (float*)malloc( endProjNumber*projImgWidth*projImgHeight*sizeof(float) );
    
     for(int ke = 0; ke < NE; ke++)
     {

    //for(int ke = 0; ke < 10; ke++)
    //{

       float*    projVolume_0 = (float*)malloc( endProjNumber*projImgWidth*projImgHeight*sizeof(float) );

    
       

	reprojectionKernel<<<dimGrid, dimBlk>>>(devprojVolume,
						devreconImageshort, xfactor,
						devAttenuation, devpdf_pe, num_mat, NE, ke);
	//printf("%d \n", ke);
	//CHECK(cudaMemcpy( (void*)&(reconImagefloat_before[0]), devreconImagefloat_before, sizeof(float)*reconSize*reconSize*num_mat, cudaMemcpyDeviceToHost ), err);
	//if (err != cudaSuccess){
	//  printf("devreconImagefloat_before cudaMemcpy DeviceToHost  error\n");
	//  exit(0);
	//}

	
	
	
	CHECK(cudaMemcpy( (void*)&(projVolume_0[0]), devprojVolume, sizeof(float)*endProjNumber*projImgWidth*projImgHeight, cudaMemcpyDeviceToHost ), err);
	if (err != cudaSuccess)
	  {
	  printf("devprojVolume cudaMemcpy DeviceToHost error\n");
	  exit(0);
	  }

       printf("ke%d/NE%d  %f %f %f %f \n", ke, NE,projVolume_0[0],projVolume_0[262144],projVolume_0[1],projVolume_0[130000]);
	
	
	for(int i = 0; i < endProjNumber*projImgWidth*projImgHeight; i++)
	  {
	    //projVolume[i] += projVolume_0[i];
	    projVolume_1[i] += projVolume_0[i];

	  
	  }


	if(ke==NE-2)//NE-2 100keV
	  {
	    for(int i = 0; i < endProjNumber*projImgWidth*projImgHeight; i++)
	      {
		projVolume_2[i] = projVolume_0[i];
		
	      }
	  }
	
     }//for ke

    
    FILE*	reproj;
    char filename[128];
    //sprintf(filename,"cuda_reprojection_polykeV_float.raw");
    //reproj=fopen(filename,"wb");
    sprintf(filename,"%s.raw",outputname);
    reproj=fopen(filename,"wb");

    /*
    fwrite( projVolume,
	    getByteSizeOfUnit(UNIT_FLOAT32),
	    endProjNumber*projImgWidth*projImgHeight,reproj);
    
    */
    
    
    fwrite( projVolume_1,
	    getByteSizeOfUnit(UNIT_FLOAT32),
	    endProjNumber*projImgWidth*projImgHeight,reproj);
   
    
    fclose(reproj);

    /*
     FILE*	reproj2;
    char filename2[128];
    sprintf(filename2,"100keV_cuda_reprojection_float.raw");
    reproj2=fopen(filename2,"wb");
    
    fwrite( projVolume_2,
	    getByteSizeOfUnit(UNIT_FLOAT32),
	    endProjNumber*projImgWidth*projImgHeight,reproj2);
   
    
    fclose(reproj2);
*/
    
//    cudaEventRecord(stop,0); //to here
//    cudaEventSynchronize(stop);
//    cudaEventElapsedTime(&elapsedTime,start,stop);
//    printf( "time: %8.2f ms\n", elapsedTime );

}



