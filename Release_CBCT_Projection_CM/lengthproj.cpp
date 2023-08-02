#include <stdio.h>
#include <string.h>
#include <math.h>

#include "physParams.h"
#include "mallocMD.h"
#include "reprojection.h"

void
lengthproj( int ite, int reconSize, int reconSlices, short* reconImageshort, unsigned short* ProjImageshort, float* npf_float, unsigned short* projVolume,
	      double reconScale, double reconStep,
	      int projImgWidth, int projImgHeight, int startProjNumber, int endProjNumber,
	      double* ProjAngle, double* cosTable, double* sinTable, double* XcenterShift, double* YcenterShift, double xfactor, double* xalpha,
	      double DETECTOR_PITCH_CM, double DETECTOR_PITCH_CM_AT_ISO_CENT,
	      double DIST_BTWN_SRC_AND_ISOCENT_CM, double DETECTOR_PIXEL_NUM,
	      double x_reg, double y_reg, double z_reg )
{
  if(ite != 0) 
    {
      x_reg = 0.0;
      y_reg = 0.0;
      z_reg = 0.0;
    }

  double rr;

  for(int i = 0; i < reconSlices*reconSize*reconSize; i++)
    {
      npf_float[i] = 0.0;
    }
#ifdef _OPENMP					// 20150828	OPENMP
#pragma omp parallel for private(rr)		// 20150828	OPENMP
#endif						// 20150828	OPENMP
  for(int jz=0;jz<reconSlices;jz++)
    {
      double rz = reconSlices*reconStep*0.5 - (jz+0.5)*reconStep;
      for(int jy=0;jy<reconSize;jy++)
	{
	  double ry = reconSize*reconScale*0.5 - (jy+0.5)*reconScale;
	  for(int jx=0;jx<reconSize;jx++)
	    {
	      double rx = -reconSize*reconScale*0.5 + (jx+0.5)*reconScale;

	      for(int jjj=startProjNumber;jjj<endProjNumber;jjj++) // projection angle;
		{
		  // Source location in cartesian coodinate
		  double sx = -DIST_BTWN_SRC_AND_ISOCENT_CM*sinTable[jjj] + x_reg;
		  double sy = DIST_BTWN_SRC_AND_ISOCENT_CM*cosTable[jjj] + y_reg;
		  double sz = z_reg;
		  
		  double zslope = fabs(rz-sz);
		  double yslope = fabs(ry-sy);
		  double xslope = fabs(rx-sx);

		  if(xslope > yslope)
		    {
		      double xx0 = rx - 0.5*reconScale;
		      double yy0 = (xx0-sx)/(rx-sx)*(ry-sy)+sy;
		      double zz0 = (xx0-sx)/(rx-sx)*(rz-sz)+sz;
		      double xx1 = rx + 0.5*reconScale;
		      double yy1 = (xx1-sx)/(rx-sx)*(ry-sy)+sy;
		      double zz1 = (xx1-sx)/(rx-sx)*(rz-sz)+sz;
		      rr = sqrt((xx0-xx1)*(xx0-xx1)+(yy0-yy1)*(yy0-yy1)+(zz0-zz1)*(zz0-zz1));
		    }
		  else
		    {
		      double yy0 = ry - 0.5*reconScale;
		      double xx0 = (yy0-sy)/(ry-sy)*(rx-sx)+sx;
		      double zz0 = (yy0-sy)/(ry-sy)*(rz-sz)+sz;
		      double yy1 = ry + 0.5*reconScale;
		      double xx1 = (yy1-sy)/(ry-sy)*(rx-sx)+sx;
		      double zz1 = (yy1-sy)/(ry-sy)*(rz-sz)+sz;
		      rr = sqrt((xx0-xx1)*(xx0-xx1)+(yy0-yy1)*(yy0-yy1)+(zz0-zz1)*(zz0-zz1));
		    }

		  double ox = sx+DETECTOR_PITCH_CM/DETECTOR_PITCH_CM_AT_ISO_CENT*sinTable[jjj]*DIST_BTWN_SRC_AND_ISOCENT_CM;
		  double oy = sy-DETECTOR_PITCH_CM/DETECTOR_PITCH_CM_AT_ISO_CENT*cosTable[jjj]*DIST_BTWN_SRC_AND_ISOCENT_CM;
		  double oz = sz;
		  double py = ( ox - rx + (sx-rx)/(sy-ry)*ry + (oy-sy)/(ox-sx)*oy )/((sx-rx)/(sy-ry) + (oy-sy)/(ox-sx));
		  double px = ((sx-rx)/(sy-ry)*(py-ry)+rx);
		  double pz = ((sz-rz)/(sy-ry)*(py-ry)+rz);
		  double xshift = -(XcenterShift[jjj] + projImgHeight/2.0)*(DETECTOR_PITCH_CM * (DETECTOR_PIXEL_NUM / projImgHeight)); // panel shift in [CM]
		  double yshift = -(YcenterShift[jjj] + projImgHeight/2.0)*(DETECTOR_PITCH_CM * (DETECTOR_PIXEL_NUM / projImgHeight));  // panel shift in [cm]
		  //if(jjj == 0 && jx == 0 && jy == 0) fprintf(stderr,"%f %f %f %f %f %f %d %d %d %lf %lf %lf %lf\n", sx, sy, sz, rx, ry, rz, jx, jy, jz, px, py, pz, rr);
		  //if(jjj == 0 && jx == 0 && jy == 0) fprintf(stderr,"%f %f %f %f %f %f %d %d %d %lf %lf %lf\n", sx, sy, sz, rx, ry, rz, jx, jy, jz, px, py, pz);
		  // In misc.cpp, 
		  //(*XcenterShift)[i] = -projImgHeight/2.0 - xshift[i]*10.0 / (DETECTOR_PITCH_CM * (DETECTOR_PIXEL_NUM / projImgHeight));
		  //(*YcenterShift)[i] = -projImgHeight/2.0 - yshift[i]*10.0 / (DETECTOR_PITCH_CM * (DETECTOR_PIXEL_NUM / projImgHeight));
		  // In reprojection.cpp
		  //double pr = DETECTOR_PITCH_CM_AT_ISO_CENT * (XcenterShift[jjj] + jj) * (DETECTOR_PIXEL_NUM / projImgHeight);
		  //double pz =  -(YcenterShift[jjj] + j) * DETECTOR_PITCH_CM_AT_ISO_CENT * (DETECTOR_PIXEL_NUM / projImgHeight) + z_reg*DETECTOR_PITCH_CM_AT_ISO_CENT/DETECTOR_PITCH_CM;
		  //double vv = -(pz-z_reg + yshift);
		  double vv = -(pz-z_reg - yshift);
		  double uu = (px-x_reg)*cosTable[jjj]+(py-x_reg)*sinTable[jjj] + xshift - (projImgWidth-projImgHeight)/2*(DETECTOR_PITCH_CM * (DETECTOR_PIXEL_NUM / projImgHeight));
		  double vfac = 0.5;
		  double ufac = 0.5;
		  if(vv < 0) vfac = -0.5;
		  if(uu < 0) ufac = -0.5;
		  int v1 = (int)(vv/(DETECTOR_PITCH_CM*(DETECTOR_PIXEL_NUM / projImgHeight))+vfac);
		  int u1 = (int)(uu/(DETECTOR_PITCH_CM*(DETECTOR_PIXEL_NUM / projImgHeight))+ufac);
		  double dv1 = vv/(DETECTOR_PITCH_CM*(DETECTOR_PIXEL_NUM / projImgHeight))-v1;
		  double du1 = uu/(DETECTOR_PITCH_CM*(DETECTOR_PIXEL_NUM / projImgHeight))-u1;
		  int v2 = v1 + int((vv-v1)/fabs(vv-v1));
		  int u2 = u1 + int((uu-u1)/fabs(uu-u1));
		  if(fabs(vv-v1)<0.000000001) v2 = v1;
		  if(fabs(uu-u1)<0.000000001) u2 = u1;
		  if(v1 + projImgHeight/2 >= 0 && v1 + projImgHeight/2 < projImgHeight && v2 + projImgHeight/2 >= 0 && v2 + projImgHeight/2 < projImgHeight
                     && u1 + projImgWidth/2 >= 0 && u1 + projImgWidth/2 < projImgWidth && u2 + projImgWidth/2 >= 0 && u2 + projImgWidth/2 < projImgWidth        )
		    {
		      int j11 = jjj*projImgWidth*projImgHeight + (v1 + projImgHeight/2)*projImgWidth + u1 + projImgWidth/2;
		      int j12 = jjj*projImgWidth*projImgHeight + (v2 + projImgHeight/2)*projImgWidth + u1 + projImgWidth/2;
		      int j21 = jjj*projImgWidth*projImgHeight + (v1 + projImgHeight/2)*projImgWidth + u2 + projImgWidth/2;
		      int j22 = jjj*projImgWidth*projImgHeight + (v2 + projImgHeight/2)*projImgWidth + u2 + projImgWidth/2;
		      double np_para1 = ((1-dv1)*projVolume[j11]+dv1*projVolume[j21]);
		      double np_para2 = ((1-dv1)*projVolume[j12]+dv1*projVolume[j22]);
		      double np_para3 = ((1-du1)*np_para1+du1*np_para2);
		      double np_para = (MAX_PIXVALUE - np_para3);
		      //double np_para = (MAX_PIXVALUE - projVolume[j11]);
		      //double np_para = 1.0;
		      int vnum = jx + jy*reconSize + jz*reconSize*reconSize;
		      npf_float[vnum] += (float)(np_para * rr);
		      //if(jx == 205 && jy == 205 && jz == 59) fprintf(stderr,"%f %f %f %f %f %f %d %d %d %f\n", sx, sy, sz, rx, ry, rz, jx, jy, jz, npf_float[vnum]);
		      //if(np_para * rr < 0) fprintf(stderr,"%f %f %f %f %f %f %d %d %d %f %lf %lf\n", sx, sy, sz, rx, ry, rz, jx, jy, jz, npf_float[vnum], np_para1, rr);
		    }
		} // projection angle

	    }
	}
    }
  return;
}

