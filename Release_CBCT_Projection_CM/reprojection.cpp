#include <stdio.h>
#include <string.h>
#include <math.h>

#include "physParams.h"
#include "mallocMD.h"
#include "reprojection.h"
//
void
reprojection( int ite, int reconSize, int reconSlices, unsigned short* reconImageshort, float* reprojection_float, float* npf_float, float* projVolume,
	      double reconScale, double reconStep,
	      int projImgWidth, int projImgHeight, int startProjNumber, int endProjNumber,
	      double* ProjAngle, double* cosTable, double* sinTable, double* XcenterShift, double* YcenterShift, double xfactor, double* xalpha,
	      double DETECTOR_PITCH_MM, double DETECTOR_PITCH_MM_AT_ISO_CENT,
	      double DIST_BTWN_SRC_AND_ISOCENT_MM, double DETECTOR_PIXEL_NUM,
	      double x_reg, double y_reg, double z_reg, float* Attenuation, float* pdf_pe, int num_mat, int NE, int ke )
{//////
  double gf=0, ff=0;
  if(ite != 0) 
    {//
      x_reg = 0.0;
      y_reg = 0.0;
      z_reg = 0.0;
    }
  int nxm = reconSize*0.5;
  int nym = reconSize*0.5;
  int nzm = reconSlices*0.5;
  int pri_size = reconSize*reconSize;
  //double actual_recon_z = (DIST_BTWN_SRC_AND_ISOCENT_MM_kV*DETECTOR_PITCH_MM_kV/DETECTOR_PITCH_MM_AT_ISO_CENT_kV)/(DIST_BTWN_SRC_AND_ISOCENT_MM_kV+reconSize*0.5*reconScale)*reconSlices*0.5*reconStep;
  //int projection_range = (int)(actual_recon_z/DETECTOR_PITCH_MM_kV) + 2;
  //double reslice = reconSlices*0.5*reconStep*(DIST_BTWN_SRC_AND_ISOCENT_MM_kV-reconSize*0.5*reconScale)/(DIST_BTWN_SRC_AND_ISOCENT_MM_kV+reconSize*0.5*reconScale);
  //double z_panel_range = (DIST_BTWN_SRC_AND_ISOCENT_MM_kV*DETECTOR_PITCH_MM_kV/DETECTOR_PITCH_MM_AT_ISO_CENT_kV)/(DIST_BTWN_SRC_AND_ISOCENT_MM_kV-reconSize*0.5*reconScale)*reslice;
  //int projection_range = (int)(z_panel_range/(DETECTOR_PITCH_MM_kV*(DETECTOR_PIXEL_NUM / projImgHeight))) + 15;
  fprintf(stderr, "projection range = %lf\n", z_reg);
  //for(int i = 0; i < reconSlices*reconSize*reconSize; i++)
  //  {
  //    npf_float[i] = 0.0;
  //  }
  //
  double rmax = reconSize*0.5 *reconScale;
#ifdef _OPENMP				// 20150828	OPENMP
#pragma omp parallel for		// 20150828	OPENMP
#endif
  //for(int jjj=startProjNumber;jjj<endProjNumber;jjj++)
  for(int jjj=0;jjj<1;jjj++)
  {
    // Source location in cartesian coodinate
    double sx = -DIST_BTWN_SRC_AND_ISOCENT_MM*sinTable[jjj] + x_reg;
    double sy = DIST_BTWN_SRC_AND_ISOCENT_MM*cosTable[jjj] + y_reg;
    double sz = z_reg;
    //
    for(int j=0;j<projImgHeight;j++) // Hight ////important!!!!!!
    //for(int j=projImgHeight/2-projection_range +1 ;j<projImgHeight/2+projection_range +1 ;j++) // Hight
    //for(int j=232;j<233;j++) // Hight
    {

      double pz =  -(YcenterShift[jjj] + j) * DETECTOR_PITCH_MM_AT_ISO_CENT * (DETECTOR_PIXEL_NUM / projImgHeight) + z_reg*DETECTOR_PITCH_MM_AT_ISO_CENT/DETECTOR_PITCH_MM;

      for(int jj=0;jj<projImgWidth;jj++)	// Width
      //for(int jj=359;jj<361;jj++)	        // Width
      {
      int jjjj = jjj*projImgWidth*projImgHeight + j*projImgWidth + jj;
	// Detector location at IC
      //double pr = DETECTOR_PITCH_MM_AT_ISO_CENT * (XcenterShift[jjj] + jj) * (DETECTOR_PIXEL_NUM / projImgWidth);
	double pr = DETECTOR_PITCH_MM_AT_ISO_CENT * (XcenterShift[jjj] + jj) * (DETECTOR_PIXEL_NUM / projImgHeight);
	double px = pr*cosTable[jjj] + x_reg*DETECTOR_PITCH_MM_AT_ISO_CENT/DETECTOR_PITCH_MM;
	double py = pr*sinTable[jjj] + y_reg*DETECTOR_PITCH_MM_AT_ISO_CENT/DETECTOR_PITCH_MM;
	double leng = DIST_BTWN_SRC_AND_ISOCENT_MM*DIST_BTWN_SRC_AND_ISOCENT_MM/(DIST_BTWN_SRC_AND_ISOCENT_MM*DIST_BTWN_SRC_AND_ISOCENT_MM + (pr*pr+pz*pz));
	leng = leng*sqrt(leng);
	//double leng = 1.0;
	//if(jjj == 0) fprintf(stderr,"%lf %lf %lf %lf %lf %lf %lf %lf %lf\n", ProjAngle[jjj]*180/PI, sx, sy, sz, px, py, pz, YcenterShift[jjj], XcenterShift[jjj] );
	//double np_para = (MAX_PIXVALUE - projVolume[jjjj]) * log(MAX_AIR_PIXDIFF/(MAX_PIXVALUE - projVolume[jjjj]));
	double np_para = (MAX_PIXVALUE - projVolume[jjjj]);
	//np_para = 1.0;
	double xx[4], yy[4], zz[4];
	int m = 0;
	yy[m] = double(nym)*reconScale;
	if(fabs((yy[m]-sy)/(py-sy)*(px-sx)+sx) < double(nxm)*reconScale) 
	{
	  xx[m] = (yy[m]-sy)/(py-sy)*(px-sx)+sx;
	  zz[m] = (yy[m]-sy)/(py-sy)*(pz-sz)+sz;
	  m++;
	}
	yy[m] = -double(nym)*reconScale;
	if(fabs((yy[m]-sy)/(py-sy)*(px-sx)+sx) < double(nxm)*reconScale)
	{
	  xx[m] = (yy[m]-sy)/(py-sy)*(px-sx)+sx;
	  zz[m] = (yy[m]-sy)/(py-sy)*(pz-sz)+sz;
	  m++;
	}
	xx[m] = double(nxm)*reconScale;
	if(fabs((xx[m]-sx)/(px-sx)*(py-sy)+sy) < double(nym)*reconScale) 
	{
	  yy[m] = (xx[m]-sx)/(px-sx)*(py-sy)+sy;
	  zz[m] = (xx[m]-sx)/(px-sx)*(pz-sz)+sz;
	  m++;
	}
	xx[m] = -double(nxm)*reconScale;
	if(fabs((xx[m]-sx)/(px-sx)*(py-sy)+sy) < double(nym)*reconScale)
	{
	  yy[m] = (xx[m]-sx)/(px-sx)*(py-sy)+sy;
	  zz[m] = (xx[m]-sx)/(px-sx)*(pz-sz)+sz;
	  m++;
	}
	double xi, yi, zi, xf, yf, zf; // Initial and Final Voxels, in which X-ray passes through
	if((xx[0]-sx)*(xx[0]-sx)+(yy[0]-sy)*(yy[0]-sy)+(zz[0]-sz)*(zz[0]-sz) < (xx[1]-sx)*(xx[1]-sx)+(yy[1]-sy)*(yy[1]-sy)+(zz[1]-sz)*(zz[1]-sz))
	{
	  xi=xx[0]; yi=yy[0]; zi=zz[0]; xf=xx[1]; yf=yy[1]; zf=zz[1];
	}
	else
	{
	  xi=xx[1]; yi=yy[1]; zi=zz[1]; xf=xx[0]; yf=yy[0]; zf=zz[0];
	}
	//if(jj == 255 | jj == 254) fprintf(stderr,"%d %d %d %lf %lf %lf %lf %lf %lf\n", jj, j, jjj, xi, yi, zi, xf, yf, zf);
	double xb0 = xi, yb0 = yi, zb0 = zi, xb1, yb1, zb1;
	double xii, yii, zii, xiii, yiii, ziii, xiiii, yiiii, ziiii, yj, xj, zj, rr, transp, total_length;
	int xma, yma, zma, xmb, ymb, zmb, xmc, ymc, zmc, icheck, vnum;
	//if(fabs(xi-xf) > fabs(yi-yf) && fabs(zf) < reconSlices*reconStep/2 ) // dx is larger than dy   *****************************************************


	
	//if(fabs(xi-xf) >= fabs(yi-yf) ) // dx is larger than dy   *****************************************************
	 if(fabs(xi-xf) >= fabs(yi-yf) ) // dx is larger than dy   ***************************************************** 2021/6/11 shimo
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
	  while( xii <= double(nxm)*reconScale && xii >= -double(nxm)*reconScale &&  yii <= double(nym)*reconScale && yii >= -double(nym)*reconScale)
	  {
	    xii = xin + pha * ik*reconScale;
	    ymb=yma;
	    yj = yii;
	    yii = (xii-sx)/(px-sx)*(py-sy)+sy;
	    yma=(int)(yii/reconScale);
	    zmb=zma;
	    zj = zii;
	    zii = (xii-sx)/(px-sx)*(pz-sz)+sz;
	    zma=(int)(zii/reconStep);
	    icheck = 0;
	    chy = 0;
	    chz = 0;
	    //if(jj==256 || jj==255) fprintf(stderr,"%e %lf \n", yii, yf-yi);
	    if((ik != 0 && yma != ymb) || (ik != 0 && (yii/fabs(yii) != yj/fabs(yj))))
	    {
	      if((yii < 0 && (yf-yi) < 0) || (yii > 0 && (yf-yi) > 0)) yiii = (double)(yma*reconScale);
	      else yiii = (double)(yma*reconScale) - (yf-yi)/fabs(yf-yi)*reconScale;
	      xiii = (yiii-sy)/(py-sy)*(px-sx)+sx;
	      ziii = (yiii-sy)/(py-sy)*(pz-sz)+sz;
	      chy = (xii-xiii)*(xii-xiii)+(yii-yiii)*(yii-yiii)+(zii-ziii)*(zii-ziii);
	      icheck++;
	      //if(jj==256 || jj==255) fprintf(stderr,"%lf %lf %lf\n", xiii, sx, (yiii-sy)/(py-sy));
	      //if(aaa < 40 && aaa > 70) fprintf(stderr,"%lf %lf %lf %lf\n", aaa, xiii, yiii, ziii);
            }
            //if((ik != 0 && zma != zmb) || (ik != 0 && (zii/fabs(zii) != zj/fabs(zj))))
	    if(ik != 0 && zma != zmb)   //2021/6/11
            {
	      if((zii < 0 && (zf-zi) < 0) || (zii > 0 && (zf-zi) > 0)) ziiii = (double)(zma*reconStep);
	      else ziiii = (double)(zma*reconStep) - (zf-zi)/fabs(zf-zi)*reconStep;
              //  ziiii = (double)zma;
	      xiiii = (ziiii-sz)/(pz-sz)*(px-sx)+sx;
	      yiiii = (ziiii-sz)/(pz-sz)*(py-sy)+sy;
	      chz = (xii-xiiii)*(xii-xiiii)+(yii-yiiii)*(yii-yiiii)+(zii-ziiii)*(zii-ziiii);
	      icheck++;
	      //fprintf(stderr,"%lf %lf %lf %lf\n", aaa, xiii, yiii, ziii);
	      //if(aaa < 40 && aaa > 70) fprintf(stderr,"%lf %lf %lf %lf\n", aaa, xiiii, yiiii, ziiii);
            }
            if(icheck == 2 && chz > chy)
            {
                iks++;
                xb1=xiiii;
                yb1=yiiii;
                zb1=ziiii;
                rr = sqrt((xb0-xb1)*(xb0-xb1)+(yb0-yb1)*(yb0-yb1)+(zb0-zb1)*(zb0-zb1));
                xmc = (int)((xb0+xb1)*0.5/reconScale+nxm);
                ymc = (int)((-yb0-yb1)*0.5/reconScale+nym);
                zmc = (int)((-zb0-zb1)*0.5/reconStep+nzm);
		//if(j==0 || j==1) fprintf(stderr,"%lf %d %d %d %d %d\n", rr,xmc,ymc,zmc,j,jj);/////////////////////
                if(xmc >= 0 && ymc >= 0 && zmc >= 0 && xmc <= reconSize-1 && ymc <= reconSize-1 && zmc <= reconSlices-1)
                {
 		    //double r = sqrt((xb0+xb1)*(xb0+xb1)+(-yb0-yb1)*(-yb0-yb1)+(-zb0-zb1)*(-zb0-zb1))*0.5;
		    double r = sqrt((xb0+xb1)*(xb0+xb1)+(-yb0-yb1)*(-yb0-yb1))*0.5;
		    double ph = 1.0;
		    for(int mat = 0; mat < num_mat; mat++)
		      {
			double ph = 1.0;
			vnum = xmc + ymc*reconSize + (zmc*num_mat+mat)*reconSize*reconSize;
			if(reconImageshort[vnum] < 0 || r >= rmax) ph = 0.0;
			transp += reconImageshort[vnum] * Attenuation[mat*NE + ke] * rr * ph;
		      }
		    //npf_float[vnum] += (float)(np_para * rr);
 		    total_length += rr;
		    //fprintf(stderr,"%d %d %d %d %lf %lf \n", xmc, ymc, zmc, vnum, npf_float[vnum], exp(-transp));
                }
		xb0=xb1;
		yb0=yb1;
		zb0=zb1;
		iks++;
		xb1=xiii;
		yb1=yiii;
		zb1=ziii;
		rr = sqrt((xb0-xb1)*(xb0-xb1)+(yb0-yb1)*(yb0-yb1)+(zb0-zb1)*(zb0-zb1));
                xmc = (int)((xb0+xb1)*0.5/reconScale+nxm);
                ymc = (int)((-yb0-yb1)*0.5/reconScale+nym);
                zmc = (int)((-zb0-zb1)*0.5/reconStep+nzm);
		//if(jj==256 || jj==255) fprintf(stderr,"%lf %d %d %d %d %d\n", rr,xmc,ymc,zmc,j,jj);/////////////////////
		if(xmc >= 0 && ymc >= 0 && zmc >= 0 && xmc <= reconSize-1 && ymc <= reconSize-1 && zmc <= reconSlices-1)
		{
 		  //double r = sqrt((xb0+xb1)*(xb0+xb1)+(-yb0-yb1)*(-yb0-yb1)+(-zb0-zb1)*(-zb0-zb1))*0.5;
		  double r = sqrt((xb0+xb1)*(xb0+xb1)+(-yb0-yb1)*(-yb0-yb1))*0.5;
		  double ph = 1.0;
		  for(int mat = 0; mat < num_mat; mat++)
		    {
		      double ph = 1.0;
		      vnum = xmc + ymc*reconSize + (zmc*num_mat+mat)*reconSize*reconSize;
		      if(reconImageshort[vnum] < 0 || r >= rmax) ph = 0.0;
		      transp += reconImageshort[vnum] * Attenuation[mat*NE + ke] * rr * ph;
		    }
		  //npf_float[vnum] += (float)(np_para * rr);
 		  total_length += rr;
		  //fprintf(stderr,"%d %d %d %d %lf %lf \n", xmc, ymc, zmc, vnum, npf_float[vnum], exp(-transp));
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
	      xmc = (int)((xb0+xb1)*0.5/reconScale+nxm);
	      ymc = (int)((-yb0-yb1)*0.5/reconScale+nym);
	      zmc = (int)((-zb0-zb1)*0.5/reconStep+nzm);
	      //if(jj==256 || jj==255) fprintf(stderr,"%lf %d %d %d %d %d \n", rr,xmc,ymc,zmc,j,jj);/////////////////////
	      if(xmc >= 0 && ymc >= 0 && zmc >= 0 && xmc <= reconSize-1 && ymc <= reconSize-1 && zmc <= reconSlices-1)
	      {

		//double r = sqrt((xb0+xb1)*(xb0+xb1)+(-yb0-yb1)*(-yb0-yb1)+(-zb0-zb1)*(-zb0-zb1))*0.5;
		double r = sqrt((xb0+xb1)*(xb0+xb1)+(-yb0-yb1)*(-yb0-yb1))*0.5;
		double ph = 1.0;
		for(int mat = 0; mat < num_mat; mat++)
		  {
		    double ph = 1.0;
		    vnum = xmc + ymc*reconSize + (zmc*num_mat+mat)*reconSize*reconSize;
		    if(reconImageshort[vnum] < 0 || r >= rmax) ph = 0.0;
		    transp += reconImageshort[vnum] * Attenuation[mat*NE + ke] * rr * ph;
		  }
		//npf_float[vnum] += (float)(np_para * rr);
		total_length += rr;
		//fprintf(stderr,"%d %d %d %d %lf %lf \n", xmc, ymc, zmc, vnum, npf_float[vnum], exp(-transp));
	      }
	      xb0=xb1;
	      yb0=yb1;
	      zb0=zb1;
	      iks++;
	      xb1=xiiii;
              yb1=yiiii;
              zb1=ziiii;
              rr = sqrt((xb0-xb1)*(xb0-xb1)+(yb0-yb1)*(yb0-yb1)+(zb0-zb1)*(zb0-zb1));
	      xmc = (int)((xb0+xb1)*0.5/reconScale+nxm);
	      ymc = (int)((-yb0-yb1)*0.5/reconScale+nym);
	      zmc = (int)((-zb0-zb1)*0.5/reconStep+nzm);
	      //if(jj==256 || jj==255)  fprintf(stderr,"%lf %d %d %d %d %d \n", rr,xmc,ymc,zmc,j,jj);/////////////////////
              if(xmc >= 0 && ymc >= 0 && zmc >= 0 && xmc <= reconSize-1 && ymc <= reconSize-1 && zmc <= reconSlices-1)
              {
    		      //double r = sqrt((xb0+xb1)*(xb0+xb1)+(-yb0-yb1)*(-yb0-yb1)+(-zb0-zb1)*(-zb0-zb1))*0.5;
		      double r = sqrt((xb0+xb1)*(xb0+xb1)+(-yb0-yb1)*(-yb0-yb1))*0.5;
		      double ph = 1.0;
		      for(int mat = 0; mat < num_mat; mat++)
			{
			  double ph = 1.0;
			  vnum = xmc + ymc*reconSize + (zmc*num_mat+mat)*reconSize*reconSize;
			  if(reconImageshort[vnum] < 0 || r >= rmax) ph = 0.0;
			  transp += reconImageshort[vnum] * Attenuation[mat*NE + ke] * rr * ph;
			}
    		      //npf_float[vnum] += (float)(np_para * rr);
    		      total_length += rr;
		      //fprintf(stderr,"%d %d %d %d %lf %lf \n", xmc, ymc, zmc, vnum, npf_float[vnum], exp(-transp));
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
		xmc = (int)((xb0+xb1)*0.5/reconScale+nxm);
		ymc = (int)((-yb0-yb1)*0.5/reconScale+nym);
		zmc = (int)((-zb0-zb1)*0.5/reconStep+nzm);
		//if(jj==256 || jj==255) fprintf(stderr,"%lf %d %d %d %d %d %lf %lf \n", rr,xmc,ymc,zmc,j,jj,xb0,xb1);/////////////////////
		if(xmc >= 0 && ymc >= 0 && zmc >= 0 && xmc <= reconSize-1 && ymc <= reconSize-1 && zmc <= reconSlices-1)
                  {
		    //double r = sqrt((xb0+xb1)*(xb0+xb1)+(-yb0-yb1)*(-yb0-yb1)+(-zb0-zb1)*(-zb0-zb1))*0.5;
		    double r = sqrt((xb0+xb1)*(xb0+xb1)+(-yb0-yb1)*(-yb0-yb1))*0.5;
		    double ph = 1.0;
		    for(int mat = 0; mat < num_mat; mat++)
		      {
			double ph = 1.0;
			vnum = xmc + ymc*reconSize + (zmc*num_mat+mat)*reconSize*reconSize;
			if(reconImageshort[vnum] < 0 || r >= rmax) ph = 0.0;
			transp += reconImageshort[vnum] * Attenuation[mat*NE + ke] * rr * ph;
		      }
		    //npf_float[vnum] += (float)(np_para * rr);
		    total_length += rr;
		    //fprintf(stderr,"%d %d %d %d %lf %lf \n", xmc, ymc, zmc, vnum, npf_float[vnum], exp(-transp));
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
		xmc = (int)((xb0+xb1)*0.5/reconScale+nxm);
		ymc = (int)((-yb0-yb1)*0.5/reconScale+nym);
		zmc = (int)((-zb0-zb1)*0.5/reconStep+nzm);
		//if(jj==256 || jj==255) fprintf(stderr,"%lf %d %d %d %d %d %lf %lf \n", rr,xmc,ymc,zmc,j,jj, xb0, xb1);/////////////////////
		if(xmc >= 0 && ymc >= 0 && zmc >= 0 && xmc <= reconSize-1 && ymc <= reconSize-1 && zmc <= reconSlices-1)
                  {
		    //double r = sqrt((xb0+xb1)*(xb0+xb1)+(-yb0-yb1)*(-yb0-yb1)+(-zb0-zb1)*(-zb0-zb1))*0.5;
		    double r = sqrt((xb0+xb1)*(xb0+xb1)+(-yb0-yb1)*(-yb0-yb1))*0.5;
		    for(int mat = 0; mat < num_mat; mat++)
		      {
			double ph = 1.0;
			vnum = xmc + ymc*reconSize + (zmc*num_mat+mat)*reconSize*reconSize;
			if(reconImageshort[vnum] < 0 || r >= rmax) ph = 0.0;
			transp += reconImageshort[vnum] * Attenuation[mat*NE + ke] * rr * ph;
		      }
		    //npf_float[vnum] += (float)(np_para * rr);
		    total_length += rr;
		    //if(jj==256 || jj==255) fprintf(stderr,"%d %d %d %d %lf %lf \n", xmc, ymc, zmc, vnum, npf_float[vnum], exp(-transp));
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
	    xmc = (int)((xb0+xb1)*0.5/reconScale+nxm);
	    ymc = (int)((-yb0-yb1)*0.5/reconScale+nym);
	    zmc = (int)((-zb0-zb1)*0.5/reconStep+nzm);
	    //if(j==0 || j==1) fprintf(stderr,"%lf %d %d %d %d %d \n", rr,xmc,ymc,zmc,j,jj);/////////////////////
	    if(xmc >= 0 && ymc >= 0 && zmc >= 0 && xmc <= reconSize-1 && ymc <= reconSize-1 && zmc <= reconSlices-1)
              {
		//double r = sqrt((xb0+xb1)*(xb0+xb1)+(-yb0-yb1)*(-yb0-yb1)+(-zb0-zb1)*(-zb0-zb1))*0.5;
		double r = sqrt((xb0+xb1)*(xb0+xb1)+(-yb0-yb1)*(-yb0-yb1))*0.5;
		for(int mat = 0; mat < num_mat; mat++)
		  {
		    double ph = 1.0;
		    vnum = xmc + ymc*reconSize + (zmc*num_mat+mat)*reconSize*reconSize;
		    if(reconImageshort[vnum] < 0 || r >= rmax) ph = 0.0;
		    transp += reconImageshort[vnum] * Attenuation[mat*NE + ke] * rr * ph;
		  }
		//npf_float[vnum] += (float)(np_para * rr);
		total_length += rr;
		//if(j==0 || j==1) fprintf(stderr,"%d %d %d %d %lf %lf \n", xmc, ymc, zmc, vnum, reconImageshort[vnum], rr);
		//fprintf(stderr,"%d %d %d %d %lf %lf \n", xmc, ymc, zmc, vnum, reconImageshort[vnum], transp);
              }
              //fprintf(stderr,"%lf %d %d %lf %lf %lf %lf *\n", aaa, iks, vnum, xii, yii, zii, total_length);
             
              xb0=xb1;
              yb0=yb1;
              zb0=zb1;
              ik++;
              //if(jj==256 || jj==255) fprintf(stderr,"%d %d %lf %lf %lf \n", iks, nzm, xb0, xb1, rr);
	      //fprintf(stderr,"%lf %d %d %d %d %d %d %lf\n", aaa, iks, vnum, xmc, ymc, zmc, reconImageshort[vnum], transp);
          } // while
    	  //double intens = MAX_PIXVALUE - (MAX_AIR_PIXDIFF * exp(-transp/xfactor)*leng); // * atten;
	  double intens = pdf_pe[ke] * exp(-transp/xfactor);
	  if(ite == 0 && intens > 57000)
	  //if(ite == 0)
	    {
	      gf += log(np_para/MAX_AIR_PIXDIFF/leng) * transp/xfactor;
	      ff += transp/xfactor * transp/xfactor;
	      //xalpha += -gf/ff; 
	    }
	  reprojection_float[jjjj] = (float)(intens);
	  //if(j == projImgHeight/2-projection_range || j == projImgHeight/2+projection_range-1) ProjImageshort[jjjj] = 25000;
	  //double length = sqrt((xi-xf)*(xi-xf)+(yi-yf)*(yi-yf)+(zi-zf)*(zi-zf));
	  //if(j==0 || j==1) fprintf(stderr,"%d %d %lf %d %lf %lf %lf %lf %lf %lf %e %e \n", j, jj, transp, ProjImageshort[jjjj], xi,yi,zi,xf,yf,zf, leng, total_length);
	  //fprintf(stderr,"%d %d %lf %lf %lf %lf\n", j, jj, intens, gf, ff, -gf/ff);
	  //if(aaa > 75 && aaa < 76) ProjImageshort[j][jj] = (unsigned short)(intens);
	  //if(aaa > 75 && aaa < 76) fprintf(stderr,"%lf %d %d %d %lf\n", aaa, j, jj, ProjImageshort[j][jj], atten);
        } // dx dy 


	

	// dy is larger than dx *********************************************************************
	// dy is larger than dx *********************************************************************
	// dy is larger than dx *********************************************************************
	// dy is larger than dx *********************************************************************
	//else if(fabs(xi-xf) < fabs(yi-yf) && fabs(zf) < reconSlices*reconStep/2) 

	//if(fabs(xi-xf) < fabs(yi-yf) )

	 else //2021/6/11
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
	  // fprintf(stderr,"%lf %lf %lf %lf\n", aaa, xii, yii, zii);
	  while( xii <= double(nxm)*reconScale && xii >= -double(nxm)*reconScale &&  yii <= double(nym)*reconScale && yii >= -double(nym)*reconScale)
	  {
	    yii = yin + pha * ik*reconScale;
	    xmb=xma;
	    xj = xii;
	    xii = (yii-sy)/(py-sy)*(px-sx)+sx;
	    xma=(int)(xii/reconScale);
	    zmb=zma;
	    zj = zii;
	    zii = (yii-sy)/(py-sy)*(pz-sz)+sz;
	    zma=(int)(zii/reconStep);
	    icheck = 0;
	    chy = 0;
	    chz = 0;
	    if((ik != 0 && xma != xmb) || (ik != 0 && (xii/fabs(xii) != xj/fabs(xj))))
	    {
	      if((xii < 0 && (xf-xi) < 0) || (xii > 0 && (xf-xi) > 0)) xiii = (double)(xma*reconScale);
	      else xiii = (double)(xma*reconScale) - (xf-xi)/fabs(xf-xi)*reconScale;
	      yiii = (xiii-sx)/(px-sx)*(py-sy)+sy;
	      ziii = (xiii-sx)/(px-sx)*(pz-sz)+sz;
	      chy = (xii-xiii)*(xii-xiii)+(yii-yiii)*(yii-yiii)+(zii-ziii)*(zii-ziii);
	      icheck++;
	      //if(aaa > 160 && aaa < 161) fprintf(stderr,"%lf %lf %lf %lf\n", aaa, xiii, yiii, ziii);
	    }
	    
	    //if((ik != 0 && zma != zmb) || (ik != 0 && (zii/fabs(zii) != zj/fabs(zj))))
	    if(ik != 0 && zma != zmb)   //2021/6/11
	    {
	      if((zii < 0 && (zf-zi) < 0) || (zii > 0 && (zf-zi) > 0)) ziiii = (double)(zma*reconStep);
	      else ziiii = (double)(zma*reconStep) - (zf-zi)/fabs(zf-zi)*reconStep;
	      //ziiii = (double)zma*reconStep;
	      xiiii = (ziiii-sz)/(pz-sz)*(px-sx)+sx;
	      yiiii = (ziiii-sz)/(pz-sz)*(py-sy)+sy;
	      chz = (xii-xiiii)*(xii-xiiii)+(yii-yiiii)*(yii-yiiii)+(zii-ziiii)*(zii-ziiii);
	      icheck++;
	      //if(aaa > 160 && aaa < 161) fprintf(stderr,"%lf %lf %lf %lf\n", aaa, xiiii, yiiii, ziiii);
	      //fprintf(stderr,"%lf %lf %lf %d %d\n", zii, zj, ziiii, zma, zmb);
	    }
 			  				  
	    if(icheck == 2 && chz > chy)
	    {
	      iks++;
	      xb1=xiiii;
	      yb1=yiiii;
	      zb1=ziiii;
	      rr = sqrt((xb0-xb1)*(xb0-xb1)+(yb0-yb1)*(yb0-yb1)+(zb0-zb1)*(zb0-zb1));
	      xmc = (int)((xb0+xb1)*0.5/reconScale+nxm);
	      ymc = (int)((-yb0-yb1)*0.5/reconScale+nym);
	      zmc = (int)((-zb0-zb1)*0.5/reconStep+nzm);
	      //fprintf(stderr,"%lf %d %d %d %d %d\n", rr,xmc,ymc,zmc,j,jj);/////////////////////
	      if(xmc >= 0 && ymc >= 0 && zmc >= 0 && xmc <= reconSize-1 && ymc <= reconSize-1 && zmc <= reconSlices-1)
	      {
		//double r = sqrt((xb0+xb1)*(xb0+xb1)+(-yb0-yb1)*(-yb0-yb1)+(-zb0-zb1)*(-zb0-zb1))*0.5;
		double r = sqrt((xb0+xb1)*(xb0+xb1)+(-yb0-yb1)*(-yb0-yb1))*0.5;
		for(int mat = 0; mat < num_mat; mat++)
		  {
		    double ph = 1.0;
		    vnum = xmc + ymc*reconSize + (zmc*num_mat+mat)*reconSize*reconSize;
		    if(reconImageshort[vnum] < 0 || r >= rmax) ph = 0.0;
		    transp += reconImageshort[vnum] * Attenuation[mat*NE + ke] * rr * ph;
		  }
		//npf_float[vnum] += (float)(np_para * rr);
		total_length += rr;
		//fprintf(stderr,"%d %d %d %d %lf %lf \n", xmc, ymc, zmc, vnum, npf_float[vnum], exp(-transp));
		//fprintf(stderr,"%d %d %d %d %d %lf %lf\n", xmc, ymc, zmc, vnum, reconImageshort[vnum], transp, xb0);
	      }
	      //if(aaa > 160 && aaa < 161) fprintf(stderr,"%lf %lf %lf %lf\n", aaa, xb1, yb1, zb1);
	      xb0=xb1;
	      yb0=yb1;
	      zb0=zb1;
	      iks++;
	      xb1=xiii;
	      yb1=yiii;
	      zb1=ziii;
	      rr = sqrt((xb0-xb1)*(xb0-xb1)+(yb0-yb1)*(yb0-yb1)+(zb0-zb1)*(zb0-zb1));
	      xmc = (int)((xb0+xb1)*0.5/reconScale+nxm);
	      ymc = (int)((-yb0-yb1)*0.5/reconScale+nym);
	      zmc = (int)((-zb0-zb1)*0.5/reconStep+nzm);
	      //fprintf(stderr,"%lf %d %d %d %d %d\n", rr,xmc,ymc,zmc,j,jj);/////////////////////
	      if(xmc >= 0 && ymc >= 0 && zmc >= 0 && xmc <= reconSize-1 && ymc <= reconSize-1 && zmc <= reconSlices-1)
	      {
		//double r = sqrt((xb0+xb1)*(xb0+xb1)+(-yb0-yb1)*(-yb0-yb1)+(-zb0-zb1)*(-zb0-zb1))*0.5;
		double r = sqrt((xb0+xb1)*(xb0+xb1)+(-yb0-yb1)*(-yb0-yb1))*0.5;
		for(int mat = 0; mat < num_mat; mat++)
		  {
		    double ph = 1.0;
		    vnum = xmc + ymc*reconSize + (zmc*num_mat+mat)*reconSize*reconSize;
		    if(reconImageshort[vnum] < 0 || r >= rmax) ph = 0.0;
		    transp += reconImageshort[vnum] * Attenuation[mat*NE + ke] * rr * ph;
		  }
		//npf_float[vnum] += (float)(np_para * rr);
		total_length += rr;
		//fprintf(stderr,"%d %d %d %d %lf %lf \n", xmc, ymc, zmc, vnum, npf_float[vnum], exp(-transp));
		//fprintf(stderr,"%d %d %d %d %d %lf %lf\n", xmc, ymc, zmc, vnum, reconImageshort[vnum], transp, xb0);
	      }
	      //if(aaa > 160 && aaa < 161) fprintf(stderr,"%lf %lf %lf %lf\n", aaa, xb1, yb1, zb1);
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
	      xmc = (int)((xb0+xb1)*0.5/reconScale+nxm);
	      ymc = (int)((-yb0-yb1)*0.5/reconScale+nym);
	      zmc = (int)((-zb0-zb1)*0.5/reconStep+nzm);
	      
	      if(xmc >= 0 && ymc >= 0 && zmc >= 0 && xmc <= reconSize-1 && ymc <= reconSize-1 && zmc <= reconSlices-1)
	      {
  		//double r = sqrt((xb0+xb1)*(xb0+xb1)+(-yb0-yb1)*(-yb0-yb1)+(-zb0-zb1)*(-zb0-zb1))*0.5;
		double r = sqrt((xb0+xb1)*(xb0+xb1)+(-yb0-yb1)*(-yb0-yb1))*0.5;
		for(int mat = 0; mat < num_mat; mat++)
		  {
		    double ph = 1.0;
		    vnum = xmc + ymc*reconSize + (zmc*num_mat+mat)*reconSize*reconSize;
		    if(reconImageshort[vnum] < 0 || r >= rmax) ph = 0.0;
		    transp += reconImageshort[vnum] * Attenuation[mat*NE + ke] * rr * ph;
		  }
		//npf_float[vnum] += (float)(np_para * rr);
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
	      xmc = (int)((xb0+xb1)*0.5/reconScale+nxm);
	      ymc = (int)((-yb0-yb1)*0.5/reconScale+nym);
	      zmc = (int)((-zb0-zb1)*0.5/reconStep+nzm);
	      //fprintf(stderr,"%lf %d %d %d %d %d\n", rr,xmc,ymc,zmc,j,jj);/////////////////////
	      if(xmc >= 0 && ymc >= 0 && zmc >= 0 && xmc <= reconSize-1 && ymc <= reconSize-1 && zmc <= reconSlices-1)
	      {
		//double r = sqrt((xb0+xb1)*(xb0+xb1)+(-yb0-yb1)*(-yb0-yb1)+(-zb0-zb1)*(-zb0-zb1))*0.5;
		double r = sqrt((xb0+xb1)*(xb0+xb1)+(-yb0-yb1)*(-yb0-yb1))*0.5;
		for(int mat = 0; mat < num_mat; mat++)
		  {
		    double ph = 1.0;
		    vnum = xmc + ymc*reconSize + (zmc*num_mat+mat)*reconSize*reconSize;
		    if(reconImageshort[vnum] < 0 || r >= rmax) ph = 0.0;
		    transp += reconImageshort[vnum] * Attenuation[mat*NE + ke] * rr * ph;
		  }
	        //npf_float[vnum] += (float)(np_para * rr);
		total_length += rr;
		//fprintf(stderr,"%d %d %d %d %lf %lf \n", xmc, ymc, zmc, vnum, npf_float[vnum], exp(-transp));
		//fprintf(stderr,"%d %d %d %d %d %lf %lf\n", xmc, ymc, zmc, vnum, reconImageshort[vnum], transp, xb0);
	      }
	      //if(aaa > 160 && aaa < 161) fprintf(stderr,"%lf %lf %lf %lf\n", aaa, xb1, yb1, zb1);
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
	      xmc = (int)((xb0+xb1)*0.5/reconScale+nxm);
	      ymc = (int)((-yb0-yb1)*0.5/reconScale+nym);
	      zmc = (int)((-zb0-zb1)*0.5/reconStep+nzm);
	      //fprintf(stderr,"%lf %d %d %d %d %d\n", rr,xmc,ymc,zmc,j,jj);/////////////////////
	      if(xmc >= 0 && ymc >= 0 && zmc >= 0 && xmc <= reconSize-1 && ymc <= reconSize-1 && zmc <= reconSlices-1)
	      {
		//double r = sqrt((xb0+xb1)*(xb0+xb1)+(-yb0-yb1)*(-yb0-yb1)+(-zb0-zb1)*(-zb0-zb1))*0.5;
		double r = sqrt((xb0+xb1)*(xb0+xb1)+(-yb0-yb1)*(-yb0-yb1))*0.5;
		for(int mat = 0; mat < num_mat; mat++)
		  {
		    double ph = 1.0;
		    vnum = xmc + ymc*reconSize + (zmc*num_mat+mat)*reconSize*reconSize;
		    if(reconImageshort[vnum] < 0 || r >= rmax) ph = 0.0;
		    transp += reconImageshort[vnum] * Attenuation[mat*NE + ke] * rr * ph;
		  }
		//npf_float[vnum] += (float)(np_para * rr);
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
	      xmc = (int)((xb0+xb1)*0.5/reconScale+nxm);
	      ymc = (int)((-yb0-yb1)*0.5/reconScale+nym);
	      zmc = (int)((-zb0-zb1)*0.5/reconStep+nzm);
	      //fprintf(stderr,"%lf %d %d %d %d %d\n", rr,xmc,ymc,zmc,j,jj);/////////////////////
	      if(xmc >= 0 && ymc >= 0 && zmc >= 0 && xmc <= reconSize-1 && ymc <= reconSize-1 && zmc <= reconSlices-1)
	      {
		//double r = sqrt((xb0+xb1)*(xb0+xb1)+(-yb0-yb1)*(-yb0-yb1)+(-zb0-zb1)*(-zb0-zb1))*0.5;
		double r = sqrt((xb0+xb1)*(xb0+xb1)+(-yb0-yb1)*(-yb0-yb1))*0.5;
		for(int mat = 0; mat < num_mat; mat++)
		  {
		    double ph = 1.0;
		    vnum = xmc + ymc*reconSize + (zmc*num_mat+mat)*reconSize*reconSize;
		    if(reconImageshort[vnum] < 0 || r >= rmax) ph = 0.0;
		    transp += reconImageshort[vnum] * Attenuation[mat*NE + ke] * rr * ph;
		  }
		//npf_float[vnum] += (float)(np_para * rr);
		total_length += rr;
	
	      }
	      //if(aaa > 160 && aaa < 161) fprintf(stderr,"%lf %lf %lf %lf\n", aaa, xb1, yb1, zb1);
	      xb0=xb1;
	      yb0=yb1;
	      zb0=zb1;
	    }
	    iks++;
	    xb1=xii;
	    yb1=yii;
	    zb1=zii;
	    rr = sqrt((xb0-xb1)*(xb0-xb1)+(yb0-yb1)*(yb0-yb1)+(zb0-zb1)*(zb0-zb1));
	    xmc = (int)((xb0+xb1)*0.5/reconScale+nxm);
	    ymc = (int)((-yb0-yb1)*0.5/reconScale+nym);
	    zmc = (int)((-zb0-zb1)*0.5/reconStep+nzm);
	    //if(jj==256 || jj==255) fprintf(stderr,"%lf %d %d %d %d %d\n", rr,xmc,ymc,zmc,j,jj);/////////////////////
	    //if(ik == 0) rr = 0.0;
	    if(xmc >= 0 && ymc >= 0 && zmc >= 0 && xmc <= reconSize-1 && ymc <= reconSize-1 && zmc <= reconSlices-1)
	    {
	      //double r = sqrt((xb0+xb1)*(xb0+xb1)+(-yb0-yb1)*(-yb0-yb1)+(-zb0-zb1)*(-zb0-zb1))*0.5;
	      double r = sqrt((xb0+xb1)*(xb0+xb1)+(-yb0-yb1)*(-yb0-yb1))*0.5;
	      for(int mat = 0; mat < num_mat; mat++)
		{
		  double ph = 1.0;
		  vnum = xmc + ymc*reconSize + (zmc*num_mat+mat)*reconSize*reconSize;
		  if(reconImageshort[vnum] < 0 || r >= rmax) ph = 0.0;
		  transp += reconImageshort[vnum] * Attenuation[mat*NE + ke] * rr * ph;
		}
	      //npf_float[vnum] += (float)(np_para * rr);
	      // total_length += rr;
	      //if(jj==256 || jj==255) fprintf(stderr,"%d %d %d %d %d %lf \n", xmc, ymc, zmc, vnum, reconImageshort[vnum], exp(-transp));
	      //fprintf(stderr,"%d %d %d %d %d %lf %lf\n", xmc, ymc, zmc, vnum, reconImageshort[vnum], transp, xb0);
	    }
	    //fprintf(stderr,"%lf %lf %lf %lf\n", aaa, xb1, yb1, zb1);
	    xb0=xb1;
	    yb0=yb1;
            zb0=zb1;
            ik++;
	  } //while
 	  //double intens = MAX_PIXVALUE - (MAX_AIR_PIXDIFF * exp(-transp/xfactor)*leng); // * atten;
	  double intens = pdf_pe[ke] * exp(-transp/xfactor);
	  if(ite == 0 && intens > 57000)
	  //if(ite == 0)
	    {
	      gf += log(np_para/MAX_AIR_PIXDIFF/leng) * transp/xfactor;
	      ff += transp/xfactor * transp/xfactor;
	      //xalpha += -gf/ff;
	    } 
 	  reprojection_float[jjjj] = (float)(intens);

	  
 	} //dx dy
 	  //if(jj == 256) fprintf(stderr,"%d %d %d %d \n", jjj, j, jj, ProjImageshort[jjjj]);
      } //Width
    } //Hight
  }
  
  //if(ite == 0) 
  //  {
  //    xalpha[0] = -gf/ff;
  //    xfactor = xfactor/xalpha[0];
  //    //fprintf(stderr,"xa = %lf xfac = %lf\n", xalpha[0], xfactor);
  //  }
  return;
}
