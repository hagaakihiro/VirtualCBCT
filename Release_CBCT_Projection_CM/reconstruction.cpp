#include <stdio.h>
#include <string.h>
#include <math.h>

#include "physParams.h"
#include "mallocMD.h"
#include "reconstruction.h"


void
reconstruction( int reconSize, int reconSlices, short* reconImageshort, unsigned short* ProjImageshort, float* npf_float, unsigned short* projVolume,
	      int projImgWidth, int projImgHeight, int usedProjNumber,
	      double* ProjAngle, double* cosTable, double* sinTable, double* XcenterShift, double* YcenterShift, double xfactor,
	      double DETECTOR_PITCH_MM, double DETECTOR_PITCH_MM_AT_ISO_CENT,
	      double DIST_BTWN_SRC_AND_ISOCENT_MM, double DETECTOR_PIXEL_NUM )
{
  int nxm = reconSize/2;
  int nym = reconSize/2;
  int nzm = reconSlices/2;

  double rmax = reconSize*0.5 * 0.9;
  //for(int jjj=0;jjj<usedProjNumber;jjj++)
  for(int jjj=0;jjj<1;jjj++)
  {
    //	int jjj=jjj0*70;
    //double aaa = ProjAngle[jjj];
    //aaa = aaa * 180 / PI;
    // Source location in cartesian coodinate
    double sx = -DIST_BTWN_SRC_AND_ISOCENT_MM*sinTable[jjj];
    double sy = DIST_BTWN_SRC_AND_ISOCENT_MM*cosTable[jjj];
    double sz = 0.0;
              
    //for(int j=135;j<projImgHeight-136;j++) // Hight
    for(int j=0;j<projImgHeight;j++) // Hight
    {
      double pz =  -(YcenterShift[jjj] + j) * DETECTOR_PITCH_MM_AT_ISO_CENT * (DETECTOR_PIXEL_NUM / projImgHeight);
      for(int jj=0;jj<projImgWidth;jj++)	// Width
      //for(int jj=511;jj<512;jj++)	// Width
      {
	int jjjj = jjj*projImgWidth*projImgHeight + j*projImgWidth + jj;
	// Detector location at IC
	double pr = DETECTOR_PITCH_MM_AT_ISO_CENT * (XcenterShift[jjj] + jj) * (DETECTOR_PIXEL_NUM / projImgWidth);
	double px = pr*cosTable[jjj];
	double py = pr*sinTable[jjj];
	//double leng = DIST_BTWN_SRC_AND_ISOCENT_MM*DIST_BTWN_SRC_AND_ISOCENT_MM/(DIST_BTWN_SRC_AND_ISOCENT_MM*DIST_BTWN_SRC_AND_ISOCENT_MM + (pr*pr+pz*pz));
	double leng = 1.0;
	//fprintf(stderr,"%lf %lf %lf %lf %lf %lf\n", ProjAngle[jjj]*180/PI, sx, sy, sz, px,py,pz  );
	double np_para = (MAX_PIXVALUE - projVolume[jjjj]) * log(MAX_AIR_PIXDIFF/(MAX_PIXVALUE - projVolume[jjjj]));
	double xx[4], yy[4], zz[4];
	int m = 0;
	yy[m] = double(nym);
	if(fabs((yy[m]-sy)/(py-sy)*(px-sx)+sx) <double(nxm)) 
	{
	  xx[m] = (yy[m]-sy)/(py-sy)*(px-sx)+sx;
	  zz[m] = (yy[m]-sy)/(py-sy)*(pz-sz)+sz;
	  m++;
	}
	yy[m] = -double(nym);
	if(fabs((yy[m]-sy)/(py-sy)*(px-sx)+sx) <double(nxm))
	{
	  xx[m] = (yy[m]-sy)/(py-sy)*(px-sx)+sx;
	  zz[m] = (yy[m]-sy)/(py-sy)*(pz-sz)+sz;
	  m++;
	}
	xx[m] = double(nxm);
	if(fabs((xx[m]-sx)/(px-sx)*(py-sy)+sy) <double(nym)) 
	{
	  yy[m] = (xx[m]-sx)/(px-sx)*(py-sy)+sy;
	  zz[m] = (xx[m]-sx)/(px-sx)*(pz-sz)+sz;
	  m++;
	}
	xx[m] = -double(nxm);
	if(fabs((xx[m]-sx)/(px-sx)*(py-sy)+sy) < double(nym))
	{
	  yy[m] = (xx[m]-sx)/(px-sx)*(py-sy)+sy;
	  zz[m] = (xx[m]-sx)/(px-sx)*(pz-sz)+sz;
	  m++;
	}
	double xi, yi, zi, xf, yf, zf; // Initial and Final Voxels, which X-ray passes through
	if((xx[0]-sx)*(xx[0]-sx)+(yy[0]-sy)*(yy[0]-sy)+(zz[0]-sz)*(zz[0]-sz) < (xx[1]-sx)*(xx[1]-sx)+(yy[1]-sy)*(yy[1]-sy)+(zz[1]-sz)*(zz[1]-sz))
	{
	  xi=xx[0]; yi=yy[0]; zi=zz[0]; xf=xx[1]; yf=yy[1]; zf=zz[1];
	}
	else
	{
	  xi=xx[1]; yi=yy[1]; zi=zz[1]; xf=xx[0]; yf=yy[0]; zf=zz[0];
	}
	//fprintf(stderr,"%lf %lf %lf %lf %lf %lf\n", aaa, xi, yi, zi, xf,yf  );
	double xb0 = xi, yb0 = yi, zb0 = zi, xb1, yb1, zb1;
	double xii, yii, zii, xiii, yiii, ziii, xiiii, yiiii, ziiii, yj, xj, rr, transp, total_length;
	int xma, yma, zma, xmb, ymb, zmb, xmc, ymc, zmc, icheck, vnum;
	if(fabs(xi-xf) > fabs(yi-yf)) // dx is larger than dy   *****************************************************
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
	  while( xii <= double(nxm) && xii >= -double(nxm) &&  yii <= double(nym) && yii >= -double(nym))
	  {
	    xii = xin + pha * ik;
	    ymb=yma;
	    yj = yii;
	    yii = (xii-sx)/(px-sx)*(py-sy)+sy;
	    yma=(int)yii;
	    zmb=zma;
	    zii = (xii-sx)/(px-sx)*(pz-sz)+sz;
	    zma=(int)zii;
	    icheck = 0;
	    chy = 0;
	    chz = 0;
	    //fprintf(stderr,"%lf %lf %lf %lf\n", aaa, xii, yii, zii);
	    if((ik != 0 && yma != ymb) || (ik != 0 && (yii/fabs(yii) != yj/fabs(yj))))
	    {
	      if((yii < 0 && (yf-yi) < 0) || (yii > 0 && (yf-yi) > 0)) yiii = (double)yma;
	      else yiii = (double)yma - (yf-yi)/fabs(yf-yi);
	      xiii = (yiii-sy)/(py-sy)*(px-sx)+sx;
	      ziii = (yiii-sy)/(py-sy)*(pz-sz)+sz;
	      chy = (xii-xiii)*(xii-xiii)+(yii-yiii)*(yii-yiii)+(zii-ziii)*(zii-ziii);
	      icheck++;
	      //fprintf(stderr,"%lf %lf %lf %lf\n", aaa, xiii, yiii, ziii);
	      //if(aaa < 40 && aaa > 70) fprintf(stderr,"%lf %lf %lf %lf\n", aaa, xiii, yiii, ziii);
            }
            if(ik != 0 && zma != zmb)
            {
                ziiii = (double)zma;
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
                xmc = (int)((xb0+xb1)*0.5+nxm);
                ymc = (int)((-yb0-yb1)*0.5+nym);
                zmc = (int)((-zb0-zb1)*0.5+nzm);
		//fprintf(stderr,"%lf %lf %lf %d %d %d \n", (xb0+xb1)*0.5,-(-yb0-yb1)*0.5,-(-zb0-zb1)*0.5,xmc,ymc,zmc);/////////////////////
                if(xmc >= 0 && ymc >= 0 && zmc >= 0 && xmc <= reconSize-1 && ymc <= reconSize-1 && zmc <= reconSlices-1)
                {
                    vnum = xmc + ymc*reconSize + zmc*reconSize*reconSize;
 		    double r = sqrt((xb0+xb1)*(xb0+xb1)+(-yb0-yb1)*(-yb0-yb1)+(-zb0-zb1)*(-zb0-zb1))*0.5;
                    if(reconImageshort[vnum] < 0 || r >= rmax) reconImageshort[vnum]=0;
                    transp += reconImageshort[vnum] * rr / xfactor;
		    npf_float[vnum] += (float)(np_para * rr);
 		    total_length += rr;
		    //fprintf(stderr,"%d %d %d %d %lf %lf \n", xmc, ymc, zmc, vnum, transp, total_length);
                }
		xb0=xb1;
		yb0=yb1;
		zb0=zb1;
		iks++;
		xb1=xiii;
		yb1=yiii;
		zb1=ziii;
		rr = sqrt((xb0-xb1)*(xb0-xb1)+(yb0-yb1)*(yb0-yb1)+(zb0-zb1)*(zb0-zb1));
		xmc = (int)((xb0+xb1)*0.5+nxm);
		ymc = (int)((-yb0-yb1)*0.5+nym);
		zmc = (int)((-zb0-zb1)*0.5+nzm);
		//fprintf(stderr,"%lf %lf %lf %d %d %d \n", (xb0+xb1)*0.5,-(-yb0-yb1)*0.5,-(-zb0-zb1)*0.5,xmc,ymc,zmc);/////////////////////
		if(xmc >= 0 && ymc >= 0 && zmc >= 0 && xmc <= reconSize-1 && ymc <= reconSize-1 && zmc <= reconSlices-1)
		{
		  vnum = xmc + ymc*reconSize + zmc*reconSize*reconSize;
 		  double r = sqrt((xb0+xb1)*(xb0+xb1)+(-yb0-yb1)*(-yb0-yb1)+(-zb0-zb1)*(-zb0-zb1))*0.5;
		  if(reconImageshort[vnum] < 0 || r >= rmax) reconImageshort[vnum]=0;
		  transp += reconImageshort[vnum] * rr / xfactor;
		  npf_float[vnum] += (float)(np_para * rr);
 		  total_length += rr;
		  //fprintf(stderr,"%d %d %d %d %lf %lf \n", xmc, ymc, zmc, vnum, transp, total_length);
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
	      xmc = (int)((xb0+xb1)*0.5+nxm);
	      ymc = (int)((-yb0-yb1)*0.5+nym);
	      zmc = (int)((-zb0-zb1)*0.5+nzm);
	      //fprintf(stderr,"%lf %lf %lf %d %d %d \n", (xb0+xb1)*0.5,-(-yb0-yb1)*0.5,-(-zb0-zb1)*0.5,xmc,ymc,zmc);/////////////////////
	      if(xmc >= 0 && ymc >= 0 && zmc >= 0 && xmc <= reconSize-1 && ymc <= reconSize-1 && zmc <= reconSlices-1)
	      {
		vnum = xmc + ymc*reconSize + zmc*reconSize*reconSize;
		double r = sqrt((xb0+xb1)*(xb0+xb1)+(-yb0-yb1)*(-yb0-yb1)+(-zb0-zb1)*(-zb0-zb1))*0.5;
		if(reconImageshort[vnum] < 0 || r >= rmax) reconImageshort[vnum]=0;
		transp += reconImageshort[vnum] * rr / xfactor;
		npf_float[vnum] += (float)(np_para * rr);
		total_length += rr;
		//fprintf(stderr,"%d %d %d %d %lf %lf \n", xmc, ymc, zmc, vnum, transp, total_length);
	      }
	      xb0=xb1;
	      yb0=yb1;
	      zb0=zb1;
	      iks++;
	      xb1=xiiii;
              yb1=yiiii;
              zb1=ziiii;
              rr = sqrt((xb0-xb1)*(xb0-xb1)+(yb0-yb1)*(yb0-yb1)+(zb0-zb1)*(zb0-zb1));
              xmc = (int)((xb0+xb1)*0.5+nxm);
              ymc = (int)((-yb0-yb1)*0.5+nym);
              zmc = (int)((-zb0-zb1)*0.5+nzm);
	      //fprintf(stderr,"%lf %lf %lf %d %d %d \n", (xb0+xb1)*0.5,-(-yb0-yb1)*0.5,-(-zb0-zb1)*0.5,xmc,ymc,zmc);/////////////////////
              if(xmc >= 0 && ymc >= 0 && zmc >= 0 && xmc <= reconSize-1 && ymc <= reconSize-1 && zmc <= reconSlices-1)
              {
                      vnum = xmc + ymc*reconSize + zmc*reconSize*reconSize;
    		      double r = sqrt((xb0+xb1)*(xb0+xb1)+(-yb0-yb1)*(-yb0-yb1)+(-zb0-zb1)*(-zb0-zb1))*0.5;
                      if(reconImageshort[vnum] < 0 || r >= rmax) reconImageshort[vnum]=0;
                      transp += reconImageshort[vnum] * rr / xfactor;
    		      npf_float[vnum] += (float)(np_para * rr);
    		      total_length += rr;
		      //fprintf(stderr,"%d %d %d %d %lf %lf \n", xmc, ymc, zmc, vnum, transp, total_length);
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
                  xmc = (int)((xb0+xb1)*0.5+nxm);
                  ymc = (int)((-yb0-yb1)*0.5+nym);
                  zmc = (int)((-zb0-zb1)*0.5+nzm);
		  //fprintf(stderr,"%lf %lf %lf %d %d %d \n", (xb0+xb1)*0.5,-(-yb0-yb1)*0.5,-(-zb0-zb1)*0.5,xmc,ymc,zmc);/////////////////////
                  if(xmc >= 0 && ymc >= 0 && zmc >= 0 && xmc <= reconSize-1 && ymc <= reconSize-1 && zmc <= reconSlices-1)
                  {
                      vnum = xmc + ymc*reconSize + zmc*reconSize*reconSize;
    		      double r = sqrt((xb0+xb1)*(xb0+xb1)+(-yb0-yb1)*(-yb0-yb1)+(-zb0-zb1)*(-zb0-zb1))*0.5;
                      if(reconImageshort[vnum] < 0 || r >= rmax) reconImageshort[vnum]=0;
                      transp += reconImageshort[vnum] * rr / xfactor;
		      npf_float[vnum] += (float)(np_para * rr);
    		      total_length += rr;
		      //fprintf(stderr,"%d %d %d %d %lf %lf \n", xmc, ymc, zmc, vnum, transp, total_length);
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
                  xmc = (int)((xb0+xb1)*0.5+nxm);
                  ymc = (int)((-yb0-yb1)*0.5+nym);
                  zmc = (int)((-zb0-zb1)*0.5+nzm);
		  //fprintf(stderr,"%lf %lf %lf %d %d %d \n", (xb0+xb1)*0.5,-(-yb0-yb1)*0.5,-(-zb0-zb1)*0.5,xmc,ymc,zmc);/////////////////////
                  if(xmc >= 0 && ymc >= 0 && zmc >= 0 && xmc <= reconSize-1 && ymc <= reconSize-1 && zmc <= reconSlices-1)
                  {
                      vnum = xmc + ymc*reconSize + zmc*reconSize*reconSize;
    		      double r = sqrt((xb0+xb1)*0.5*(xb0+xb1)*0.5+(-yb0-yb1)*0.5*(-yb0-yb1)*0.5+(-zb0-zb1)*0.5*(-zb0-zb1)*0.5);
                      if(reconImageshort[vnum] < 0 || r >= rmax) reconImageshort[vnum]=0;
                      transp += reconImageshort[vnum] * rr / xfactor;
		      npf_float[vnum] += (float)(np_para * rr);
    		      total_length += rr;
		      //fprintf(stderr,"%d %d %d %d %lf %lf \n", xmc, ymc, zmc, vnum, transp, total_length);
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
              xmc = (int)((xb0+xb1)*0.5+nxm);
              ymc = (int)((-yb0-yb1)*0.5+nym);
              zmc = (int)((-zb0-zb1)*0.5+nzm);
	      //fprintf(stderr,"%lf %lf %lf %d %d %d \n", (xb0+xb1)*0.5,-(-yb0-yb1)*0.5,-(-zb0-zb1)*0.5,xmc,ymc,zmc);/////////////////////
              if(xmc >= 0 && ymc >= 0 && zmc >= 0 && xmc <= reconSize-1 && ymc <= reconSize-1 && zmc <= reconSlices-1)
              {
                  vnum = xmc + ymc*reconSize + zmc*reconSize*reconSize;
    					  double r = sqrt((xb0+xb1)*0.5*(xb0+xb1)*0.5+(-yb0-yb1)*0.5*(-yb0-yb1)*0.5+(-zb0-zb1)*0.5*(-zb0-zb1)*0.5);
                  if(reconImageshort[vnum] < 0 || r >= rmax) reconImageshort[vnum] = 0;
                  transp += reconImageshort[vnum] * rr / xfactor;
		  npf_float[vnum] += (float)(np_para * rr);
    		  total_length += rr;
		  //fprintf(stderr,"%d %d %d %d %lf %lf \n", xmc, ymc, zmc, vnum, transp, total_length);
              }
              //fprintf(stderr,"%lf %d %d %lf %lf %lf %lf *\n", aaa, iks, vnum, xii, yii, zii, total_length);
              //if(aaa > 75 && aaa < 76) fprintf(stderr,"%lf %d %d %d %d %d %d %lf\n", aaa, iks, vnum, xmc, ymc, zmc, reconImageshort[vnum], transp);
              xb0=xb1;
              yb0=yb1;
              zb0=zb1;
              ik++;
              //fprintf(stderr,"%lf %d %d %lf %lf %lf %lf\n", aaa, iks, nzm, xii, yii, zii, rr);
	      //fprintf(stderr,"%lf %d %d %d %d %d %d %lf\n", aaa, iks, vnum, xmc, ymc, zmc, reconImageshort[vnum], transp);
          } // while
    	  double intens = MAX_PIXVALUE - (MAX_AIR_PIXDIFF * exp(-transp)); // * atten;
	  ProjImageshort[jjjj] = (unsigned short)(intens*leng);
	  //double length = sqrt((xi-xf)*(xi-xf)+(yi-yf)*(yi-yf)+(zi-zf)*(zi-zf));
	  //fprintf(stderr,"%d %d %lf %lf %lf %lf %lf %lf %lf %e %e \n", j, jj, transp, xi,yi,zi,xf,yf,zf, leng, total_length);
	  //if(aaa > 75 && aaa < 76) ProjImageshort[j][jj] = (unsigned short)(intens);
	  //if(aaa > 75 && aaa < 76) fprintf(stderr,"%lf %d %d %d %lf\n", aaa, j, jj, ProjImageshort[j][jj], atten);
        } // dx dy 
	else if(fabs(xi-xf) < fabs(yi-yf)) // dy is larger than dx *********************************************************************
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
	  while( xii <= double(nxm) && xii >= -double(nxm) &&  yii <= double(nym) && yii >= -double(nym))
	  {
	    yii = yin + pha * ik;
	    xmb=xma;
	    xj = xii;
	    xii = (yii-sy)/(py-sy)*(px-sx)+sx;
	    xma=(int)xii;
	    zmb=zma;
	    zii = (yii-sy)/(py-sy)*(pz-sz)+sz;
	    zma=(int)zii;
	    icheck = 0;
	    chy = 0;
	    chz = 0;
	    if((ik != 0 && xma != xmb) || (ik != 0 && (xii/fabs(xii) != xj/fabs(xj))))
	    {
	      if((xii < 0 && (xf-xi) < 0) || (xii > 0 && (xf-xi) > 0)) xiii = (double)xma;
	      else xiii = (double)xma - (xf-xi)/fabs(xf-xi);
	      yiii = (xiii-sx)/(px-sx)*(py-sy)+sy;
	      ziii = (xiii-sx)/(px-sx)*(pz-sz)+sz;
	      chy = (xii-xiii)*(xii-xiii)+(yii-yiii)*(yii-yiii)+(zii-ziii)*(zii-ziii);
	      icheck++;
	      //if(aaa > 160 && aaa < 161) fprintf(stderr,"%lf %lf %lf %lf\n", aaa, xiii, yiii, ziii);
	    }
	    if(ik != 0 && zma != zmb)
	    {
	      ziiii = (double)zma;
	      xiiii = (ziiii-sz)/(pz-sz)*(px-sx)+sx;
	      yiiii = (ziiii-sz)/(pz-sz)*(py-sy)+sy;
	      chz = (xii-xiiii)*(xii-xiiii)+(yii-yiiii)*(yii-yiiii)+(zii-ziiii)*(zii-ziiii);
	      icheck++;
	      //if(aaa > 160 && aaa < 161) fprintf(stderr,"%lf %lf %lf %lf\n", aaa, xiiii, yiiii, ziiii);
	    }
 			  				  
	    if(icheck == 2 && chz > chy)
	    {
	      iks++;
	      xb1=xiiii;
	      yb1=yiiii;
	      zb1=ziiii;
	      rr = sqrt((xb0-xb1)*(xb0-xb1)+(yb0-yb1)*(yb0-yb1)+(zb0-zb1)*(zb0-zb1));
	      xmc = (int)((xb0+xb1)*0.5+nxm);
	      ymc = (int)((-yb0-yb1)*0.5+nym);
	      zmc = (int)((-zb0-zb1)*0.5+nzm);
	      if(xmc >= 0 && ymc >= 0 && zmc >= 0 && xmc <= reconSize-1 && ymc <= reconSize-1 && zmc <= reconSlices-1)
	      {
		vnum = xmc + ymc*reconSize + zmc*reconSize*reconSize;
		double r = sqrt((xb0+xb1)*0.5*(xb0+xb1)*0.5+(-yb0-yb1)*0.5*(-yb0-yb1)*0.5+(-zb0-zb1)*0.5*(-zb0-zb1)*0.5);
		if(reconImageshort[vnum] < 0 || r >= rmax) reconImageshort[vnum]=0;
		transp += reconImageshort[vnum] * rr / xfactor;
		npf_float[vnum] += (float)(np_para * rr);
		total_length += rr;
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
	      xmc = (int)((xb0+xb1)*0.5+nxm);
	      ymc = (int)((-yb0-yb1)*0.5+nym);
	      zmc = (int)((-zb0-zb1)*0.5+nzm);
	      if(xmc >= 0 && ymc >= 0 && zmc >= 0 && xmc <= reconSize-1 && ymc <= reconSize-1 && zmc <= reconSlices-1)
	      {
		vnum = xmc + ymc*reconSize + zmc*reconSize*reconSize;
		double r = sqrt((xb0+xb1)*0.5*(xb0+xb1)*0.5+(-yb0-yb1)*0.5*(-yb0-yb1)*0.5+(-zb0-zb1)*0.5*(-zb0-zb1)*0.5);
		if(reconImageshort[vnum] < 0 || r >= rmax) reconImageshort[vnum]=0;
		transp += reconImageshort[vnum] * rr / xfactor;
		npf_float[vnum] += (float)(np_para * rr);
		total_length += rr;
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
  	      xmc = (int)((xb0+xb1)*0.5+nxm);
  	      ymc = (int)((-yb0-yb1)*0.5+nym);
	      zmc = (int)((-zb0-zb1)*0.5+nzm);
	      if(xmc >= 0 && ymc >= 0 && zmc >= 0 && xmc <= reconSize-1 && ymc <= reconSize-1 && zmc <= reconSlices-1)
	      {
		vnum = xmc + ymc*reconSize + zmc*reconSize*reconSize;
  		double r = sqrt((xb0+xb1)*0.5*(xb0+xb1)*0.5+(-yb0-yb1)*0.5*(-yb0-yb1)*0.5+(-zb0-zb1)*0.5*(-zb0-zb1)*0.5);
		if(reconImageshort[vnum] < 0 || r >= rmax) reconImageshort[vnum]=0;
		transp += reconImageshort[vnum] * rr / xfactor;
		npf_float[vnum] += (float)(np_para * rr);
		total_length += rr;
	      }
	      //if(aaa > 160 && aaa < 161) fprintf(stderr,"%lf %lf %lf %lf\n", aaa, xb1, yb1, zb1);
	      xb0=xb1;
	      yb0=yb1;
	      zb0=zb1;
	      iks++;
	      xb1=xiiii;
	      yb1=yiiii;
	      zb1=ziiii;
	      rr = sqrt((xb0-xb1)*(xb0-xb1)+(yb0-yb1)*(yb0-yb1)+(zb0-zb1)*(zb0-zb1));
	      xmc = (int)((xb0+xb1)*0.5+nxm);
	      ymc = (int)((-yb0-yb1)*0.5+nym);
	      zmc = (int)((-zb0-zb1)*0.5+nzm);
	      if(xmc >= 0 && ymc >= 0 && zmc >= 0 && xmc <= reconSize-1 && ymc <= reconSize-1 && zmc <= reconSlices-1)
	      {
		vnum = xmc + ymc*reconSize + zmc*reconSize*reconSize;
		double r = sqrt((xb0+xb1)*0.5*(xb0+xb1)*0.5+(-yb0-yb1)*0.5*(-yb0-yb1)*0.5+(-zb0-zb1)*0.5*(-zb0-zb1)*0.5);
		if(reconImageshort[vnum] < 0 || r >= rmax) reconImageshort[vnum]=0;
		transp += reconImageshort[vnum] * rr / xfactor;
	        npf_float[vnum] += (float)(np_para * rr);
		total_length += rr;
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
	      xmc = (int)((xb0+xb1)*0.5+nxm);
	      ymc = (int)((-yb0-yb1)*0.5+nym);
	      zmc = (int)((-zb0-zb1)*0.5+nzm);
	      if(xmc >= 0 && ymc >= 0 && zmc >= 0 && xmc <= reconSize-1 && ymc <= reconSize-1 && zmc <= reconSlices-1)
	      {
		vnum = xmc + ymc*reconSize + zmc*reconSize*reconSize;
		double r = sqrt((xb0+xb1)*0.5*(xb0+xb1)*0.5+(-yb0-yb1)*0.5*(-yb0-yb1)*0.5+(-zb0-zb1)*0.5*(-zb0-zb1)*0.5);
		if(reconImageshort[vnum] < 0 || r >= rmax) reconImageshort[vnum]=0;
		transp += reconImageshort[vnum] * rr / xfactor;
		npf_float[vnum] += (float)(np_para * rr);
		total_length += rr;
	      }
	      //fprintf(stderr,"%lf %lf %lf %lf\n", aaa, xb1, yb1, zb1);
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
	      xmc = (int)((xb0+xb1)*0.5+nxm);
	      ymc = (int)((-yb0-yb1)*0.5+nym);
	      zmc = (int)((-zb0-zb1)*0.5+nzm);
	      if(xmc >= 0 && ymc >= 0 && zmc >= 0 && xmc <= reconSize-1 && ymc <= reconSize-1 && zmc <= reconSlices-1)
	      {
		vnum = xmc + ymc*reconSize + zmc*reconSize*reconSize;
		double r = sqrt((xb0+xb1)*0.5*(xb0+xb1)*0.5+(-yb0-yb1)*0.5*(-yb0-yb1)*0.5+(-zb0-zb1)*0.5*(-zb0-zb1)*0.5);
		if(reconImageshort[vnum] < 0 || r >= rmax) reconImageshort[vnum]=0;
		transp += reconImageshort[vnum] * rr / xfactor;
		npf_float[vnum] += (float)(np_para * rr);
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
	    xmc = (int)((xb0+xb1)*0.5+nxm);
	    ymc = (int)((-yb0-yb1)*0.5+nym);
	    zmc = (int)((-zb0-zb1)*0.5+nzm);
	    //if(ik == 0) rr = 0.0;
	    if(xmc >= 0 && ymc >= 0 && zmc >= 0 && xmc <= reconSize-1 && ymc <= reconSize-1 && zmc <= reconSlices-1)
	    {
	      vnum = xmc + ymc*reconSize + zmc*reconSize*reconSize;
	      double r = sqrt((xb0+xb1)*0.5*(xb0+xb1)*0.5+(-yb0-yb1)*0.5*(-yb0-yb1)*0.5+(-zb0-zb1)*0.5*(-zb0-zb1)*0.5);
	      if(reconImageshort[vnum] < 0 || r >= rmax) reconImageshort[vnum]=0;
	      transp += reconImageshort[vnum] * rr / xfactor;
	      npf_float[vnum] += (float)(np_para * rr);
	      total_length += rr;
	    }
	    //fprintf(stderr,"%lf %lf %lf %lf\n", aaa, xb1, yb1, zb1);
	    xb0=xb1;
	    yb0=yb1;
            zb0=zb1;
            ik++;
	  } //while
 	  double intens = MAX_PIXVALUE - (MAX_AIR_PIXDIFF * exp(-transp)); // * atten;
 	  ProjImageshort[jjjj] = (unsigned short)(intens*leng);
 	  //double length = sqrt((xi-xf)*(xi-xf)+(yi-yf)*(yi-yf)+(zi-zf)*(zi-zf));
 	  //fprintf(stderr,"%d %d %lf %lf %lf %lf %lf %lf %lf %e %e\n", j, jj, transp, xi,yi,zi,xf,yf,zf, leng, total_length);
 	  //fprintf(stderr,"%lf %d %d %d %lf\n", aaa, j, jj, ProjImageshort[j][jj], transp);
 	  //if(aaa > 5 && aaa < 6) ProjImageshort[j][jj] = (unsigned short)(intens);
 	  //if(aaa > 5 && aaa < 6) fprintf(stderr,"%lf %d %d %d %lf\n", aaa, j, jj, ProjImageshort[j][jj], atten);
 	} //dx dy
      }
    }
  }

}
