#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>

#include "physParams.h"
#include "mallocMD.h"
#include "virtual_projection.h"
#include "reprojection.h"

#define LINEBUF_SIZE	512

  // Uniform distribution
  double Uniform( void ){
    return ((double)rand()+1.0)/((double)RAND_MAX+2.0);
  }

  double rand_gauss( double mu, double sigma ){
    //double M_PI = 3.14159265;
    double z = sqrt( -2.0*log(Uniform()) ) * sin( 2.0*PI*Uniform() );
    return mu + sigma*z;
  }

  int gaussiannoise(float* ProjImagefloat, int NN)
  {
    double mu=0;
    for(int i = 0; i < NN; i++)
      {
	double nstar = MAX_PIXVALUE - ProjImagefloat[i];
	double sigma=sqrt(nstar);
	double psn = rand_gauss(0,sigma);
	nstar = nstar + psn;
	ProjImagefloat[i] = MAX_PIXVALUE - nstar;
	if(MAX_PIXVALUE - nstar < 0)  ProjImagefloat[i] = 0;
	if(MAX_PIXVALUE - nstar > MAX_PIXVALUE)  ProjImagefloat[i] = (float)MAX_PIXVALUE;
      }
    return 0;
  }


void
prior_weight_production( int pri_size, int reconSize, int reconSlices, int num_material, 
			 float* ProjImagefloat, float* npf_float, float* npf_float_re, float* projVolume,
			 int projImgWidth, int projImgHeight, int usedProjNumber, int startProjNumber, int endProjNumber,
			 double* anglesRad, double* cosTable, double* sinTable, double* XcenterShift, double* YcenterShift, 
			 short* req_reconImageshort, short* Pri_reconImageshort, short* Pri_reconImageshort_0, double* Attenuation, double* pdf_pe, int NE, double xfactor)
{
  ///////////////  Prior material image production ///////////////////
  FILE*	fp;
  FILE* fp3;
  char filename[128];
  int mnum=0, jmax;
  float dens[76], hdens[76], cdens[76], ndens[76], odens[76], pdens[76], cadens[76];
  if((fp = fopen("material_density_mod.txt","r")) == NULL ) {
                printf("Can not open\n");
                return;
        }
  while(fscanf(fp, "%d %f %f %f %f %f %f %f\n", &jmax, &dens[mnum], &hdens[mnum], &cdens[mnum], 
	       &ndens[mnum], &odens[mnum], &pdens[mnum], &cadens[mnum]) != EOF)
    {
      printf("%d %f %f \n", mnum, dens[mnum], hdens[mnum]);
      mnum++;
    }
  fclose(fp);
  
  int mat0[3] = {1,2,3};
  int mat11[6] = {4,9,14,26,39,49};
  int mat12[6] = {5,10,18,30,42,50};
  int mat13[6] = {6,11,22,35,45,51};
  int mat21[3] = {52,59,73};
  int mat22[3] = {53,63,75};
  int mat23[3] = {54,68,76};

  for( int i = 0; i < 3; i++)
    for(int j = 0; j < 6; j++)
      for(int k = 0; k < 3; k++)
	{
	  sprintf(filename,"ph3D300SL16bit.dat");///////////////////////////////////////////////////////////////////////////
	  int ict1 = 1310;
	  int ict2 = 1965;
	  if( (fp=fopen(filename,"rb")) == NULL )
	    {
	      fprintf(stderr,"ERROR: No output CT file \n");
	      return;
	    }
	  int a = fread( Pri_reconImageshort_0, 1, pri_size * sizeof(short), fp );
	  fclose(fp);
	  
	  int material1 = mat0[i]-1;
	  int material2 = mat11[j]-1;
	  if(2 == mat0[i]) material2 = mat12[j]-1;
	  if(3 == mat0[i]) material2 = mat13[j]-1;
	  int material3 = mat21[k]-1;
	  if(2 == mat0[i]) material3 = mat22[k]-1;
	  if(3 == mat0[i]) material3 = mat23[k]-1;
	  sprintf(filename,"ph3D300SL16bit_%d_%d_%d.dat", material1, material2, material3);///////////////////////////////////////////////////////////////////////////
	  fprintf(stderr,"%s %f %f %f \n", filename, dens[material1],dens[material2],dens[material3]);
	  fp=fopen(filename,"wb");
	  for(int ii = 0; ii < pri_size; ii++) 
	    {
	      int i0 = ii/reconSize/reconSize;
	      int i1 = (ii-i0*reconSize*reconSize)/reconSize;
	      int i2 = ii-i1*reconSize-i0*reconSize*reconSize;
	      if(Pri_reconImageshort_0[ii] < 0) Pri_reconImageshort_0[ii] = 0;
	      if(Pri_reconImageshort_0[ii] == ict1) Pri_reconImageshort_0[ii] = (short)(dens[material2]*1000);
	      else if(Pri_reconImageshort_0[ii] == ict2) Pri_reconImageshort_0[ii] = (short)(dens[material3]*1000);
	      else if(Pri_reconImageshort_0[ii] == 0 && (i1 <= 210 && i1 >= 80) && (i2 <= 195 && i2 >= 75) && (i0 <= 194 && i0 >= 108))  
		Pri_reconImageshort_0[ii] = (short)(dens[material1]*1000);
	      else Pri_reconImageshort_0[ii] = 0;
	      
	    }
	  fwrite( Pri_reconImageshort_0,
		  getByteSizeOfUnit(UNIT_SINT16),pri_size,fp);
	  fclose(fp);
	  
	  float rho1 = dens[material1]*1000;
	  float rho2 = dens[material2]*1000;
	  float rho3 = dens[material3]*1000;
	  float wei1[6], wei2[6], wei3[6];
	  wei1[0] = hdens[material1]; wei1[1] = cdens[material1]; wei1[2] = ndens[material1]; wei1[3] = odens[material1]; wei1[4] = pdens[material1], wei1[5] = cadens[material1];  
	  wei2[0] = hdens[material2]; wei2[1] = cdens[material2]; wei2[2] = ndens[material2]; wei2[3] = odens[material2]; wei2[4] = pdens[material2], wei2[5] = cadens[material2];  
	  wei3[0] = hdens[material3]; wei3[1] = cdens[material3]; wei3[2] = ndens[material3]; wei3[3] = odens[material3]; wei3[4] = pdens[material3], wei3[5] = cadens[material3];  
	  
	  //exit(0);
	  //
	  
	  
	  ///////// Material density (True) //////////
	  for(int iphase = 0; iphase < num_material; iphase++)
	    {
	      for(int i = 0; i < pri_size; i++) 
		{
		  if(Pri_reconImageshort_0[i] < 0) Pri_reconImageshort_0[i] = 0;
		  Pri_reconImageshort[i+pri_size*iphase] = 0;
		  if(Pri_reconImageshort_0[i] == (short)(dens[material1]*1000)) Pri_reconImageshort[i+pri_size*iphase] = (short)(wei1[iphase]*1000);
		  if(Pri_reconImageshort_0[i] == (short)(dens[material2]*1000)) Pri_reconImageshort[i+pri_size*iphase] = (short)(wei2[iphase]*1000);
		  if(Pri_reconImageshort_0[i] == (short)(dens[material3]*1000)) Pri_reconImageshort[i+pri_size*iphase] = (short)(wei3[iphase]*1000);
		}  
	    }
	  /*
	    sprintf(filename,"Weight_input.raw");
	    fp3=fopen(filename,"wb");	      
	    fwrite( Pri_reconImageshort,
	    getByteSizeOfUnit(UNIT_SINT16),pri_size*num_material,fp3);
	    fclose(fp3);
	  */
		      
	  int ire0=0;
	  for(int i1 = 0; i1 < num_material; i1++)
	    for(int i2 = 149; i2 < 150; i2++)
	      for(int i3 = 0; i3 < reconSize; i3++)
		for(int i4 = 0; i4 < reconSize; i4++)
		  {
		    req_reconImageshort[ire0] = Pri_reconImageshort[i1*reconSlices*reconSize*reconSize
								    +i2*reconSize*reconSize+i3*reconSize+i4]; 
		    ire0++;
		  }
	  // Pri_Weight: True Material Weight distribution
	  sprintf(filename,"Pri_Weight_input_%d_%d_%d.raw", material1, material2, material3);
	  if((fp3 = fopen(filename,"wb")) == NULL ) {
	    printf("Can not open\n");
	    return;
	  }
	  fwrite( req_reconImageshort,
		  getByteSizeOfUnit(UNIT_SINT16),
		  reconSize*reconSize*num_material,fp3);
	  fclose(fp3);
	  
	  //////////////////////// Sinogram production //////////////////
	  double xfac0, xalpha[1];
	  xalpha[0]= 1;
	  /*
	  reprojection( 1, reconSize, reconSlices, Pri_reconImageshort, 
			ProjImagefloat, npf_float, npf_float_re, projVolume,
			1, 1,
			projImgWidth, projImgHeight, startProjNumber, endProjNumber,
			anglesRad, cosTable, sinTable, XcenterShift, YcenterShift, 
			xfactor, xalpha,
			DETECTOR_PITCH_CM_kV, DETECTOR_PITCH_CM_AT_ISO_CENT_kV,
			DIST_BTWN_SRC_AND_ISOCENT_CM_kV, DETECTOR_PIXEL_NUM_kV,
			0, 0, 0,
			Attenuation, pdf_pe, num_material, NE);
	  */
	  // Gaussian Noise Input
	  gaussiannoise(ProjImagefloat,projImgWidth*projImgHeight*usedProjNumber);
	  fprintf(stderr,"ReprojImagefloat_data write \n");
	  sprintf(filename,"ReprojImagefloat_%d_%d_%d.raw", material1, material2, material3);
	  fp3=fopen(filename,"wb");
	  fwrite( ProjImagefloat,
		  getByteSizeOfUnit(UNIT_FLOAT32),
		  projImgWidth*projImgHeight*usedProjNumber,fp3); // Data
	  fclose(fp3);
	  //for(int i = 0; i < projImgWidth*projImgHeight; i++)
	  //  projData->projVolume[i] = ProjImagefloat[i];
 	  //////////////////////// Virtual projection image production -- END -- ////////////////////////////	  
	}
  return;
}

void
attenuation_coefficient( char* spectrum_data_name, int num_material, float* Attenuation, float* pdf_pe, int *NE, int nthin, char* mat_type)
{
  ///////////////  Attenuation data (using materials and x-ray energies) ///////////////////
  //// Photon Energy Spectrum
  FILE*	fp;
  FILE* fp3;
  float data0, data1, data2, sum=0;
  char  row[1024];
  int   NNE=0;
  double photonenergy[1024];
  if( (fp=fopen(spectrum_data_name,"r")) == NULL )
    {
      fprintf(stderr,"ERROR: no photon energy file \n");
      return;
    }
  if(fgets(row, sizeof(row), fp) == NULL)
    {
      fprintf(stderr,"ERROR: no photon energy file \n");
    }
  while( fgets(row, sizeof(row), fp) != NULL)
    {
      //sscanf(row, "%f %f %f", &data0, &data1, &data2);//txt
      sscanf(row, "%f,%f,%f", &data0, &data1, &data2);//csv
      photonenergy[NNE] = (data1+data0)*0.5;
      pdf_pe[NNE] = data2;
      sum += data2;
      NNE++;
    }
  fclose(fp);
  //
  double psum = 0.0, esum = 0.0, sum_thinout = 0.0;
  int NNE_thinout;
  float pdf_thinout[1024], pe_thinout[1024];
  for (int i = 0; i < NNE; i++)
    {
      psum += pdf_pe[i];
      esum += photonenergy[i];
      if(i-nthin+1 == (i-nthin+1)/nthin*nthin) 
	{
	  pdf_thinout[(i-nthin+1)/nthin] = psum;
	  pe_thinout[(i-nthin+1)/nthin] = esum/nthin;
	  sum_thinout += psum; 
	  psum = 0.0;
	  esum = 0.0;
	  //printf( "%d %f %f \n", (i-nthin+1)/nthin, pe_thinout[(i-nthin+1)/nthin], pdf_thinout[(i-nthin+1)/nthin]);
	  NNE_thinout = (i-nthin+1)/nthin+1;
	}
    }
  //

  

  
  NNE = NNE_thinout;
  sum = sum_thinout;
  for (int i = 0; i < NNE; i++) 
    {
      pdf_pe[i] = pdf_thinout[i];
      photonenergy[i] = pe_thinout[i];
    }
  float norm = 0.0;
  float mean_energy;
  for (int i = 0; i < NNE; i++)
    {
      pdf_pe[i] = pdf_pe[i]/sum;
      norm += pdf_pe[i];
      mean_energy += photonenergy[i]*pdf_pe[i];
      printf( "%d %f %f %f \n", i, photonenergy[i], pdf_pe[i], norm);
    }
  printf("mean x-ray energy = %f \n",mean_energy);
 
  ////////////////////////////
  //  Unit [cm^2/g]; Mass attenuation coefficient
  //  H(1,1), C(6,12), N(7,14), O(8,16), P(15,31), Ca(20,40)

  double ZZ[num_material], AA[num_material];


     
  if(strcmp(mat_type,"watercylinder")==0)
    {
     
     ZZ[0] = 1; AA[0] = 1;//H
     ZZ[1] = 8; AA[1] = 16;//O

    }

     
     
  else
    {
  
  ZZ[0] = 1; AA[0] = 1;
  ZZ[1] = 6; AA[1] = 12;
  ZZ[2] = 7; AA[2] = 14;
  ZZ[3] = 8; AA[3] = 16;
  ZZ[4] = 15; AA[4] = 31;
  ZZ[5] = 20; AA[5] = 40;
  if(strcmp(mat_type,"gammex")==0)
    {
      //ZZ[4] = 9; AA[4] = 19;
      //ZZ[5] = 18; AA[5] = 40;
      ZZ[6] = 12; AA[6] = 24.31;
      ZZ[7] = 14; AA[7] = 28.09;
    }
  else if(strcmp(mat_type,"catphan")==0)
    {
      ZZ[4] = 9; AA[4] = 19;
      ZZ[5] = 18; AA[5] = 40;
    }


    }


  

  //X-ray data base version!!!!!
  //int allene_db=23001;//5keVstart(interval_5eV) 5.00keV->5.005keV->5.010keV..to120keV
  int allene_db=29001;//5keVstart(interval_5eV) 5.00keV->5.005keV->5.010keV...to150keV
 
  //double*pairdb;
  double* energydb = (double*)malloc( allene_db*6*sizeof(double) );
  double* photodb = (double*)malloc( allene_db*6*sizeof(double) );
  double* comptondb = (double*)malloc( allene_db*6*sizeof(double) );
  double* rayleighdb = (double*)malloc( allene_db*6*sizeof(double) );
  
  FILE*	fpphoto;
  FILE*	fpcomp;
  FILE*	fpraylei;
  FILE*	fppair;
  //if( (fpphoto=fopen("attenuation_xraydb/photo.csv","rb")) == NULL )//photo 5keVto120keV
  if( (fpphoto=fopen("attenuation_xraydb/photo_5keVto150keV.csv","rb")) == NULL )//photo
   {
     fprintf(stderr,"data file: xraydb not found\n");
     return;
    }
  //if( (fpcomp=fopen("attenuation_xraydb/comp.csv","rb")) == NULL )//comp  5keVto120keV
  if( (fpcomp=fopen("attenuation_xraydb/comp_5keVto150keV.csv","rb")) == NULL )//comp 
   {
     fprintf(stderr,"data file: xraydb not found\n");
     return;
    }
  //if( (fpraylei=fopen("attenuation_xraydb/rayl.csv","rb")) == NULL )// rayl 5keVto120keV
  if( (fpraylei=fopen("attenuation_xraydb/rayl_5keVto150keV.csv","rb")) == NULL )// rayl
   {
     fprintf(stderr,"data file: xraydb not found\n");
     return;
    }
  
  //if( (fppair=fopen("pair.csv","rb")) == NULL )//pair
   //{
     //fprintf(stderr,"data file: xraydb not found\n");
     //return;
    //}
  
  static char	lineBuf[LINEBUF_SIZE];
 
  fgets(lineBuf,LINEBUF_SIZE,fpphoto);//skip header H,C,N,O,P,Ca,energy
  fgets(lineBuf,LINEBUF_SIZE,fpcomp);//skip header H,C,N,O,P,Ca,energy
  fgets(lineBuf,LINEBUF_SIZE,fpraylei);//skip header H,C,N,O,P,Ca,energy
  //fgets(lineBuf,LINEBUF_SIZE,fppair);//skip header H,C,N,O,P,Ca,energy

  printf("%s\n", lineBuf);
  
  double energy_db;
  double photodb_H,photodb_C,photodb_N,photodb_O,photodb_P,photodb_Ca;
  double comptondb_H,comptondb_C,comptondb_N,comptondb_O,comptondb_P,comptondb_Ca;
  double rayleighdb_H,rayleighdb_C,rayleighdb_N,rayleighdb_O,rayleighdb_P,rayleighdb_Ca;
  //double pairdb_H,pairdb_C,pairdb_N,pairdb_O,pairdb_P,pairdb_Ca;
  
  for(int enei=0;enei<allene_db;enei++)//photoelectron
      {
	//printf("enei number %d\n", enei);
	fgets(lineBuf,LINEBUF_SIZE,fpphoto);
       	//printf(" %s\n", lineBuf);
	sscanf(lineBuf,"%lf,%lf,%lf,%lf,%lf,%lf,%lf",&photodb_H,&photodb_C,&photodb_N,&photodb_O,&photodb_P,&photodb_Ca,&energy_db);


	energydb[enei]=energy_db/1000/1000;//MeV
	energydb[1*allene_db+enei]=energy_db/1000/1000;//MeV
	energydb[2*allene_db+enei]=energy_db/1000/1000;//MeV
	energydb[3*allene_db+enei]=energy_db/1000/1000;//MeV
	energydb[4*allene_db+enei]=energy_db/1000/1000;//MeV
	energydb[5*allene_db+enei]=energy_db/1000/1000;//MeV
	
	photodb[enei]=photodb_H;
	photodb[1*allene_db+enei]=photodb_C;
	photodb[2*allene_db+enei]=photodb_N;
	photodb[3*allene_db+enei]=photodb_O;
	photodb[4*allene_db+enei]=photodb_P;
	photodb[5*allene_db+enei]=photodb_Ca;
	
      }
 
    for(int enei=0;enei<allene_db;enei++)//compton
      {
	fgets(lineBuf,LINEBUF_SIZE,fpcomp);
	sscanf(lineBuf,"%lf,%lf,%lf,%lf,%lf,%lf,%lf",&comptondb_H,&comptondb_C,&comptondb_N,&comptondb_O,&comptondb_P,&comptondb_Ca,&energy_db);
      
	comptondb[enei]=comptondb_H;
	comptondb[1*allene_db+enei]=comptondb_C;
	comptondb[2*allene_db+enei]=comptondb_N;
	comptondb[3*allene_db+enei]=comptondb_O;
	comptondb[4*allene_db+enei]=comptondb_P;
	comptondb[5*allene_db+enei]=comptondb_Ca;
      }

    for(int enei=0;enei<allene_db;enei++)//rayleigh
      {
	fgets(lineBuf,LINEBUF_SIZE,fpraylei);
	sscanf(lineBuf,"%lf,%lf,%lf,%lf,%lf,%lf,%lf",&rayleighdb_H,&rayleighdb_C,&rayleighdb_N,&rayleighdb_O,&rayleighdb_P,&rayleighdb_Ca,&energy_db);
      
	rayleighdb[enei]=rayleighdb_H;
	rayleighdb[1*allene_db+enei]=rayleighdb_C;
	rayleighdb[2*allene_db+enei]=rayleighdb_N;
	rayleighdb[3*allene_db+enei]=rayleighdb_O;
	rayleighdb[4*allene_db+enei]=rayleighdb_P;
	rayleighdb[5*allene_db+enei]=rayleighdb_Ca;
      }
    
   

    
    
  
  double E, Z, A, Ek=0;
  double photo, compton, pair;
  double rayleigh;
  for(int iz = 0; iz < num_material; iz++)
    {
      printf("iz number %d\n", iz);
      Z=ZZ[iz]; A=AA[iz];
      double fact = Avogadro*Z/A*(2.81794*2.81794)*0.001;//Avogadro=6.0221  //0.001 is  Avogadro(6.02x10^??) *r0*r0  ;r0(2.817x10^??) 
      //printf("Z = %lf, A = %lf \n", Z, A);
      //printf("Energy[MeV], photo[cm^2/g], compton[cm^2/g], pair[cm^2/g], Attenuation[cm^2/g], pdf_pe[i]\n");
      for(int i = 0; i < NNE; i++)
	{
	  
	  

	  /*
	  //Yao formula
	 //Yao formula
	    //Yao formula
	    //E=0.01*i/0.511;
	    E=photonenergy[i]/0.511;//E[MeV]
	    photo = fact*0.00000345*(1.0+0.008*Z)*(Z*Z*Z)/(E*E*E)
	      *(1.0-Ek/(4.0*E)-Ek*Ek/(1.21*E));
	    compton = fact*2*PI*((1.0+E)/(E*E)*(2.0*(1.0+E)/(1.0+2*E)
				 -log(1.0+2*E)/E)
				 +log(1.0+2*E)/(2*E)
				 -(1.0+3*E)/((1.0+2*E)*(1.0+2*E)));
	    pair = fact*0.2545*(E-2.332)*Z/137;
	    if(pair < 0) pair = 0;
	 
	 
	    //printf("energy%lf photo%lf, compton%lf\n", photonenergy[i], photo, compton);
	    //exit(0);

	   
	    */

	  

	  
	 
	     //X-ray Data Base attenuation
	      //X-ray Data Base attenuation
	      //X-ray Data Base attenuation
	   for(int enei=0;enei<allene_db;  enei++)//forx-ray DB
	    {
	    E=photonenergy[i];//E,photonenergy[i]=spectre energy [MeV]
	    double aaa=energydb[iz*allene_db+enei]*1000*1000*1000-photonenergy[i]*1000*1000*1000;//(X-raydb -spectre)*1000*1000*1000,not 0,approximately 1
	    
	   
	    //if(energydb[iz*allene_db+enei]==E)//this is not matching because of float calculation of spectre(float vs short)
	    if(-30<aaa && aaa<30)//aaa is approximately -1 - 1, when  the energy is matching between "spectre" and "xraydb"
	    //if(-5<aaa && aaa<5)//aaa is approximately -1 - 1, when  the energy is matching between "spectre" and "xraydb"
	      {
		printf("energy %lf %lf \n",photonenergy[i]*1000*1000*1000,energydb[iz*allene_db+enei]*1000*1000*1000);
		photo =photodb[iz*allene_db+enei];
		compton =comptondb[iz*allene_db+enei];
		rayleigh=rayleighdb[iz*allene_db+enei];
		pair =0 ;
		//pair = fact*0.2545*(E-2.332)*Z/137;
		if(pair < 0) pair = 0;
		
		//printf("energy%lf photo%lf, compton%lf\n", E, photo, compton);
		//exit(0);
	      }
	 
	    }//enei forx-ray DB

	   
	  

	   
	  
	   //Attenuation[iz*NNE+i] = (photo + compton + pair);//original Yao
	   Attenuation[iz*NNE+i] = (photo + compton +rayleigh + pair);//X-ray data base 2022/1/5
	    
	    //printf("%lf %lf %lf %lf %lf %f\n",E*0.511, photo, compton, pair, Attenuation[NE*iz+i], pdf_pe[i]);

	   
	    
	}
    }//iz
  
  *NE = NNE;

  free(photodb);
  free(comptondb);
  free(rayleighdb);
  free(energydb);

  
  fclose(fpphoto); 
  fclose(fpcomp);
  fclose(fpraylei);
  //fclose(fppair);
  
}


