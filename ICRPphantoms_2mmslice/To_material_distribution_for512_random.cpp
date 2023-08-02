//Transform into Material Distribution ICRP110 ver1.2
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <stdio.h>
#include <string.h>
#include <math.h>
//#include <cmath>
//#include <complex>
//#include <vector>
//#include <iostream>


double Uniform0( void )
{
    return ((double)rand()+1.0)/((double)RAND_MAX+2.0);
}

double rand_gauss0( double mu, double sigma )
{
    double PI = 3.14159265;
    double z = sqrt( -2.0*log(Uniform0()) ) * sin( 2.0*PI*Uniform0() );
    return mu + sigma*z;
}



int
main(int argc, char* argv[])
{
    int igender = 0; // =0; Male, =1; Female
    
    char filename[128], filename2[128], filename3[128];
    if(igender == 0)
    {
        //char *str2 = "AM_rescale_512x512_222_mod.raw";
        char *str2 = "AM_skin_543x271x879.raw";
        sprintf(filename, "%s", str2);
    }
    else
    {
        //char *str2 = "AF_rescale_512x512_348_mod.raw";
        char *str2 = "AF_skin_531x243x838.raw";
        sprintf(filename, "%s", str2);
    }
    
    printf("%s \n", filename);

    FILE    *fp, *fp2;
    if((fp=fopen(filename,"rb")) == NULL) {
        printf("OPEN FAILED ***.raw \n");
        exit(0);
    }
    
    // for Male
    int ncol_M=543;
    int nrow_M=271;
    int nsli_M=879;
    // for Female
    int ncol_F=531;
    int nrow_F=243;
    int nsli_F=838;
    
    short  *outputimg0, *outputimg2;
    unsigned char *outputimg1;
    float  *outputimg3;
    int raw_size, nsize, raw_size0;
    if(igender == 0)
    {
        raw_size = ncol_M*nrow_M*nsli_M*sizeof(short);
        raw_size0 = ncol_M*nrow_M*nsli_M*sizeof(unsigned char);
        nsize = ncol_M*nrow_M*nsli_M;
    }
    else
    {
        raw_size = ncol_F*nrow_F*nsli_F*sizeof(short);
        raw_size0 = ncol_F*nrow_F*nsli_F*sizeof(unsigned char);
        nsize = ncol_F*nrow_F*nsli_F;
    }
    //outputimg1 = (short*)malloc(raw_size);
    outputimg1 = (unsigned char*)malloc(raw_size);
    fread(outputimg1, 1, raw_size0, fp);
    fclose(fp);
    //  read raw data -- end
    
    if(igender == 0)
    {
        char *str2 = "AM_media.csv";
        sprintf(filename, "%s", str2);
        char *str3 = "AM_organs.csv";
        sprintf(filename2, "%s", str3);
    }
    else
    {
        char *str2 = "AF_media.csv";
        sprintf(filename, "%s", str2);
        char *str3 = "AF_organs.csv";
        sprintf(filename2, "%s", str3);
    }
    if((fp=fopen(filename,"r")) == NULL) {
        printf("OPEN FAILED ***.dat \n");
        exit(0);
    }
    
    printf("%s \n", filename);
    
    int i = 0;
    float Hw[100],Cw[100],Nw[100],Ow[100],Naw[100],Maw[100],
    Pw[100],Sw[100],Clw[100],Kw[100],Caw[100],Few[100],Iw[100];
    int ID[100];
    double sigma, mu;

    while(fscanf(fp,"%d,%f,%f,%f,%f,%f,%f,%f,%f,%f,%f,%f,%f,%f",&ID[i],&Hw[i],&Cw[i],&Nw[i],&Ow[i],&Naw[i],&Maw[i],&Pw[i],&Sw[i],&Clw[i],&Kw[i],&Caw[i],&Few[i],&Iw[i]) != EOF)
    {
        if(ID[i]!=53)
        {
            sigma = double(Hw[i]*0.1);
            mu = double(Hw[i]);
            Hw[i] = (float)(rand_gauss0( mu, sigma ));
            if(Hw[i]<0.0) Hw[i] = 0.0;
            
            sigma = double(Cw[i]*0.1);
            mu = double(Cw[i]);
            Cw[i] = (float)(rand_gauss0( mu, sigma ));
            if(Cw[i]<0.0) Cw[i] = 0.0;
            
            sigma = double(Nw[i]*0.1);
            mu = double(Nw[i]);
            Nw[i] = (float)(rand_gauss0( mu, sigma ));
            if(Nw[i]<0.0) Nw[i] = 0.0;
            
            if(Pw[i]>0.09)
            {
                sigma = double(Pw[i]*0.1);
                mu = double(Pw[i]);
                Pw[i] = (float)(rand_gauss0( mu, sigma ));
                if(Pw[i]<0.0) Pw[i] = 0.0;
            }
            if(Caw[i]>0.09)
            {
                sigma = double(Caw[i]*0.1);
                mu = double(Caw[i]);
                Caw[i] = (float)(rand_gauss0( mu, sigma ));
                if(Caw[i]<0.0) Caw[i] = 0.0;
            }
        }
        Ow[i] = 100 - Hw[i] - Cw[i] - Nw[i] - Pw[i] - Caw[i];
        if(Ow[i]<0.0)
        {
            double s = (Hw[i] + Cw[i] + Nw[i] + Pw[i] + Caw[i])*0.01;
            Ow[i]=0.0;
            Hw[i]=Hw[i]/s;
            Cw[i]=Cw[i]/s;
            Nw[i]=Nw[i]/s;
            Pw[i]=Pw[i]/s;
            Caw[i]=Caw[i]/s;
            
        }
        printf("%d %f %f %f %f %f %f\n", ID[i], Hw[i], Cw[i], Nw[i], Ow[i], Pw[i], Caw[i]);
        i++;
    }
    fclose(fp);
    //exit(0);
    if((fp=fopen(filename2,"r")) == NULL) {
        printf("OPEN FAILED ***.dat \n");
        exit(0);
    }
    
    printf("%s \n", filename2);
    int Organ[150], OrganID[150];
    float Dens[150];
    float dens;
    i = 0;
    while(fscanf(fp,"%d,%d,%f",&Organ[i],&OrganID[i],&Dens[i]) != EOF)
    {
        if(OrganID[i]==53)
        {
            dens = Dens[i];
        }
        else if(OrganID[i]<=2)
        {
            sigma = double(Dens[i]*0.05);
            mu = double(Dens[i]*0.9);
            dens = (float)(rand_gauss0( mu, sigma ));
        }
        else
        {
            sigma = double(Dens[i]*0.05);
            mu = double(Dens[i]);
            dens = (float)(rand_gauss0( mu, sigma ));
        }
        printf("%d %d %d %f %f\n", i, Organ[i], OrganID[i], Dens[i], dens);
        i++;
    }
    int imax = i;
    printf("imax = %d \n",i);
    fclose(fp);
    //  read organ and material data -- end
    //exit(0);

    
    if((fp=fopen("test.raw","wb")) == NULL) {
        printf("OPEN FAILED ***.raw \n");
        exit(0);
    }
    outputimg2 = (short*)malloc(raw_size*6);
    // Material includes H, C, N, O, P, and Ca
    for(int n=0; n<nsize; n++)
    {
        for(int j=0; j<imax; j++)
        {
            if(Organ[j]==outputimg1[n])
            {
                outputimg2[n]=(short)(Hw[OrganID[j]-1]*Dens[j]*10);
                outputimg2[n+nsize]=(short)(Cw[OrganID[j]-1]*Dens[j]*10);
                outputimg2[n+nsize*2]=(short)(Nw[OrganID[j]-1]*Dens[j]*10);
                outputimg2[n+nsize*3]=(short)(Ow[OrganID[j]-1]*Dens[j]*10);
                outputimg2[n+nsize*4]=(short)(Pw[OrganID[j]-1]*Dens[j]*10);
                outputimg2[n+nsize*5]=(short)(Caw[OrganID[j]-1]*Dens[j]*10);
            }
        }
    }
    fwrite(outputimg2, 1, raw_size*6, fp);
    fclose(fp);
    
    double AA[6], ZZ[6];
    AA[0] = 1.0; ZZ[0] = 1.0;
    AA[1] = 12.0; ZZ[1] = 6.0;
    AA[2] = 14.0; ZZ[2] = 7.0;
    AA[3] = 16.0; ZZ[3] = 8.0;
    AA[4] = 31.0; ZZ[4] = 15.0;
    AA[5] = 40.0; ZZ[5] = 20.0;
    
    outputimg0 = (short*)malloc(ncol_M*nrow_M*6*sizeof(short));
    outputimg3 = (float*)malloc(ncol_M*nrow_M*sizeof(float));
    if(igender == 0)
    {
        for(int iz = 0; iz < nsli_M; iz++)
        {
            char *str = "AM_MD/AM_MD";
            sprintf(filename2, "%s_%d_%d_%03d.raw", str, ncol_M, nrow_M, iz);
            char *str2 = "AM_ED/AM_ED";
            sprintf(filename3, "%s_%d_%d_%03d.raw", str2, ncol_M, nrow_M, iz);
            
            if((fp=fopen(filename2,"wb")) == NULL) {
                printf("OPEN FAILED ***.raw \n");
                exit(0);
            }
            if((fp2=fopen(filename3,"wb")) == NULL) {
                printf("OPEN FAILED ***.raw \n");
                exit(0);
            }
            int k = 0;
            for(int m = 0; m < 6; m++)
            {
                int kk=0;
                for(int i=0; i<ncol_M*nrow_M; i++)
                {
                    int ii = m*ncol_M*nrow_M*nsli_M + iz*ncol_M*nrow_M + i;
                    outputimg0[k] = outputimg2[ii];
                    outputimg3[kk] += outputimg0[k]*1.0;
                    k++;
                    kk++;
                }
            }
            fwrite(outputimg0, 1, ncol_M*nrow_M*6*sizeof(short), fp);
            fwrite(outputimg3, 1, ncol_M*nrow_M*sizeof(float), fp2);
            fclose(fp);
            fclose(fp2);
        }
    }
    else
    {
        for(int iz = 0; iz < nsli_F; iz++)
        {
            char *str = "AF_MD/AF_MD";
            sprintf(filename2, "%s_%d_%d_%03d.raw", str, ncol_F, nrow_F, iz);
            char *str2 = "AF_ED/AF_ED";
            sprintf(filename3, "%s_%d_%d_%03d.raw", str2, ncol_F, nrow_F, iz);
            
            if((fp=fopen(filename2,"wb")) == NULL) {
                printf("OPEN FAILED ***.raw \n");
                exit(0);
            }
            if((fp2=fopen(filename3,"wb")) == NULL) {
                printf("OPEN FAILED ***.raw \n");
                exit(0);
            }
            int k = 0;
            for(int m = 0; m < 6; m++)
            {
                int kk=0;
                for(int i=0; i<ncol_F*nrow_F; i++)
                {
                    int ii = m*ncol_F*nrow_F*nsli_F + iz*ncol_F*nrow_F + i;
                    outputimg0[k] = outputimg2[ii];
                    outputimg0[kk] += outputimg0[k]*1.0;
                    k++;
                    kk++;
                }
            }
            fwrite(outputimg0, 1, ncol_F*nrow_F*6*sizeof(short), fp);
            fwrite(outputimg3, 1, ncol_F*nrow_F*sizeof(float), fp2);
            fclose(fp);
            fclose(fp2);
        }
    }
    
    
    free(outputimg0);
    free(outputimg1);
    free(outputimg2);
    free(outputimg3);
}
    
