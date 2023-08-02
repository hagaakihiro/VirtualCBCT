
TARGET := proj.exe
gC := g++
CC := nvcc -D_FORCE_INLINES -arch=sm_35
#CC := nvcc -D_FORCE_INLINES -arch=sm_75
#CC := nvcc
#CC := icc

$(TARGET):		main.o	\
			mallocMD.o \
                        misc.o \
			projectionData.o \
			reprojection.o \
			IR_ImageCBCT.o \
			prior_weight_production.o \
			lengthproj.o \
			reprojectionOnGPU.o 

	$(gC)	main.o	\
		mallocMD.o \
                misc.o \
		projectionData.o \
		reprojection.o \
		IR_ImageCBCT.o \
		prior_weight_production.o \
		lengthproj.o \
		reprojectionOnGPU.o \
		-fopenmp -O2 -lm -lcudart\
		-o $(TARGET)



main.o :		main.cpp \
			physParams.h \
			mallocMD.h \
                        misc.h \
			projectionData.h \
			reprojection.h \
			IR_ImageCBCT.h \
			virtual_projection.h \
			lengthproj.h \
			reprojectionOnGPU.h
	$(gC)	main.cpp \
		-fopenmp -O2 -c

mallocMD.o:		mallocMD.cpp \
			mallocMD.h
	$(gC)	mallocMD.cpp \
		-fopenmp -O2 -c

misc.o:		misc.cpp \
	        mallocMD.h \
		misc.h \
		physParams.h 
	$(gC)	misc.cpp \
		-fopenmp -O2 -c

projectionData.o:	projectionData.cpp \
			projectionData.h \
			physParams.h \
			mallocMD.h
	$(gC)	projectionData.cpp \
		-fopenmp -O2 -c

reprojection.o:	reprojection.cpp \
			physParams.h \
			mallocMD.h \
			reprojection.h
	$(gC)	reprojection.cpp \
		-fopenmp -O2 -c

IR_ImageCBCT.o:	IR_ImageCBCT.cpp \
			physParams.h \
			mallocMD.h \
			IR_ImageCBCT.h
	$(gC)	IR_ImageCBCT.cpp \
		-fopenmp -O2 -c

lengthproj.o: lengthproj.cpp \
		physParams.h \
		mallocMD.h \
		lengthproj.h
	$(gC)	lengthproj.cpp \
		-fopenmp -O2 -c

prior_weight_production.o: prior_weight_production.cpp \
		physParams.h \
		mallocMD.h \
		virtual_projection.h
	$(gC)	prior_weight_production.cpp \
		-fopenmp -O2 -c

reprojectionOnGPU.o:	reprojectionOnGPU.cu \
			physParams.h \
			mallocMD.h \
			reprojectionOnGPU.h
	$(CC)	reprojectionOnGPU.cu \
		-O2 -c

clean:
	rm -f *.o $(TARGET)
