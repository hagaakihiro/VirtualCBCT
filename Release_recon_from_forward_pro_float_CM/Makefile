
TARGET := cbct2
CC := g++ -no-pie
#CC := nvcc
#CC := icc

$(TARGET):		main.o	\
			mallocMD.o \
			misc.o \
			projectionDataMVkV.o \
			filterSinogram.o \
			reconstructImage.o
	$(CC)	main.o	\
		mallocMD.o \
		misc.o \
		projectionDataMVkV.o \
		filterSinogram.o \
		reconstructImage.o \
		-fopenmp -O2 -lm \
		-o $(TARGET)


main.o :		main.cpp \
			physParams.h \
			mallocMD.h \
			misc.h \
			projectionDataMVkV.h \
			filterSinogram.h \
			reconstructImage.h
	$(CC)	main.cpp \
		-fopenmp -O2 -c

mallocMD.o:		mallocMD.cpp \
			mallocMD.h
	$(CC)	mallocMD.cpp \
		-fopenmp -O2 -c

misc.o:		misc.cpp \
			misc.h \
			mallocMD.h
	$(CC)	misc.cpp \
		-fopenmp -O2 -c

projectionDataMVkV.o:	projectionDataMVkV.cpp \
			projectionDataMVkV.h \
			physParams.h \
			mallocMD.h
	$(CC)	projectionDataMVkV.cpp \
		-fopenmp -O2 -c

filterSinogram.o:	filterSinogram.cpp \
			filterSinogram.h \
			physParams.h \
			mallocMD.h
	$(CC)	filterSinogram.cpp \
		-fopenmp -O2 -c

reconstructImage.o:	reconstructImage.cpp \
			reconstructImage.h \
			physParams.h \
			mallocMD.h
	$(CC)	reconstructImage.cpp \
		-fopenmp -O2 -c

clean:
	rm -f *.o $(TARGET)
