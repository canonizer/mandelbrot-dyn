NAME=mandelbrot-dyn
SRC=$(NAME).cu
TGT=$(NAME)
TMP=*~ $(TGT)

build:	$(TGT)
$(TGT):	$(SRC) makefile
	nvcc -O3 -arch=sm_35 -rdc=true -lcudadevrt -Xcompiler -fopenmp -lpng $(SRC) -o $(TGT)

run:	$(TGT)
	./$(TGT)

clean:
	rm -f $(TMP)
