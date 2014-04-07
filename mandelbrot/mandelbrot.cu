/** @file histo-global.cu histogram with global memory atomics */

#include <png.h>
#include <omp.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <time.h>

/** CUDA check macro */
#define cucheck(call) \
	{\
	cudaError_t res = (call);\
	if(res != cudaSuccess) {\
	const char* err_str = cudaGetErrorString(res);\
	fprintf(stderr, "%s (%d): %s in %s", __FILE__, __LINE__, err_str, #call);	\
	exit(-1);\
	}\
	}

/** time spent in device */
double gpu_time = 0;

/** a useful function to compute the number of threads */
int divup(int x, int y) { return x / y + (x % y ? 1 : 0); }

/** gets the color, given the dwell */
void dwell_color(int *r, int *g, int *b, int dwell);

/** save the dwell into a PNG file 
		@remarks: code to save PNG file taken from here 
		  (error handling is removed):
		http://www.labbookpages.co.uk/software/imgProc/libPNG.html
 */
void save_image(const char *filename, int *dwells, int w, int h) {
	png_bytep row;
	
	FILE *fp = fopen(filename, "wb");
	png_structp png_ptr = png_create_write_struct(PNG_LIBPNG_VER_STRING, 0, 0, 0);
	png_infop info_ptr = png_create_info_struct(png_ptr);
	// exception handling
	setjmp(png_jmpbuf(png_ptr));
	png_init_io(png_ptr, fp);
	// write header (8 bit colour depth)
	png_set_IHDR(png_ptr, info_ptr, w, h,
							 8, PNG_COLOR_TYPE_RGB, PNG_INTERLACE_NONE,
							 PNG_COMPRESSION_TYPE_BASE, PNG_FILTER_TYPE_BASE);
	// set title
	png_text title_text;
	title_text.compression = PNG_TEXT_COMPRESSION_NONE;
	title_text.key = "Title";
	title_text.text = "Mandelbrot set, per-pixel";
	png_set_text(png_ptr, info_ptr, &title_text, 1);
	png_write_info(png_ptr, info_ptr);

	// write image data
	row = (png_bytep) malloc(3 * w * sizeof(png_byte));
	for (int y = 0; y < h; y++) {
		for (int x = 0; x < w; x++) {
			int r, g, b;
			dwell_color(&r, &g, &b, dwells[y * w + x]);
			row[3 * x + 0] = (png_byte)r;
			row[3 * x + 1] = (png_byte)g;
			row[3 * x + 2] = (png_byte)b;
		}
		png_write_row(png_ptr, row);
	}
	png_write_end(png_ptr, NULL);

  fclose(fp);
  png_free_data(png_ptr, info_ptr, PNG_FREE_ALL, -1);
  png_destroy_write_struct(&png_ptr, (png_infopp)NULL);
  free(row);
}  // save_image


/** a simple complex type */
struct complex {
	__host__ __device__ complex(float re, float im = 0) {
		this->re = re;
		this->im = im;
	}
	/** real and imaginary part */
	float re, im;
}; // struct complex

// operator overloads for complex numbers
inline __host__ __device__ complex operator+
(const complex &a, const complex &b) {
	return complex(a.re + b.re, a.im + b.im);
}
inline __host__ __device__ complex operator-
(const complex &a) { return complex(-a.re, -a.im); }
inline __host__ __device__ complex operator-
(const complex &a, const complex &b) {
	return complex(a.re - b.re, a.im - b.im);
}
inline __host__ __device__ complex operator*
(const complex &a, const complex &b) {
	return complex(a.re * b.re - a.im * b.im, a.im * b.re + a.re * b.im);
}
inline __host__ __device__ float abs2(const complex &a) {
	return a.re * a.re + a.im * a.im;
}
inline __host__ __device__ complex operator/
(const complex &a, const complex &b) {
	float invabs2 = 1 / abs2(b);
	return complex((a.re * b.re + a.im * b.im) * invabs2,
								 (a.im * b.re - b.im * a.re) * invabs2);
}  // operator/

#define MAX_DWELL 256
#define BS 256
/** computes the dwell for a single pixel */
__device__ int pixel_dwell
(int w, int h, complex cmin, complex cmax, int x, int y) {
	complex dc = cmax - cmin;
	float fx = (float)x / w, fy = (float)y / h;
	complex c = cmin + complex(fx * dc.re, fy * dc.im);
	int dwell = 0;
	complex z = c;
	while(dwell < MAX_DWELL && abs2(z) < 2 * 2) {
		z = z * z + c;
		dwell++;
	}
	return dwell;
}  // pixel_dwell

/** computes the dwells for Mandelbrot image 
		@param dwells the output array
		@param w the width of the output image
		@param h the height of the output image
		@param cmin the complex value associated with the left-bottom corner of the
		image
		@param cmax the complex value associated with the right-top corner of the
		image
 */
__global__ void mandelbrot_k
(int *dwells, int w, int h, complex cmin, complex cmax) {
	// complex value to start iteration (c)
	int x = threadIdx.x + blockIdx.x * blockDim.x;
	int y = threadIdx.y + blockIdx.y * blockDim.y;
	int dwell = pixel_dwell(w, h, cmin, cmax, x, y);
	dwells[y * w + x] = dwell;
}  // mandelbrot_k

/** gets the color, given the dwell (on host) */
#define CUT_DWELL (MAX_DWELL / 4)
void dwell_color(int *r, int *g, int *b, int dwell) {
	// black for the Mandelbrot set
	if(dwell >= MAX_DWELL) {
		*r = *g = *b = 0;
	} else {
		// cut at zero
		if(dwell < 0)
			dwell = 0;
		if(dwell <= CUT_DWELL) {
			// from black to blue the first half
			*r = *g = 0;
			*b = 128 + dwell * 127 / (CUT_DWELL);
		} else {
			// from blue to white for the second half
			*b = 255;
			*r = *g = (dwell - CUT_DWELL) * 255 / (MAX_DWELL - CUT_DWELL);
		}
	}
}  // dwell_color

/** data size */
#define H (8 * 1024)
#define W (8 * 1024)
#define IMAGE_PATH "./mandelbrot.png"

int main(int argc, char **argv) {
	// allocate memory
	int w = W, h = H;
	size_t dwell_sz = w * h * sizeof(int);
	int *h_dwells, *d_dwells;
	cucheck(cudaMalloc((void**)&d_dwells, dwell_sz));
	h_dwells = (int*)malloc(dwell_sz);

	// compute the dwells, copy them back
	double t1 = omp_get_wtime();
	dim3 bs(64, 4), grid(divup(w, bs.x), divup(h, bs.y));
	mandelbrot_k<<<grid, bs>>>
		(d_dwells, w, h, complex(-1.5, -1), complex(0.5, 1));
	cucheck(cudaThreadSynchronize());
	double t2 = omp_get_wtime();
	cucheck(cudaMemcpy(h_dwells, d_dwells, dwell_sz, cudaMemcpyDeviceToHost));
	gpu_time = t2 - t1;
	
	// save the image to PNG 
	save_image(IMAGE_PATH, h_dwells, w, h);

	// print performance
	printf("Mandelbrot set computed in %.3lf s, at %.3lf Mpix/s\n", gpu_time, 
				 h * w * 1e-6 / gpu_time);

	// free data
	cudaFree(d_dwells);
	free(h_dwells);
	return 0;
}  // main
