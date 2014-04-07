Mandelbrot fractal with CUDA dynamic parallelism (sample)


INTRO

This is a set of two samples that both create an image of the Mandelbrot fractal
on NVidia GPUs with and without dynamic parallelism. Each sample has its own
directory:

mandelbrot/                 - without dynamic parallelism
mandelbrot-dyn/             - with dynamic parallelism

The sample without dynamic parallelism uses escape time algorithm
(http://en.wikipedia.org/wiki/Mandelbrot_fractal#Escape_time_algorithm) to
compute the dwell on the GPU. The dynamic parallelism sample uses the
hierarchical Mariani-Silver algorithm (see, e.g.,
http://mrob.com/pub/muency/marianisilveralgorithm.html). In both cases, final
coloring is done on the CPU, just to make the code clearer. Depending on the
image size and maximum dwell, the dynamic parallelism version is 1.3-6x faster
than the version where each pixel is evaluated.

This code has been developed by Andrew V. Adinetz, who is employed as researcher
in Juelich Supercomputing Centre, Forschungszentrum Juelich. This code is
property of Forschungszentrum Juelich. Users can modify and/or distribute code
free of charge, provided that this notice is retained. 

These samples have been primarily developed for educational purposes. They have
been used, e.g., in a lecture on CUDA dynamic parallelism in Advanced GPU
Programming course in Forschungszentrum Juelich
(http://www.fz-juelich.de/ias/jsc/EN/Expertise/Services/Documentation/presentations/presentation-adv-gpu_table.html). Other
users can freely use this code for teaching, illustrative or copy-paste
purposes, if the link to the original source is provided.


REQUIREMENTS

For both: 

- Unix-like OS (tested on Linux)
- libpng (to save the generated image)

For Mandelbrot without dynamic parallelism:

- CUDA 4.0 or higher 
- NVidia GPU with Compute Capability 2.0 or higher (= Fermi or higher)

For Mandelbrot with dynamic parallelism:

- CUDA 5.0 or higher
- NVidia GPU with Compute Capability 3.5 or higher (= Kepler K20 or higher)


COMPILING AND RUNNING

The following commands should be executed in the directory of the corresponding
sample. 

To compile a sample, execute the following command:

	 make

To compile and run a sample, execute the following command:

	 make run

Upon running successfully, the sample will create the PNG image of Mandelbrot
fractal (mandelbrot.png). The resolution of the image and the maximum dwell used
to generate it are set in the source code of the sample using W, H and MAX_DWELL
constants, respectively. 
