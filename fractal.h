/*
 * This code may be freely redistributed unter the
 * terms of the GPL
 *
 * James McCarty 
 *
 * 2012
 */

/* default values */
#define XRES 1920
#define YRES 1080
#define XMIN -1.777
#define XMAX 1.777
#define NUMV 16
#define YMIN -1.0
#define YMAX 1.0
#define SAMPLES 20000
#define ITT 1000
#define SUPER 1
#define GAMMA 2.2
#define SEED 1
#define PATH "/tmp/fractal.tif"
#define THREADGROUPSIZE 9

/* structs */

typedef struct {
	double ac, bc, cc, dc, ec, fc; /* set of coefficients a-f */
	double pa1, pa2, pa3, pa4;
	unsigned char r, g, b; /* color content of a pixel: RBG channels */
	} coeff;

typedef union {
	unsigned int counter; /* number of times pixel has been incremented */
	float normal; /* normalized value at pixel */
	} hitcounter;

typedef struct {
	hitcounter value; 
	unsigned char r, g, b; /* color content of a pixel: RGB channels */
	} pixel;

typedef struct { /* all components of a flame fractal */
	int xres, yres; /* x and y resolution of image */
	double xmin,ymin,xmax,ymax; /* axis bounds */
	double ranx, rany; /* numerical range of x/y axis */
	int R,G,B;		/* fixed color channels */
	int n;          /* number of equations */
	int sup;		/* super sample value  */
	int samples;	/* number of flame samples */
	long int iterations; /* number of iterations per sample */
	int invert;		/* use inverse colors? 0 false, else true */
	int symmetry;	/* use symmetrical rotation axis? set to greater than 1 */
	int seed;		/* random seed */
  int num_threads; /* number of threads to use in render */
	double gamma;	/* gamma correction factor */
	char *file;		/* output file path */
	coeff *coarray; /* array of coefficients */
	pixel **pixels; /* image buffer */
	int *choice;	/* transformations to use */
	int count;		/* number of tranformations available */
	FILE *pallete; /* file containing RGB values for coeff array */
	FILE *cofile;  /* file containing coefficients for coeff array */
  pthread_mutex_t *lock; /* lock so that only one section of memory buffer being written at a time */
	} flame;

/* prototypes */
int random_bit(void);
double modulus(double a, double b);
void ContractiveMapping(coeff *coeff);
void print_usage();
void reduce(flame *fractal);
void fractal_init(flame* fractal);
void coeff_init(flame* fractal);
void gamma_log(flame* fractal);
void parse_args(int argc, char **argv, flame* fractal);
void buffer_init(flame* fractal);
void write_to_tiff(flame* fractal);
void *render(void *fract);
