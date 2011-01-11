/* defaults */
#define DXRES 3500
#define DYRES 2750
#define DXMIN -1.272727
#define DXMAX 1.272727
#define DNUMV 16
#define DYMIN -1.0
#define DYMAX 1.0
#define DSAMPLES 20000
#define DITT 1000
#define DQUALITY 100
#define DSUPER 1
#define DGAMMA 2.2
#define PATH "/tmp/fractal.tif"

/* structs */

typedef struct { /* color content of a pixel: R G B channels" */
	unsigned char r, g, b;
	} trio;

typedef union { 
	unsigned int count;
	float normal;
	} values;

typedef struct { /* set of arrays of coefficients a-f */
	double *ac, *bc, *cc, *dc, *ec, *fc;
	trio *colors;
	} coeff;

typedef struct { /* basic contents of a pixel: color, hit counter, normalized value */
	values vals;
	trio color;
	} pixel;

typedef struct { /* fractal: pixel resolutions, x and y coordinate ranges, pixel array */
	int xres, yres;
	double xmax, ymax, ymin, xmin, ranx, rany;
	pixel** pixels;
	} fract;

/* prototypes */
double modulus( double, double);
void V(int k,  double x,  double y, int i,  double* ret_x,  double* ret_y, coeff* co);
unsigned char channel();
double sign();
double choose( double a,  double b);
void fractal_set(int xres, int yres,  double xmin,  double xmax,  double ymin,  double ymax, fract* fractal);
void reduce(fract* fractal, int sample);
void print_usage();
FILE* my_file_open(char* filename);
void write_to_jpeg(fract* fractal, int quality, char* filename, int invert);
void write_to_tiff(fract* fractal, char* filename, int invert);
void coeffs_init(coeff *coeffs, int NUMV, int R, int G, int B);
