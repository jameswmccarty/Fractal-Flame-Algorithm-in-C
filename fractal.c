#include<stdio.h>
#include<stdlib.h>
#include<jpeglib.h>
#include<tiffio.h>
#include<math.h>
#include<malloc.h>
#include<string.h>
#include "fractal.h"

void V(int k,  double x,  double y, int i,  double* ret_x,  double* ret_y, coeff* co)
{

	double a, b, c, d, e, f, r, theta, phi;
	double P0, P1, omega, prefix, t;
	a = (*co).ac[i];
	b = (*co).bc[i];
	c = (*co).cc[i];
	d = (*co).dc[i];
	e = (*co).ec[i];
	f = (*co).fc[i];
	r = sqrt( pow(x, 2.0) + pow (y, 2.0) );
	theta = atan2(y,x);
	phi   = atan2(x,y);

	switch(k){
		case 0: /* Linear */
			*ret_x = x;
			*ret_y = y;
			break;
		case 1: /* Sinusoidal */
			*ret_x = sin(x);
			*ret_y = sin(y);
			break;
		case 2: /* Spherical */
			*ret_x = (1.0 / pow(r, 2.0)) * x;
			*ret_y = (1.0 / pow(r, 2.0)) * y;
			break;
		case 3: /* Swirl */
			*ret_x = x * sin( pow (r, 2.0) ) - y * cos( pow (r, 2.0) );
			*ret_y = x * cos( pow (r, 2.0) ) + y * sin( pow (r, 2.0) );
			break;
		case 4: /* Horseshoe */
			*ret_x = (1.0 / r) * (x - y) * (x + y);
			*ret_y = (1.0 / r) * 2.0 * x * y;
			break;
		case 5: /* Polar */
			*ret_x = theta / M_PI;
			*ret_y = r - 1.0;
			break;
		case 6: /* Handkerchief */
			*ret_x = r * sin( theta + r );
			*ret_y = r * cos( theta - r );
			break;
		case 7: /* Heart */
			*ret_x = r * sin( theta * r );
			*ret_y = r * cos( theta * r ) * (-1.0);
			break;
		case 8: /* Disk */
			*ret_x = (theta / M_PI) * sin(M_PI * r);
			*ret_y = (theta / M_PI) * cos(M_PI * r);
			break;
		case 9: /* Spiral */
			*ret_x = (1.0 / r) * (cos(theta) + sin(r));
			*ret_y = (1.0 / r) * (sin(theta) - cos(r));
			break;
		case 10: /* Hyperbolic */
			*ret_x = sin(theta) / r;
			*ret_y = r * cos(theta);
			break;
		case 11: /* Diamond */
			*ret_x = sin(theta) * cos(r);
			*ret_y = cos(theta) * sin(r);
			break;
		case 12: /* Ex */
			P0 = sin( theta + r );
			P1 = cos( theta - r );
			*ret_x = r * ( pow( P0, 3.0 ) + pow( P1, 3.0 ) );
			*ret_y = r * ( pow( P0, 3.0 ) - pow( P1, 3.0 ) );
			break;
		case 13: /* Julia */
			omega = choose(M_PI, 0.0);
			*ret_x = sqrt(r) * cos((theta / 2.0) + omega);
			*ret_y = sqrt(r) * sin((theta / 2.0) + omega);
			break;
		case 14: /* Bent */
			if( x >= 0.0 && y >= 0.0 ) {
				*ret_x = x;
				*ret_y = y; 
			}
			else if( x < 0.0 && y >= 0.0) {
				*ret_x = 2.0 * x;
				*ret_y = y;
			}
			else if( x >= 0.0 && y < 0.0) {
				*ret_x = x;
				*ret_y = y / 2.0;
			}
			else if( x < 0.0 && y < 0.0) {
				*ret_x = 2.0 * x;
				*ret_y = y / 2.0;
			}
			break;
		case 15: /* Waves */
			*ret_x = x + b * sin( y / pow( c, 2.0) );
			*ret_y = y + e * sin( x / pow( f, 2.0) );
			break;
		case 16: /* Fisheye */
			*ret_x = (2.0 / (r + 1.0)) * y;
			*ret_y = (2.0 / (r + 1.0)) * x;
			break;
		case 17: /* Popcorn */
			*ret_x = x + c*sin(tan(3.0*y));
			*ret_y = y + f*sin(tan(3.0*x));
			break;
		case 18: /* Exponential */
			*ret_x = exp(x - 1.0) * cos(M_PI*y);
			*ret_y = exp(x - 1.0) * sin(M_PI*y);
			break;
		case 19: /* Power */
			*ret_x = pow(r, sin(theta)) * cos(theta);
			*ret_y = pow(r, sin(theta)) * sin(theta);
			break;
		case 20: /* Cosine */
			*ret_x = cos(M_PI*x)*cosh(y);
			*ret_y = -1.0 * sin(M_PI * x) * sinh(y);
			break;
		case 21: /* Rings */
			prefix = modulus((r+c*c), (2.0*c*c))-(c*c)+(r*(1.0-c*c));
			*ret_x = prefix * cos(theta);
			*ret_y = prefix * sin(theta);
			break;
		case 22: /* Fan */
			t = M_PI * pow(c, 2.0);
			if(modulus((theta + f), t) > (t / 2.0)) {
				*ret_x = r * cos(theta - (t / 2.0));
				*ret_y = r * sin(theta - (t / 2.0));
			} else {
				*ret_x = r * cos(theta + (t / 2.0));
				*ret_y = r * sin(theta + (t / 2.0));
			}
			break;
		case 23: /* Eyefish */
			*ret_x = (2.0 / (r + 1.0)) * x;
			*ret_y = (2.0 / (r + 1.0)) * y;
			break;
		case 24: /* Bubble */
			*ret_x = ((4.0 * x) / (pow(r,2.0) + 4.0));
			*ret_y = ((4.0 * y) / (pow(r,2.0) + 4.0));
			break;
		case 25: /* Cylinder */
			*ret_x = sin(x);
			*ret_y = y;
			break;
		case 26: /* Tangent */
			*ret_x = sin(x) / cos(y); 
			*ret_y = tan(y);
			break;
		case 27: /* Cross */
			*ret_x = x * sqrt( 1.0 / pow(x*x-y*y, 2.0) );
			*ret_y = y * sqrt( 1.0 / pow(x*x-y*y, 2.0) );
			break;
		case 28: /* Collatz */
			*ret_x = (1.0 / 4.0) * (1.0 + 4.0 * x - (1.0 + 2.0 * x) * cos(M_PI * x));

			*ret_y = (1.0 / 4.0) * (1.0 + 4.0 * y - (1.0 + 2.0 * y) * cos(M_PI * y));
			break;
		default:
			break;
	}
}

double modulus( double a,  double b)
{
	int cast;
	cast =  (int) (a / b);
	return a - ((double) cast * b);
}


/* random value 0-255 */
unsigned char channel()
{
	return (unsigned char) ((double)rand() / (((double)RAND_MAX + (double)1) / (double)256)) ;
}

/* random -1 or 1 */
 double sign()
{
	return ((((double)rand() / (((double)(RAND_MAX)+(double)1)/(double)100)) >= (double) 50)) ? -1.0 : 1.0;
}

/* return either a or b */
 double choose(double a,  double b)
{
	return (sign()==1.0) ? a : b;
}

FILE* my_file_open(char* filename)
{
	FILE* outfile;
	if((outfile = fopen(filename, "wb")) == NULL) {
		fprintf(stderr, "can't open %s\n", filename);
		exit(1);
	}
	return outfile;
}

/* initialize new fractal */
void fractal_set(int xres, int yres,  double xmin,  double xmax,  double ymin,  double ymax, fract* fractal)
{
	int y;
	(*fractal).xres = xres <= 0 ? DXRES : xres;
	(*fractal).yres = yres <= 0 ? DYRES : yres;
	(*fractal).xmin = xmax < xmin ? xmax : xmin;
	(*fractal).xmax = xmin > xmax ? xmin : xmax;
	(*fractal).ymin = ymax < ymin ? ymax : ymin;
	(*fractal).ymax = ymin > ymax ? ymin : ymax;
	(*fractal).ranx = (*fractal).xmax - (*fractal).xmin;
	(*fractal).rany = (*fractal).ymax - (*fractal).ymin;

	/* malloc new memory array for each row */
	/* then memory can be freed after being sent to JPEG */

	printf("Size of pixel structure is %ld bytes.\n", sizeof(pixel));
	printf("Attempting to allocate %ld MiB of RAM.\n", ((*fractal).yres * sizeof( pixel* ) + ((*fractal).yres * (*fractal).xres * sizeof(pixel))) / (1024 * 1024) );

	(*fractal).pixels = malloc((*fractal).yres * sizeof( pixel* ));
	if((*fractal).pixels == NULL) {
		printf("malloc() failed");
		exit(1);
	}
	for(y=0; y<(*fractal).yres; y++) {
		(*fractal).pixels[y] = malloc((*fractal).xres * sizeof(pixel));
		if((*fractal).pixels[y] == NULL) {
			printf("malloc() failed\n");
			exit(1);
		}
		memset((*fractal).pixels[y], '\0', (*fractal).xres * sizeof(pixel));
	}
}

/* take a fractal rendered with larger bit bucket and shrink it */
void reduce(fract* fractal, int sample) {
	unsigned int R, G, B, count, y, x;
	int sx, sy;
	int old_yres;
	pixel** reduction;
	old_yres = (*fractal).yres;
	(*fractal).xres = (*fractal).xres / sample;
	(*fractal).yres = (*fractal).yres / sample;
	reduction = malloc((*fractal).yres * sizeof(pixel*));
	if(reduction == NULL) {
		printf("malloc() failed\n");
		exit(1);
	}
	for(y=0; y<(*fractal).yres; y++) {
		reduction[y] = malloc((*fractal).xres * sizeof(pixel));
		if(reduction[y] == NULL) {
			printf("malloc() failed\n");
			exit(1);
		}
	}

	/* simple grid algorithm anti-aliasing */
	/* numerical average of colors in a square region */

	for(y=0; y<(*fractal).yres; y++) {
		for(x=0; x<(*fractal).xres; x++) {
			R=0; G=0; B=0; count=0;
			for(sy=0; sy<sample; sy++) {
				for(sx=0; sx<sample; sx++) {
					R += (*fractal).pixels[y*sample+sy][x*sample+sx].color.r;
					G += (*fractal).pixels[y*sample+sy][x*sample+sx].color.g;
					B += (*fractal).pixels[y*sample+sy][x*sample+sx].color.b;
				count += (*fractal).pixels[y*sample+sy][x*sample+sx].vals.count;
				}
			}
			reduction[y][x].color.r = (unsigned char) (R / (sample * sample));
			reduction[y][x].color.g = (unsigned char) (G / (sample * sample));
			reduction[y][x].color.b = (unsigned char) (B / (sample * sample));
			reduction[y][x].vals.count = count;
		}
	}

	/* replace pixel array with new, smaller array */

	for(y=0; y < old_yres; y++) {
		free((*fractal).pixels[y]);
	}
	free((*fractal).pixels);
	(*fractal).pixels = reduction;
}	

void print_usage()
{
	/* print program use */

	printf("fractal usage:\n");
	printf("fractal [-options ...]\n\n");
	printf("options include:\n");

	printf("\t-h\t\t\tprint this screen\n");
	printf("\t-R NUM\t\t\tseed randomizer with NUM\n");
	printf("\t-f NAME [%s]\tfile to write\n", PATH);
	printf("\t-J\t\t\tUse JPEG instead of TIFF\n");
	printf("\t-q QUALITY [%d]\tJPEG image quality\n", DQUALITY);
	printf("\t-S\t\t\tenable symmetry about the line y=x\n");
	printf("\t-I\t\t\tinvert colors in final image\n");
	printf("\t-x XRES [%d]\t\timage x resolution\n", DXRES);
	printf("\t-y YRES [%d]\t\timage y resolution\n", DYRES);
	printf("\t-m XMIN [%f]\t\tgraph x minimum\n", (float) DXMIN);
	printf("\t-M XMAX [%f]\t\tgraph x maximum\n", (float) DXMAX);
	printf("\t-l YMIN [%f]\t\tgraph y minimum\n", (float) DYMIN);
	printf("\t-L YMAX [%f]\t\tgraph y maximum\n", (float) DYMAX);
	printf("\t-n NUMV [%d]\t\tnumber of random vectors to use\n", DNUMV);
	printf("\t-s SAMPLES [%d]\tnumber of image samples\n", DSAMPLES);
	printf("\t-i NUM>20 [%d]\tnumber of itterations per sample\n", DITT);
	printf("\t-r 0<=NUM<=255\t\tset static RED channel value\n");
	printf("\t-g 0<=NUM<=255\t\tset static GREEN channel value\n");
	printf("\t-b 0<=NUM<=255\t\tset static BLUE channel value\n");
	printf("\t-sup NUM [%d] \t\tsuper sample NUM^2 bit buckets\n", DSUPER);
	printf("\t-o NUM\t\t\toversample points NUM times\n");
	printf("\t-u NUM\t\t\tundersample points by factor NUM\n");
	printf("\t-G NUM [%f]\t\tcorrectional gamma factor\n", (float) DGAMMA);
	printf("\t-v NUM\t\t\tuse equation by number: see below\n");

	printf("\n\nValues for v include:\n");
	printf("\t0\t\t\tLinear\n");
	printf("\t1\t\t\tSinusoidal\n");
	printf("\t2\t\t\tSpherical\n");
	printf("\t3\t\t\tSwirl\n");
	printf("\t4\t\t\tHorseshoe\n");
	printf("\t5\t\t\tPolar\n");
	printf("\t6\t\t\tHandkerchief\n");
	printf("\t7\t\t\tHeart\n");
	printf("\t8\t\t\tDisk\n");
	printf("\t9\t\t\tSpiral\n");
	printf("\t10\t\t\tHyperbolic\n");
	printf("\t11\t\t\tDiamond\n");
	printf("\t12\t\t\tEx\n");
	printf("\t13\t\t\tJulia\n");
	printf("\t14\t\t\tBent\n");
	printf("\t15\t\t\tWaves\n");
	printf("\t16\t\t\tFisheye\n");
	printf("\t17\t\t\tPopcorn\n");
	printf("\t18\t\t\tExponential\n");
	printf("\t19\t\t\tPower\n");
	printf("\t20\t\t\tCosine\n");
	printf("\t21\t\t\tRings\n");
	printf("\t22\t\t\tFan\n");
	printf("\t23\t\t\tEyefish\n");
	printf("\t24\t\t\tBubble\n");
	printf("\t25\t\t\tCylinder\n");
	printf("\t26\t\t\tTangent\n");
	printf("\t27\t\t\tCross\n");
	printf("\t28\t\t\tCollatz\n");
	fflush(stdout);
}

void write_to_jpeg(fract* fractal, int quality, char* filename, int invert)
{
	/* write out variables */
	int row,col;

	/* basic JPEG structs */
	FILE* outfile;
	JSAMPROW row_pointer[1];
	struct jpeg_compress_struct cinfo;
	struct jpeg_error_mgr jerr;

	/* Setup JPEG for compression */
	cinfo.err = jpeg_std_error(&jerr);
	jpeg_create_compress(&cinfo);

	/* Open output file for writing */
	outfile = my_file_open(filename);

	/* finalize JPEG setup */
	jpeg_stdio_dest(&cinfo, outfile);
	cinfo.image_width = (*fractal).xres;
	cinfo.image_height = (*fractal).yres;
	cinfo.input_components = 3;
	cinfo.in_color_space = JCS_RGB;
	jpeg_set_defaults(&cinfo);
	jpeg_set_quality(&cinfo, quality, FALSE);
	jpeg_start_compress(&cinfo, TRUE);

	/* begin writing to file */
	row_pointer[0] = malloc((*fractal).xres * 3 * sizeof(unsigned char));
	if(row_pointer[0] == NULL) {
		printf("malloc() failed in write_to_jpeg.\n");
		exit(1);
	}

	for(row=0;row<(*fractal).yres;row++) {
		for(col=0;col<(*fractal).xres;col++) {
			row_pointer[0][(col*3)] = invert == 1 ? ~((*fractal).pixels[row][col].color.r): 
				(*fractal).pixels[row][col].color.r;
			row_pointer[0][(col*3)+1] = invert == 1 ? ~((*fractal).pixels[row][col].color.g): 
				(*fractal).pixels[row][col].color.g;
			row_pointer[0][(col*3)+2] = invert == 1 ? ~((*fractal).pixels[row][col].color.b): 
				(*fractal).pixels[row][col].color.b;
		}
        (void) jpeg_write_scanlines(&cinfo, row_pointer, 1);
        free((*fractal).pixels[row]);
	}

	/* finalize JPEG compression and cleanup */
	free(row_pointer[0]);
	jpeg_finish_compress(&cinfo);
	fclose(outfile);
	jpeg_destroy_compress(&cinfo);
}	

void write_to_tiff(fract* fractal, char* filename, int invert)
{
	int row, col;
	TIFF *output;
	char *raster;
	
	/* Open the output image */
	if((output = TIFFOpen(filename, "w")) == NULL) {
		fprintf(stderr, "Could not open outgoing image.\n");
		exit(1);
	}

	/* malloc space for the image lines */
	raster = malloc((*fractal).xres * 3 * sizeof(char));
	if(raster == NULL) {
		printf("malloc() failed in write_to_tiff.\n");
		exit(1);
	}

	/* Write the tiff tags to the file */

  	TIFFSetField(output, TIFFTAG_IMAGEWIDTH, (*fractal).xres);
  	TIFFSetField(output, TIFFTAG_IMAGELENGTH, (*fractal).yres);
  	TIFFSetField(output, TIFFTAG_COMPRESSION, COMPRESSION_DEFLATE);
  	TIFFSetField(output, TIFFTAG_PLANARCONFIG, PLANARCONFIG_CONTIG);
  	TIFFSetField(output, TIFFTAG_PHOTOMETRIC, PHOTOMETRIC_RGB);
  	TIFFSetField(output, TIFFTAG_BITSPERSAMPLE, 8);
  	TIFFSetField(output, TIFFTAG_SAMPLESPERPIXEL, 3);

	for(row=0;row<(*fractal).yres;row++) {
		for(col=0;col<(*fractal).xres;col++) {
			raster[(col*3)] = invert == 1 ? ~((*fractal).pixels[row][col].color.r): 
				(*fractal).pixels[row][col].color.r;
			raster[(col*3)+1] = invert == 1 ? ~((*fractal).pixels[row][col].color.g): 
				(*fractal).pixels[row][col].color.g;
			raster[(col*3)+2] = invert == 1 ? ~((*fractal).pixels[row][col].color.b): 
				(*fractal).pixels[row][col].color.b;
		}
    	if(TIFFWriteScanline(output, raster,  row , (*fractal).xres * 3) != 1) {
			fprintf(stderr, "Could not write image\n");
			exit(1);
		}
    	free((*fractal).pixels[row]);
	}

	free(raster);
	/* close the file */
	TIFFClose(output);
}

void coeffs_init(coeff *coeffs, int NUMV, int R, int G, int B)
{
	int num;	
	for(num=0; num<NUMV; num++)
	{
		(*coeffs).ac[num] = sign() * ((double) rand() / ((double)(RAND_MAX)+(double) 1));
		(*coeffs).bc[num] = sign() * ((double) rand() / ((double)(RAND_MAX)+(double) 1));
		(*coeffs).cc[num] = sign() * ((double) rand() / ((double)(RAND_MAX)+(double) 1));
		(*coeffs).dc[num] = sign() * ((double) rand() / ((double)(RAND_MAX)+(double) 1));
		(*coeffs).ec[num] = sign() * ((double) rand() / ((double)(RAND_MAX)+(double) 1));
		(*coeffs).fc[num] = sign() * ((double) rand() / ((double)(RAND_MAX)+(double) 1));
		(*coeffs).colors[num].r = R != -1 ? (unsigned char) R : (unsigned char) channel(); 
		(*coeffs).colors[num].g = G != -1 ? (unsigned char) G : (unsigned char) channel(); 
		(*coeffs).colors[num].b = B != -1 ? (unsigned char) B : (unsigned char) channel(); 
	}
}

int main(int argc, char **argv)
{
	int NUMV, SAMPLES, ITTERATIONS, QUALITY, SEED, OVERSAMPLE, UNDERSAMPLE;
	int xres, yres, R, G, B, count, invert;
	int symmetry, super, use_jpeg;
	int * choice;
	double xmin, xmax, ymin, ymax;
	char * FILENAME;
	unsigned  int num, step, row, col, i;
	double x, y;
	double gamma, max;
	char last_percent = 0;
	coeff coeffs;
	fract fractal;

	/*---------------------------------------------*/ 
	/* DEFINE PRESETS  */
	NUMV = DNUMV;
	SAMPLES = DSAMPLES;
	ITTERATIONS = DITT;
	QUALITY = DQUALITY;
	OVERSAMPLE = 1;
	UNDERSAMPLE = 1;
	SEED = 1;
	R = -1;
	G = -1;
	B = -1;
	super = DSUPER;
	gamma = DGAMMA;
	xres = DXRES;
	yres = DYRES;
	ymin = DYMIN;
	xmin = DXMIN;
	xmax = DXMAX;
	ymax = DYMAX;
	choice = 0;
	invert = 0;
	symmetry = 0;
	use_jpeg = 0;
	FILENAME = PATH;
	count = 0;
	choice = malloc(count+1 * sizeof(int));
	if(choice == NULL)
		exit(1);
	choice[0] = 0;
	/*---------------------------------------------*/

	/*---------------------------------------------*/ 
	/* PARSE USER ARGUMENTS                        */
	i=1;
	while(i < argc)
	{
		
			if(!strcmp(argv[i], "-h")){
				print_usage();
				exit(0);}
			else if(!strcmp(argv[i], "-S")){
				symmetry = 1; i++;}
			else if(!strcmp(argv[i], "-I")){
				invert = 1; i++;}
			else if(!strcmp(argv[i], "-J")){
				use_jpeg = 1; i++;}
			else if(!strcmp(argv[i], "-R")){
				SEED = atoi(argv[i+1]); i+=2;}
			else if(!strcmp(argv[i], "-f")){
				FILENAME = argv[i+1]; i+=2;}
			else if(!strcmp(argv[i], "-q")){
				QUALITY = atoi(argv[i+1]); i+=2;}
			else if(!strcmp(argv[i], "-x")){
				xres = atoi(argv[i+1]); i+=2;}
			else if(!strcmp(argv[i], "-y")){
				yres = atoi(argv[i+1]); i+=2;}
			else if(!strcmp(argv[i], "-m")){
				xmin = ( double) atof(argv[i+1]); i+=2;}
			else if(!strcmp(argv[i], "-M")){
				xmax = ( double) atof(argv[i+1]); i+=2;}
			else if(!strcmp(argv[i], "-l")){
				ymin = ( double) atof(argv[i+1]); i+=2;}	
			else if(!strcmp(argv[i], "-L")){
				ymax = ( double) atof(argv[i+1]); i+=2;}
			else if(!strcmp(argv[i], "-n")){
				NUMV = atoi(argv[i+1]);	i+=2;}	
			else if(!strcmp(argv[i], "-s")){
				SAMPLES = atoi(argv[i+1]); i+=2;}		
			else if(!strcmp(argv[i], "-i")){
				ITTERATIONS = atoi(argv[i+1]) < 20 ? 1000 : atoi(argv[i+1]); i+=2;}
			else if(!strcmp(argv[i], "-r")){
				R = atoi(argv[i+1]) % 256; i+=2;}
			else if(!strcmp(argv[i], "-g")){
				G = atoi(argv[i+1]) % 256; i+=2;}
			else if(!strcmp(argv[i], "-b")){
				B = atoi(argv[i+1]) % 256; i+=2;}
			else if(!strcmp(argv[i], "-o")){
				OVERSAMPLE = atoi(argv[i+1]); i+=2;}
			else if(!strcmp(argv[i], "-u")){
				UNDERSAMPLE = atoi(argv[i+1]);
				if(UNDERSAMPLE <= 0)
					UNDERSAMPLE=1;
				i+=2;}
			else if(!strcmp(argv[i], "-G")){
				gamma = (double) atof(argv[i+1]);
				if(gamma == 0.0)
					gamma = 2.2;
				i+=2;}
			else if(!strcmp(argv[i], "-sup")){
				super = atoi(argv[i+1]);
				if(super==0)
					super++;
				i+=2;}
			else if(!strcmp(argv[i], "-v")){
				choice[count] = atoi(argv[i+1]);
				count++;
				choice = realloc(choice, (count+1) * sizeof(int));
				if(choice == NULL)
					exit(1);
				i+=2;}
			else{
				print_usage();
				exit(1);
			}
	} /* end while */

	if(count > 0) count--;

	printf("Beginning images with transformation(s): ");
	for(i=0;i<=count;i++)
		printf("%d ", choice[i]);
	fputc('\n', stdout);	
	printf("Using random seed: %d.\n", SEED);
	srand(SEED); /* seed random */
	fractal_set(xres * super, yres * super, xmin, xmax, ymin, ymax, &fractal);
	/*fractal_set(1920, 1080, -1.5625, 1.5625, -1.0, 1.0, &fractal);*/
	/*fractal_set(1024, 1024, -1.0, 1.0, -1.0, 1.0, &fractal);*/

	/* setup random coefficeints and colors */
	coeffs.ac = malloc(NUMV * sizeof(double));
	if(coeffs.ac == NULL)
		exit(1);
	coeffs.bc = malloc(NUMV * sizeof(double));
	if(coeffs.bc == NULL)
		exit(1);
	coeffs.cc = malloc(NUMV * sizeof(double));
	if(coeffs.cc == NULL)
		exit(1);
	coeffs.dc = malloc(NUMV * sizeof(double));
	if(coeffs.dc == NULL)
		exit(1);
	coeffs.ec = malloc(NUMV * sizeof(double));
	if(coeffs.ec == NULL)
		exit(1);
	coeffs.fc = malloc(NUMV * sizeof(double));
	if(coeffs.fc == NULL)
		exit(1);
	coeffs.colors = malloc(NUMV * sizeof(trio));
	if(coeffs.colors == NULL)
		exit(1);

	coeffs_init(&coeffs, NUMV, R, G, B);

	printf("Image will require %d iterations, at %d samples per iteration.\n", SAMPLES, ITTERATIONS);

	printf("Progress:   ");
	
	for(num=0; num<SAMPLES; num++)
	{
		if( (char) ((double) num / (double) SAMPLES * 100.0) > last_percent )
		{
			last_percent = (char) ((double) num / (double) SAMPLES * 100.0);
			if(last_percent <= 10) {
				printf("\b\b%d%%", (int) last_percent);
			} else {
				printf("\b\b\b%d%%", (int) last_percent);
			}
			fflush(stdout);
		}

		x = ((double) rand()) * choose(fractal.xmin, fractal.xmax);
		y = ((double) rand()) * choose(fractal.ymin, fractal.ymax);

		for(step=0; step<ITTERATIONS; step++)
		{
			i = (unsigned int) (rand() % NUMV);
			V(choice[rand() % (count+1)], 
				coeffs.ac[i] * x + coeffs.bc[i] * y + coeffs.cc[i], 
				coeffs.dc[i] * x + coeffs.ec[i] * y + coeffs.fc[i],
				(int)i,
				&x,
				&y,
				&coeffs);

			if(step > 20) {
				unsigned  int x1, y1;
				unsigned short red, green, blue;
				int valid;

				valid = 0;
				y1=0;
				x1=0;

				if(x >= fractal.xmin && x <= fractal.xmax && y >= fractal.ymin && y <= fractal.ymax) {
					valid = 1;
					x1 = fractal.xres - (unsigned int) (((fractal.xmax - x) / fractal.ranx) * (double) fractal.xres);
					y1 = fractal.yres - (unsigned int) (((fractal.ymax - y) / fractal.rany) * (double) fractal.yres);
				}

				if(valid && x1 >= 0 && x1 < fractal.xres && y1 >= 0 && y1 < fractal.yres) {
					fractal.pixels[y1][x1].vals.count += OVERSAMPLE;
					red = (unsigned short) fractal.pixels[y1][x1].color.r;
					red += (unsigned short) coeffs.colors[i].r;
					fractal.pixels[y1][x1].color.r = (unsigned char) (red / 2);
					blue = (unsigned short) fractal.pixels[y1][x1].color.b;
					blue += (unsigned short) coeffs.colors[i].b;
					fractal.pixels[y1][x1].color.b = (unsigned char) (blue / 2);
					green = (unsigned short) fractal.pixels[y1][x1].color.g;
					green += (unsigned short) coeffs.colors[i].g;
					fractal.pixels[y1][x1].color.g = (unsigned char) (green / 2);

					if(symmetry) { /* only if we have symmetry enabled, plots the mirror */
						if(!x1)	x1++;
						if(!y1)	y1++;
						x1 = fractal.xres - x1; 
						y1 = fractal.yres - y1;
						if(x1 >= 0 && x1 < fractal.xres && y1 >= 0 && y1 < fractal.yres) {
							fractal.pixels[y1][x1].vals.count += OVERSAMPLE;
							fractal.pixels[y1][x1].color.r = (unsigned char) (red / 2);
							fractal.pixels[y1][x1].color.b = (unsigned char) (blue / 2);
							fractal.pixels[y1][x1].color.g = (unsigned char) (green / 2);
						}

					}
				}
			}
		}
	}
	
	if(super > 1)
		reduce(&fractal, super);

	max = 0.0;

	for(row=0; row<fractal.yres; row++) {
		for(col=0; col<fractal.xres; col++){
			if(fractal.pixels[row][col].vals.count != 0) {
				/*if(fractal.pixels[row][col].vals.count > max)
					max = fractal.pixels[row][col].vals.count;*/
				fractal.pixels[row][col].vals.normal = log((double) (fractal.pixels[row][col].vals.count / UNDERSAMPLE) );
				if(fractal.pixels[row][col].vals.normal > max)
					max = fractal.pixels[row][col].vals.normal;
				}
			}
		}

	for(row=0; row<fractal.yres; row++) {
		for(col=0; col<fractal.xres; col++) {
			fractal.pixels[row][col].vals.normal  = fractal.pixels[row][col].vals.normal / (float) max;
			fractal.pixels[row][col].color.r = (unsigned char) ((float)(fractal.pixels[row][col].color.r)) * pow(fractal.pixels[row][col].vals.normal, (1.0 / gamma));
			fractal.pixels[row][col].color.g = (unsigned char) ((float)(fractal.pixels[row][col].color.g)) * pow(fractal.pixels[row][col].vals.normal, (1.0 / gamma));
			fractal.pixels[row][col].color.b = (unsigned char) ((float)(fractal.pixels[row][col].color.b)) * pow(fractal.pixels[row][col].vals.normal, (1.0 / gamma));

		}
	}

	printf("\b\b\b100%%\n");
	printf("Writing out to file...\n");

	if(use_jpeg)
		write_to_jpeg(&fractal, QUALITY, FILENAME, invert);
	else
		write_to_tiff(&fractal, FILENAME, invert);
	
	printf("Cleaning up...\n");
	
	/* clean up */
	free(fractal.pixels);
	free(choice);
	free(coeffs.ac);
	free(coeffs.bc);
	free(coeffs.cc);
	free(coeffs.dc);
	free(coeffs.ec);
	free(coeffs.fc);
	free(coeffs.colors);
	printf("Finished.\n");
	return 0;
}
