
/* M1: no changes yet */

/*

modified 17 June 2012 CFG
        - eliminate r8 option (never used).
        - eliminate read-in of color map (map.ppm no longer required)
        - insert analytic function for john.pal color map
        - removed failimage capability (probably unwise, but eliminates
                global arrays for failimage) 

*/

#include "decs.h"
#include <ctype.h>

/* Local Mnemonics for image types */
#define IRHO  (0)
#define ITT   (1)
#define IBSQ  (2)
#define IGAM  (3)
#define GJSQ  (4)

/* number of types of images to make */
#define NIMG   5        // = # mnemonics

double fimage[NIMG][N1 * N2];

/* Local functions */
void image_ppm(double *f, char *fname);


/******************************************************************************/
/******************************************************************************
  image_all(): 
  -----------
   -- Main driver for generating "image" or snapshots of any desired quantity;
   -- This is esponsible for calculating the outputted quantities and setting 
       the names of the images
   -- This is the only routine in this file that a user should modify, all 
      other routines merely control how the images are created. 
   -- The image generating routines use a linear scale, so be sure to 
      take the log of any quantity you want to see in log color scale.

******************************************************************************/
void image_all(int image_count)
{

	int i, j, k, i_img;
	static char ifnam[3 * NIMG + 1][100];
	double gamma;
        double gJsq ;
	struct of_geom *geom;
	static const double fimage_logmin = 1.e-15;

  /************************************************************************
    Set the names of the image files to be generated now : 
  ************************************************************************/
	i_img = 0;
	sprintf(ifnam[i_img++], "images/im_rho_%04d.ppm", image_count);
	sprintf(ifnam[i_img++], "images/im_T_%04d.ppm", image_count);
	sprintf(ifnam[i_img++], "images/im_bsq_%04d.ppm", image_count);
	sprintf(ifnam[i_img++], "images/im_gam_%04d.ppm", image_count);
	sprintf(ifnam[i_img++], "images/im_gJ_%04d.ppm", image_count);

	sprintf(ifnam[i_img++], "images/im_lrho_%04d.ppm", image_count);
	sprintf(ifnam[i_img++], "images/im_lT_%04d.ppm", image_count);
	sprintf(ifnam[i_img++], "images/im_lbsq_%04d.ppm", image_count);
	sprintf(ifnam[i_img++], "images/im_lgam_%04d.ppm", image_count);
	sprintf(ifnam[i_img++], "images/im_lgJ_%04d.ppm", image_count);

  /************************************************************************
    Calculate the functions to be imaged : 
        -- the log versions overwrite the non-log versions;
        -- calculate only the non-log version here
  ************************************************************************/
	double Pressure_rho0_u(double rho, double u);	/* needed for calculating T */
        struct of_state q ;
        double jcov[NDIM],Jdu,Jsq,gJsum ;

        current_calc() ;
	k = 0;
        gJsum = 0. ;
	IMAGELOOP {     // 2D array in i,j written into 1D array in k
		geom = get_geometry(i, j, CENT) ;
                
		if (mhd_gamma_calc(p[i][j], geom, &gamma)) {
			gamma = 1.;
		}

		fimage[IRHO][k] = p[i][j][RHO];
		fimage[ITT][k] = Pressure_rho0_u(p[i][j][RHO], p[i][j][UU]) / p[i][j][RHO];
		fimage[IBSQ][k] = bsq_calc(p[i][j], geom);
		fimage[IGAM][k] = gamma;

                get_state(p[i][j], geom, &q);
                lower(Jcon[i][j], geom, jcov) ;
                Jsq = jcov[0]*Jcon[i][j][0] + jcov[1]*Jcon[i][j][1] + 
                        jcov[2]*Jcon[i][j][2] + jcov[3]*Jcon[i][j][3] ;
                Jdu = jcov[0]*q.ucon[0] + jcov[1]*q.ucon[1] + 
                        jcov[2]*q.ucon[2] + jcov[3]*q.ucon[3] ;
                gJsq = geom->g * (Jsq + Jdu*Jdu) ;
                gJsum += gJsq ;
                fimage[GJSQ][k] = gJsq ;
		k++;
	}
        fprintf(stderr,"\n t,gJ: %10.5g %10.5g\n",t,gJsum) ;


  /************************************************************************
    Output non-log versions:
  ************************************************************************/
	for (i_img = 0; i_img < NIMG; i_img++) {
		image_ppm(fimage[i_img], ifnam[i_img]);
	}

  /************************************************************************
    Make log version of the image functions 
  ************************************************************************/
	for (i = 0; i < NIMG * N1 * N2; i++) {
		fimage[0][i] = log(fabs(fimage[0][i]) + fimage_logmin);
	}
	for (i_img = 0; i_img < NIMG; i_img++) {
		image_ppm(fimage[i_img], ifnam[i_img + NIMG]);
	}

	return;
}

/******************************************************************************
  image_ppm(): 
  -----------
     -- generates a color mapped "image" or pixelated output file following 
        the "raw" PPM  format ("man ppm" for more details). 

     -- color map determined by get_color_map()

    CFG 14 Sept 07: modified to go from almost-min to almost-max

******************************************************************************/
void image_ppm(double *f, char *fname)
{
	int i ;
	double max, min ;
	static double q[N1 * N2];
        int red, green, blue ;
	int compare_doubles(const void *a, const void *b);
        void john_pal(double data, double min, double max, int *pRed, int *pGreen, int *pBlue) ;
	FILE *fp;

	if ((fp = fopen(fname, "w")) == NULL) {
		fflush(stderr);
		fprintf(stderr, "image(): Cannot open %s !! \n", fname);
		fflush(stderr);
		return;
	}

	/*  mapping is in 255 steps lmax and lmin */
	for (i = 0; i < N1 * N2; i++) q[i] = f[i];
                // this bit sorts the values and lops off the top and bottom
                // end of the distribution to create max and min.
		// qsort is a standard unix utility now
	qsort(q, N1 * N2, sizeof(double), compare_doubles);
	min = q[N1 * N2 / 128];
	max = q[N1 * N2 - N2 * N2 / 128];

	/* Header information: */
	fprintf(fp, "P6\n#  min=%g  , max=%g \n%d %d\n%d\n", min, max, N1, N2, 255);
	fflush(fp);

	for (i = 0; i < N1 * N2; i++) {
                john_pal(f[i], min, max ,&red,&green,&blue) ;
		fputc((char) red, fp);
		fputc((char) green, fp);
		fputc((char) blue, fp);
	}

	fclose(fp);

	return;
}

/* this is needed for qsort, which is used to set colormap min and max */
int compare_doubles(const void *a, const void *b)
{
	const double *da = (const double *) a;
	const double *db = (const double *) b;

	return (*da > *db) - (*da < *db);
}

/* 
color palette, based on john.pal

        input: integer 0-255 
        output: red, green, blue integers 

        author: Bryan M. Johnson
*/

void john_pal(double data, double min, double max, int *pRed, int *pGreen, int *pBlue)
{
  double a, b, c, d, e, f;
  double x, y;
  double max_min = max - min ;

  if(max_min > 0.0) { // trust no one
    x = (data - min)/(max_min) ;
    
    /* ========== Red ============ */
    a = 4.0*x - 1.52549019607844;
    b = 4.52941176470589 - 4.0*x;
    y = a < b ? a : b;
    *pRed = (int)(255.0*y);
    *pRed = *pRed >   0 ? *pRed :   0;
    *pRed = *pRed < 255 ? *pRed : 255;

    /* ========== Green ========== */
    a = 4.0*x - 0.521568627450979;
    b = 2.52549019607844 - 4.0*x;
    c = a < b ? a : b;
    d = 4.0*x - 1.53725490196073;
    e = 3.52941176470581 - 4.0*x;
    f = d < e ? d : e;
    y = c > f ? c : f;
    *pGreen = (int)(255.0*y);
    *pGreen = *pGreen >   0 ? *pGreen :   0;
    *pGreen = *pGreen < 255 ? *pGreen : 255;

    /* ========== Blue =========== */
    a = 4.0*x + 0.498039215686276;
    b = 2.50980392156862 - 4.0*x;
    y = a < b ? a : b;
    *pBlue = (int)(255.0*y);
    *pBlue = *pBlue >   0 ? *pBlue :   0;
    *pBlue = *pBlue < 255 ? *pBlue : 255;

  }
  else {
    *pRed = *pGreen = *pBlue = (data > max ? 255: 0) ;
  }

  return;
}

#undef IRHO
#undef ITT
#undef IBSQ
#undef IGAM
#undef GJSQ
