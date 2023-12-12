/* Texture functions for cs580 GzLib	*/
#include    "stdafx.h" 
#include	"stdio.h"
#include	"Gz.h"

GzColor	*image=NULL;
int xs, ys;
int reset = 1;

void getColor(int i, int j, GzColor* image, int xs, GzColor c) {
    c[RED] = image[(i + j * xs)][RED];
    c[GREEN] = image[(i + j * xs)][GREEN];
    c[BLUE] = image[(i + j * xs)][BLUE];
}

/* Image texture function */
int tex_fun(float u, float v, GzColor color)
{
  unsigned char		pixel[3];
  unsigned char     dummy;
  char  		foo[8];
  int   		i, j;
  FILE			*fd;

  if (reset) {          /* open and load texture file */
    fd = fopen ("texture", "rb");
    if (fd == NULL) {
      fprintf (stderr, "texture file not found\n");
      exit(-1);
    }
    fscanf (fd, "%s %d %d %c", foo, &xs, &ys, &dummy);
    image = (GzColor*)malloc(sizeof(GzColor)*(xs+1)*(ys+1));
    if (image == NULL) {
      fprintf (stderr, "malloc for texture image failed\n");
      exit(-1);
    }

    for (i = 0; i < xs*ys; i++) {	/* create array of GzColor values representing the texture image*/
      fread(pixel, sizeof(pixel), 1, fd);
      image[i][RED] = (float)((int)pixel[RED]) * (1.0 / 255.0);
      image[i][GREEN] = (float)((int)pixel[GREEN]) * (1.0 / 255.0);
      image[i][BLUE] = (float)((int)pixel[BLUE]) * (1.0 / 255.0);
      }

    reset = 0;          /* init is done */
	fclose(fd);
  }

  // scale u, v to texture size range
  u *= (xs - 1);
  v *= (ys - 1);

  // clamp to bounds test u, v
  if (u < 1) u = 1;
  if (u > xs - 2) u = xs - 2;
  if (v < 1) v = 1;
  if (v > ys - 2) v = ys - 2;

  // determine texture cell corner values
  int I = u;
  int J = v;

  GzColor A;
  getColor(I, J, image, xs, A);

  GzColor B;
  getColor(I + 1, J, image, xs, B);

  GzColor C;
  getColor(I + 1, J + 1, image, xs, C);

  GzColor D;
  getColor(I, J + 1, image, xs, D);

  float s = u - I;
  float t = v - J;

  // perform bilinear interpolation
  color[RED] = s * t * C[RED] + (1 - s) * t * D[RED] + s * (1 - t) * B[RED] + (1 - s) * (1 - t) * A[RED];
  color[GREEN] = s * t * C[GREEN] + (1 - s) * t * D[GREEN] + s * (1 - t) * B[GREEN] + (1 - s) * (1 - t) * A[GREEN];
  color[BLUE] = s * t * C[BLUE] + (1 - s) * t * D[BLUE] + s * (1 - t) * B[BLUE] + (1 - s) * (1 - t) * A[BLUE];
  
  return GZ_SUCCESS;
}

/* Procedural texture function */
int ptex_fun(float u, float v, GzColor color)
{
    int N = 100;

    // determine which interval u and v fall into
    int intervalU = (int)(u * N) % 10;  // 0 for even, 1 for odd
    int intervalV = (int)(v * N) % 10;  // 0 for even, 1 for odd

    // assign color
    if ((intervalU == 0 || intervalV == 0)) {
        // black color for same intervals (even/even or odd/odd)
        color[0] = 0.0f;
        color[1] = 0.0f;
        color[2] = 0.0f;
    }
    else {
        // white color for different intervals (even/odd or odd/even)
        color[0] = 1.0f;
        color[1] = 1.0f;
        color[2] = 1.0f;
    }

	return GZ_SUCCESS;
}

/* Free texture memory */
int GzFreeTexture()
{
	if(image!=NULL)
		free(image);
	return GZ_SUCCESS;
}

