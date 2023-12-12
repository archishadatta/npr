/* CS580 Homework 3 */

#include	"stdafx.h"
#include	"stdio.h"
#include	"math.h"
#include	"Gz.h"
#include	"rend.h"
#include <cmath>
#include <iostream>
#include <random>


#define PI (float) 3.14159265358979323846


// Helper Functions

// Matrix operation helper functions
void GzRender::normalize(GzCoord vector) {
	float mag = sqrt(vector[X] * vector[X] + vector[Y] * vector[Y] + vector[Z] * vector[Z]);
	vector[X] /= mag;
	vector[Y] /= mag;
	vector[Z] /= mag;
}

void GzRender::matrixMultiplyVector(const GzMatrix mat, const float input[4], float output[4])
{
	for (int i = 0; i < 4; i++) {
		output[i] = 0;
		for (int j = 0; j < 4; j++) {
			output[i] += mat[i][j] * input[j];
		}
	}
}

float* GzRender::edgeVector(const float v1[], const float v2[]) {
	float result[] = {
		v1[0] - v2[0],
		v1[1] - v2[1],
		v1[2] - v2[2]
	};
	return result;
}

float* GzRender::crossProduct(const float v1[], const float v2[]) {
	float result[] = {
		v1[1] * v2[2] - v1[2] * v2[1],
		v1[2] * v2[0] - v1[0] * v2[2],
		v1[0] * v2[1] - v1[1] * v2[0]
	};
	return result;
}

float GzRender::dotProduct(const float v1[], const float v2[]) {
	double result = 0.0;

	for (int i = 0; i < 3; ++i) {
		result += v1[i] * v2[i];
	}

	return result;
}

void GzRender::identity(GzMatrix mat) {
	for (int i = 0; i < 4; i++) {
		for (int j = 0; j < 4; j++) {
			if (i == j)
				mat[i][j] = 1.0f;
			else
				mat[i][j] = 0.0f;
		}
	}
}

// Interpolation and color helper functions
float GzRender::interpolate(const float vec1[], const float vec2[], const float vec3[], float n1, float n2, float n3, int x, int y) {

	float v1[] = { vec1[0], vec1[1], n1 };
	float v2[] = { vec2[0], vec2[1], n2 };
	float v3[] = { vec3[0], vec3[1], n3 };

	float* edge1 = edgeVector(v2, v1);
	float edge1Arr[] = { edge1[0], edge1[1], edge1[2] };
	float* edge2 = edgeVector(v3, v1);
	float edge2Arr[] = { edge2[0], edge2[1], edge2[2] };
	float* cross = crossProduct(edge1Arr, edge2Arr);
	float A = cross[0];
	float B = cross[1];
	float C = cross[2];
	float D = -(A * v1[X] + B * v1[Y] + C * v1[Z]);

	// Plane equation: Ax+By+Cn+D=0

	float interpolated = (-A * x - B * y - D) / C;

	return interpolated;
}

void GzRender::getColor(float normal[], float color[], float Kd[], float Ks[], float Ka[]) {

	float N[] = { normal[X], normal[Y], normal[Z] };

	// E vector is (0, 0, -1)
	float E[] = { 0.0, 0.0, -1.0 };

	// For every light
	for (int i = 0; i < numlights; i++) {
		GzLight light = lights[i];
		// L (direction of light)
		float L[] = { light.direction[X], light.direction[Y], light.direction[Z] };
		// Ie (color of light)
		float Ie[] = { light.color[0],  light.color[1], light.color[2] };

		// Check sign of N * L, N * E
		float NdL = dotProduct(N, L);
		float NdE = dotProduct(N, E);
		if (NdL < 0 && NdE < 0) {
			// flip normal
			N[X] = -N[X];
			N[Y] = -N[Y];
			N[Z] = -N[Z];
		}
		else if (NdL * NdE <= 0) {
			// skip light
			continue;
		}

		// Compute R vector
		float R[] = {
			2 * NdL * N[0] - L[0],
			2 * NdL * N[1] - L[1],
			2 * NdL * N[2] - L[2],
		};
		normalize(R);
		float RdE = dotProduct(R, E);
		float min = 0.0;
		float max = 1.0;
		RdE = clamp(RdE, min, max);

		// Add diffuse and specular term
		color[0] += Ks[0] * Ie[0] * std::pow(RdE, spec) + Kd[0] * Ie[0] * NdL;
		color[1] += Ks[1] * Ie[1] * std::pow(RdE, spec) + Kd[1] * Ie[1] * NdL;
		color[2] += Ks[2] * Ie[2] * std::pow(RdE, spec) + Kd[2] * Ie[2] * NdL;

	}

	// Add ambient term
	float Ia[] = { ambientlight.color[0],  ambientlight.color[1], ambientlight.color[2] };
	color[0] += Ka[0] * Ia[0];
	color[1] += Ka[1] * Ia[1];
	color[2] += Ka[2] * Ia[2];

	// Clamp each term to 0-1
	float min = 0.0;
	float max = 1.0;
	color[0] = clamp(abs(color[0]), min, max);
	color[1] = clamp(abs(color[1]), min, max);
	color[2] = clamp(abs(color[2]), min, max);
}

float perspectiveTransform(float P, float z) {
	float vPrime = z / (MAXINT - z);
	return P / (vPrime + 1);
}

float shadingTransform(float Ps, float z) {
	float vPrime = z / (MAXINT - z);
	return Ps * (vPrime + 1);
}

// Xform helper functions
void GzRender::setCamXiw(GzCamera* camera) {
	GzCoord camX, camY, camZ;
	float camPosX = camera->position[X];
	float camPosY = camera->position[Y];
	float camPosZ = camera->position[Z];

	// camera z a-xis
	camZ[X] = camera->lookat[X] - camPosX;
	camZ[Y] = camera->lookat[Y] - camPosY;
	camZ[Z] = camera->lookat[Z] - camPosZ;
	normalize(camZ);

	// camera y-axis
	float upZ = camera->worldup[X] * camZ[X] + camera->worldup[Y] * camZ[Y] +
		camera->worldup[Z] * camZ[Z];
	camY[X] = camera->worldup[X] - upZ * camZ[X];
	camY[Y] = camera->worldup[Y] - upZ * camZ[Y];
	camY[Z] = camera->worldup[Z] - upZ * camZ[Z];
	normalize(camY);

	// camera x-axis
	camX[X] = camY[Y] * camZ[Z] - camY[Z] * camZ[Y];
	camX[Y] = camY[Z] * camZ[X] - camY[X] * camZ[Z];
	camX[Z] = camY[X] * camZ[Y] - camY[Y] * camZ[X];
	normalize(camX);

	// set Xiw matrix

	camera->Xiw[0][0] = camX[X];
	camera->Xiw[0][1] = camX[Y];
	camera->Xiw[0][2] = camX[Z];
	camera->Xiw[0][3] = -(camX[X] * camPosX + camX[Y] * camPosY + camX[Z] * camPosZ);

	camera->Xiw[1][0] = camY[X];
	camera->Xiw[1][1] = camY[Y];
	camera->Xiw[1][2] = camY[Z];
	camera->Xiw[1][3] = -(camY[X] * camPosX + camY[Y] * camPosY + camY[Z] * camPosZ);

	camera->Xiw[2][0] = camZ[X];
	camera->Xiw[2][1] = camZ[Y];
	camera->Xiw[2][2] = camZ[Z];
	camera->Xiw[2][3] = -(camZ[X] * camPosX + camZ[Y] * camPosY + camZ[Z] * camPosZ);
}

void GzRender::setCamXpi(GzCamera* camera) {
	float dinverse = tanf((camera->FOV * PI / 180.0) / 2.0);
	camera->Xpi[2][2] = dinverse;
	camera->Xpi[3][2] = dinverse;
}

// -------------- HOMEWORK FUNCTIONS ---------------

int GzRender::GzRotXMat(float degree, GzMatrix mat)
{
	/* HW 3.1
	// Create rotate matrix : rotate along x axis
	// Pass back the matrix using mat value
	*/

	float radians = degree * PI / 180.0f;
	//identity(mat);

	mat[1][1] = cosf(radians);
	mat[1][2] = -sinf(radians);
	mat[2][1] = sinf(radians);
	mat[2][2] = cosf(radians);

	return GZ_SUCCESS;
}

int GzRender::GzRotYMat(float degree, GzMatrix mat)
{
	/* HW 3.2
	// Create rotate matrix : rotate along y axis
	// Pass back the matrix using mat value
	*/
	float radians = degree * PI / 180.0f;
	//identity(mat);

	mat[0][0] = cosf(radians);
	mat[0][2] = sinf(radians);
	mat[2][0] = -sinf(radians);
	mat[2][2] = cosf(radians);

	return GZ_SUCCESS;
}

int GzRender::GzRotZMat(float degree, GzMatrix mat)
{
	/* HW 3.3
	// Create rotate matrix : rotate along z axis
	// Pass back the matrix using mat value
	*/
	float radians = degree * PI / 180.0f;
	//identity(mat);

	mat[0][0] = cosf(radians);
	mat[0][1] = -sinf(radians);
	mat[1][0] = sinf(radians);
	mat[1][1] = cosf(radians);

	return GZ_SUCCESS;
}

int GzRender::GzTrxMat(GzCoord translate, GzMatrix mat)
{
	/* HW 3.4
	// Create translation matrix
	// Pass back the matrix using mat value
	*/
	//identity(mat);

	mat[0][3] = translate[X];
	mat[1][3] = translate[Y];
	mat[2][3] = translate[Z];

	return GZ_SUCCESS;
}

int GzRender::GzScaleMat(GzCoord scale, GzMatrix mat)
{
	/* HW 3.5
	// Create scaling matrix
	// Pass back the matrix using mat value
	*/
	//identity(mat);

	mat[0][0] = scale[X];
	mat[1][1] = scale[Y];
	mat[2][2] = scale[Z];
	return GZ_SUCCESS;
}

GzRender::GzRender(int xRes, int yRes)
{
	/* HW1.1 create a framebuffer for MS Windows display:
	 -- set display resolution
	 -- allocate memory for framebuffer : 3 bytes(b, g, r) x width x height
	 -- allocate memory for pixel buffer
	 */

	 // set display resolution
	if (xRes > MAXXRES || xRes < 1) {
		xres = MAXXRES;
	}
	else {
		xres = xRes;
	}

	if (yRes > MAXYRES || yRes < 1) {
		yres = MAXYRES;
	}
	else {
		yres = yRes;
	}

	// allocate memory for frame buffer
	framebuffer = new char[3 * xRes * yRes];

	// allocate memory for pixel buffer
	pixelbuffer = new GzPixel[xRes * yRes];

	/* HW 3.6
	- setup Xsp and anything only done once
	- init default camera
	*/

	// setup xsp
	identity(Xsp);
	Xsp[0][0] = xres / 2;
	Xsp[0][3] = xres / 2;
	Xsp[1][1] = -yres / 2;
	Xsp[1][3] = yres / 2;
	Xsp[2][2] = MAXINT;

	// init default camera
	m_camera.FOV = DEFAULT_FOV;

	m_camera.lookat[0] = 0;
	m_camera.lookat[1] = 0;
	m_camera.lookat[2] = 0;

	m_camera.position[0] = DEFAULT_IM_X;
	m_camera.position[1] = DEFAULT_IM_Y;
	m_camera.position[2] = DEFAULT_IM_Z;

	m_camera.worldup[0] = 0;
	m_camera.worldup[1] = 1;
	m_camera.worldup[2] = 0;

	// init default color
	GzColor ka = DEFAULT_AMBIENT;
	GzColor kd = DEFAULT_DIFFUSE;
	GzColor ks = DEFAULT_SPECULAR;

	Ka[RED] = ka[RED];
	Ka[GREEN] = ka[GREEN];
	Ka[BLUE] = ka[BLUE];

	Kd[RED] = kd[RED];
	Kd[GREEN] = kd[GREEN];
	Kd[BLUE] = kd[BLUE];

	Ks[RED] = ks[RED];
	Ks[GREEN] = ks[GREEN];
	Ks[BLUE] = ks[BLUE];

	spec = DEFAULT_SPEC;

	// init Xiw
	identity(m_camera.Xiw);

	// init Xpi
	identity(m_camera.Xpi);

	// init remaining variables
	matlevel = 0;
	numlights = 0;
	tex_fun = NULL;

}

GzRender::~GzRender()
{
	/* HW1.2 clean up, free buffer memory */
	delete[] pixelbuffer;
	delete[] framebuffer;
}

int GzRender::GzDefault()
{
	/* HW1.3 set pixel buffer to some default values - start a new frame */
	for (int i = 0; i < xres * yres; i++) {
		pixelbuffer[i].red = 4095;
		pixelbuffer[i].green = 4095;
		pixelbuffer[i].blue = 4095;
		pixelbuffer[i].alpha = 1;
		pixelbuffer[i].visible = false;
		pixelbuffer[i].z = MAXINT;
	}
	return GZ_SUCCESS;
}

int GzRender::GzBeginRender()
{
	/* HW 3.7
	- setup for start of each frame - init frame buffer color,alpha,z
	- compute Xiw and projection xform Xpi from camera definition
	- init Ximage - put Xsp at base of stack, push on Xpi and Xiw
	- now stack contains Xsw and app can push model Xforms when needed
	*/

	// init frame buffer color, alpha, z
	for (int i = 0; i < xres * yres; i++) {
		pixelbuffer[i].red = 3700;
		pixelbuffer[i].green = 830;
		pixelbuffer[i].blue = 3700;
		pixelbuffer[i].alpha = 1;
		pixelbuffer[i].z = MAXINT;
	}

	// compute Xiw
	setCamXiw(&m_camera);

	// compute Xpi
	setCamXpi(&m_camera);

	// push matrices
	GzPushMatrix(Xsp, true);
	GzPushMatrix(m_camera.Xpi, true);
	GzPushMatrix(m_camera.Xiw);

	return GZ_SUCCESS;
}

int GzRender::GzPutCamera(GzCamera camera)
{
	/* HW 3.8
	/*- overwrite renderer camera structure with new camera definition
	*/
	// compute Xiw
	identity(camera.Xiw);
	setCamXiw(&camera);

	// compute Xpi
	identity(camera.Xpi);
	setCamXpi(&camera);

	m_camera = camera;

	return GZ_SUCCESS;
}

int GzRender::GzPushMatrix(GzMatrix	matrix, bool normIdentity)
{
	/* HW 3.9
	- push a matrix onto the Ximage stack
	- check for stack overflow
	*/

	if (matlevel >= MATLEVELS) {
		// handle stack overflow
		return GZ_FAILURE;
	}

	for (int i = 0; i < 4; i++) {
		for (int j = 0; j < 4; j++) {
			Ximage[matlevel][i][j] = 0;
			Xnorm[matlevel][i][j] = 0;
		}
	}


	if (matlevel == 0)
		for (int i = 0; i < 4; i++)
			for (int j = 0; j < 4; j++)
				Ximage[matlevel][i][j] = matrix[i][j];
	else
		for (int i = 0; i < 4; i++)
			for (int j = 0; j < 4; j++)
				for (int m = 0; m < 4; m++)
					Ximage[matlevel][i][j] += Ximage[matlevel - 1][i][m] * matrix[m][j];

	// Push onto norm stack
	if (normIdentity) {
		identity(matrix);
	}
	else {
		// Pre-process matrix
		// remove translations
		matrix[0][3] = 0;
		matrix[1][3] = 0;
		matrix[2][3] = 0;

		// normalize
		float mag = sqrt(matrix[0][X] * matrix[0][X] + matrix[0][Y] * matrix[0][Y] + matrix[0][Z] * matrix[0][Z]);
		for (int i = 0; i < 3; i++)
			for (int j = 0; j < 3; j++)
				matrix[i][j] /= mag;
	}

	if (matlevel == 0)
		for (int i = 0; i < 4; i++)
			for (int j = 0; j < 4; j++)
				Xnorm[matlevel][i][j] = matrix[i][j];
	else
		for (int i = 0; i < 4; i++)
			for (int j = 0; j < 4; j++)
				for (int m = 0; m < 4; m++)
					Xnorm[matlevel][i][j] += Xnorm[matlevel - 1][i][m] * matrix[m][j];

	matlevel++;
	return GZ_SUCCESS;
}

int GzRender::GzPopMatrix()
{
	/* HW 3.10
	- pop a matrix off the Ximage stack
	- check for stack underflow
	*/

	if (matlevel <= 0) {
		// handle stack underflow
		return GZ_FAILURE;
	}

	matlevel--;
	return GZ_SUCCESS;
}

int GzRender::GzPut(int i, int j, GzIntensity r, GzIntensity g, GzIntensity b, GzIntensity a, GzDepth z)
{
	/* HW1.4 write pixel values into the buffer */
	if (i < 0 || i >= xres || j < 0 || j >= yres) {
		return GZ_FAILURE;
	}

	// clamp values
	r = clamp<GzIntensity>(r, 0, 4095);
	g = clamp<GzIntensity>(g, 0, 4095);
	b = clamp<GzIntensity>(b, 0, 4095);
	a = clamp<GzIntensity>(a, 0, 4095);
	z = clamp<GzDepth>(z, 0, INT_MAX);

	// update the pixel values
	int ind = ARRAY(i, j);
	pixelbuffer[ind].red = r;
	pixelbuffer[ind].green = g;
	pixelbuffer[ind].blue = b;
	pixelbuffer[ind].alpha = a;
	pixelbuffer[ind].z = z;

	return GZ_SUCCESS;
}

int GzRender::GzGet(int i, int j, GzIntensity* r, GzIntensity* g, GzIntensity* b, GzIntensity* a, GzDepth* z)
{
	/* HW1.5 retrieve a pixel information from the pixel buffer */

	if (i < 0 || i >= xres || j < 0 || j >= yres) {
		return GZ_FAILURE;
	}

	int ind = ARRAY(i, j);
	*r = pixelbuffer[ind].red;
	*g = pixelbuffer[ind].green;
	*b = pixelbuffer[ind].blue;
	*a = pixelbuffer[ind].alpha;
	*z = pixelbuffer[ind].z;

	return GZ_SUCCESS;
}

int GzRender::GzFlushDisplay2File(FILE* outfile)
{
	/* HW1.6 write image to ppm file -- "P6 %d %d 255\r" */
		// header
	fprintf(outfile, "P6 %d %d 255\r", xres, yres);

	// write each pixel
	for (int j = 0; j < yres; ++j) {
		for (int i = 0; i < xres; ++i) {
			int ind = ARRAY(i, j);

			unsigned char blue = (unsigned char)(pixelbuffer[ind].blue >> 4);
			unsigned char green = (unsigned char)(pixelbuffer[ind].green >> 4);
			unsigned char red = (unsigned char)(pixelbuffer[ind].red >> 4);

			fputc(red, outfile);
			fputc(green, outfile);
			fputc(blue, outfile);
		}
	}

	return GZ_SUCCESS;
}

int GzRender::GzFlushDisplay2FrameBuffer()
{
	/* HW1.7 write pixels to framebuffer:
		- put the pixels into the frame buffer
		- CAUTION: when storing the pixels into the frame buffer, the order is blue, green, and red
		- NOT red, green, and blue !!!
	*/
	int bufInd = 0;
	for (int j = 0; j < yres; ++j) {
		for (int i = 0; i < xres; ++i) {
			int ind = ARRAY(i, j);

			unsigned char blue = (unsigned char)(pixelbuffer[ind].blue >> 4);
			unsigned char green = (unsigned char)(pixelbuffer[ind].green >> 4);
			unsigned char red = (unsigned char)(pixelbuffer[ind].red >> 4);

			framebuffer[bufInd] = blue;
			bufInd++;
			framebuffer[bufInd] = green;
			bufInd++;
			framebuffer[bufInd] = red;
			bufInd++;
		}
	}

	return GZ_SUCCESS;
}


/***********************************************/
/* HW2 methods: implement from here */

int GzRender::GzPutAttribute(int numAttributes, GzToken* nameList, GzPointer* valueList)
{
	/* HW 2.1
	-- Set renderer attribute states (e.g.: GZ_RGB_COLOR default color)
	-- In later homeworks set shaders, interpolaters, texture maps, and lights
	*/
	for (int i = 0; i < numAttributes; ++i) {
		GzToken token = nameList[i];
		if (token == GZ_RGB_COLOR) {
			GzColor* colorPtr = (GzColor*)valueList[i];
			// set color
			flatcolor[RED] = (*colorPtr)[RED];
			flatcolor[GREEN] = (*colorPtr)[GREEN];
			flatcolor[BLUE] = (*colorPtr)[BLUE];
		}
		else if (token == GZ_INTERPOLATE) {
			int* modePtr = (int*)valueList[i];
			interp_mode = *modePtr;
		}
		else if (token == GZ_DIRECTIONAL_LIGHT) {
			if (numAttributes > MAX_LIGHTS)
				return GZ_FAILURE;
			else {
				GzLight* dirLightPtr = (GzLight*)valueList[i];
				lights[i] = *dirLightPtr;
				numlights++;
			}
		}
		else if (token == GZ_AMBIENT_LIGHT) {
			GzLight* ambLightPtr = (GzLight*)valueList[i];
			ambientlight = *ambLightPtr;
		}
		else if (token == GZ_AMBIENT_COEFFICIENT) {
			GzColor* ambPtr = (GzColor*)valueList[i];
			Ka[RED] = (*ambPtr)[RED];
			Ka[GREEN] = (*ambPtr)[GREEN];
			Ka[BLUE] = (*ambPtr)[BLUE];
		}
		else if (token == GZ_DIFFUSE_COEFFICIENT) {
			GzColor* difPtr = (GzColor*)valueList[i];
			Kd[RED] = (*difPtr)[RED];
			Kd[GREEN] = (*difPtr)[GREEN];
			Kd[BLUE] = (*difPtr)[BLUE];
		}
		else if (token == GZ_SPECULAR_COEFFICIENT) {
			GzColor* specPtr = (GzColor*)valueList[i];
			Ks[RED] = (*specPtr)[RED];
			Ks[GREEN] = (*specPtr)[GREEN];
			Ks[BLUE] = (*specPtr)[BLUE];
		}
		else if (token == GZ_DISTRIBUTION_COEFFICIENT) {
			float* sPtr = (float*)valueList[i];
			spec = *sPtr;
		}
		else if (token == GZ_TEXTURE_MAP) {
			GzTexture tex_ptr = (GzTexture)valueList[i];
			tex_fun = tex_ptr;
		}
		else {
			return GZ_FAILURE;
		}


	}

	return GZ_SUCCESS;
}

int GzRender::GzPutTriangle(int	numParts, GzToken* nameList, GzPointer* valueList)
/* numParts - how many names and values */
{
	/* HW 2.2
	-- Pass in a triangle description with tokens and values corresponding to
		  GZ_NULL_TOKEN:		do nothing - no values
		  GZ_POSITION:		3 vert positions in model space
	-- Invoke the rastrizer/scanline framework
	-- Return error code
	*/
	float v1[4], v2[4], v3[4];
	float n1[4], n2[4], n3[4];
	float t1[2], t2[2], t3[2];

	for (int i = 0; i < numParts; ++i) {
		GzToken token = nameList[i];

		if (token == GZ_NULL_TOKEN) {
			// do nothing
		}
		else if (token == GZ_POSITION) {
			// get vertex coordinates			
			GzCoord* vertices = (GzCoord*)valueList[i];
			float oldv1[] = { vertices[0][0] , vertices[0][1], vertices[0][2] , 1 };
			float oldv2[] = { vertices[1][0] , vertices[1][1], vertices[1][2] , 1 };
			float oldv3[] = { vertices[2][0] , vertices[2][1], vertices[2][2] , 1 };

			// apply top transform of matrix to vertex
			matrixMultiplyVector(Ximage[matlevel - 1], oldv1, v1);
			matrixMultiplyVector(Ximage[matlevel - 1], oldv2, v2);
			matrixMultiplyVector(Ximage[matlevel - 1], oldv3, v3);

			v1[X] /= v1[3];
			v1[Y] /= v1[3];
			v1[Z] /= v1[3];

			v2[X] /= v2[3];
			v2[Y] /= v2[3];
			v2[Z] /= v2[3];

			v3[X] /= v3[3];
			v3[Y] /= v3[3];
			v3[Z] /= v3[3];
		}
		else if (token == GZ_NORMAL) {
			GzCoord* normals = (GzCoord*)valueList[i];

			float oldn1[] = { normals[0][0] , normals[0][1], normals[0][2] , 1 };
			float oldn2[] = { normals[1][0] , normals[1][1], normals[1][2] , 1 };
			float oldn3[] = { normals[2][0] , normals[2][1], normals[2][2] , 1 };

			// convert normals from model space to image space
			matrixMultiplyVector(Xnorm[matlevel - 1], oldn1, n1);
			matrixMultiplyVector(Xnorm[matlevel - 1], oldn2, n2);
			matrixMultiplyVector(Xnorm[matlevel - 1], oldn3, n3);
		}
		else if (token == GZ_TEXTURE_INDEX) {
			GzTextureIndex* textures = (GzTextureIndex*)valueList[i];
			// Transform u, v coordinates of vertices to perspective space
			t1[0] = perspectiveTransform(textures[0][0], v1[Z]);
			t1[1] = perspectiveTransform(textures[0][1], v1[Z]);

			t2[0] = perspectiveTransform(textures[1][0], v2[Z]);
			t2[1] = perspectiveTransform(textures[1][1], v2[Z]);

			t3[0] = perspectiveTransform(textures[2][0], v3[Z]);
			t3[1] = perspectiveTransform(textures[2][1], v3[Z]);
		}
	}

	
	// Evaluate color at vertices
	float c1flat[] = { 0.0, 0.0, 0.0 };
	getColor(n1, c1flat, Kd, Ks, Ka);

	float c1[] = { 0.0, 0.0, 0.0 };
	float c2[] = { 0.0, 0.0, 0.0 };
	float c3[] = { 0.0, 0.0, 0.0 };
	float dummyK[] = { 1.0, 1.0, 1.0 };

	if (tex_fun != nullptr) {
		getColor(n1, c1, dummyK, dummyK, dummyK);
		getColor(n2, c2, dummyK, dummyK, dummyK);
		getColor(n3, c3, dummyK, dummyK, dummyK);
	}
	else {
		getColor(n1, c1, Kd, Ks, Ka);
		getColor(n2, c2, Kd, Ks, Ka);
		getColor(n3, c3, Kd, Ks, Ka);
	}

	// Perform LEE

	// find the minimum and maximum coordinates for the bounding box
	int minX = max(0, floor(min(min(v1[X], v2[X]), v3[X])));
	int minY = max(0, floor(min(min(v1[Y], v2[Y]), v3[Y])));

	int maxX = min(xres, ceil(max(max(v1[X], v2[X]), v3[X])));
	int maxY = min(yres, ceil(max(max(v1[Y], v2[Y]), v3[Y])));


	// edge equations
	float a = v1[Y] - v2[Y];
	float b = v2[X] - v1[X];
	float c = -((v1[Y] - v2[Y]) * v1[X] + (v2[X] - v1[X]) * v1[Y]);
	float d = v2[Y] - v3[Y];
	float e = v3[X] - v2[X];
	float f = -((v2[Y] - v3[Y]) * v2[X] + (v3[X] - v2[X]) * v2[Y]);
	float g = v3[Y] - v1[Y];
	float h = v1[X] - v3[X];
	float i = -((v3[Y] - v1[Y]) * v3[X] + (v1[X] - v3[X]) * v3[Y]);

	// for every pixel in bounding box
	for (int y = minY; y <= maxY; y++) {
		for (int x = minX; x <= maxX; x++) {
			// calculate edge equations
			float e0 = a * x + b * y + c;
			float e1 = d * x + e * y + f;
			float e2 = g * x + h * y + i;

			// if all edge values are same sign
			if ((e0 >= 0 && e1 >= 0 && e2 >= 0) || (e0 <= 0 && e1 <= 0 && e2 <= 0)) {
				// interpolate z
				float interpolatedZ = interpolate(v1, v2, v3, v1[Z], v2[Z], v3[Z], x, y);

				int ind = ARRAY(x, y);
				float currZ = pixelbuffer[ind].z;
				// check z-buffer
				if (interpolatedZ < pixelbuffer[ind].z) {
					pixelbuffer[ind].z = interpolatedZ;
					pixelbuffer[ind].visible = true;


					// interpolate u, v
					float interpolatedU = interpolate(v1, v2, v3, t1[0], t2[0], t3[0], x, y);
					float interpolatedV = interpolate(v1, v2, v3, t1[1], t2[1], t3[1], x, y);

					pixelbuffer[ind].u = interpolatedU;
					pixelbuffer[ind].v = interpolatedV;

					// transform u, v back to shading space
					float u = shadingTransform(interpolatedU, interpolatedZ);
					float v = shadingTransform(interpolatedV, interpolatedZ);
					
					// lookup texture color
					GzColor texColor;
					if (tex_fun != NULL) {
						tex_fun(u, v, texColor);
					}

					// Flat interpolation
					if (interp_mode == 0) {
						GzPut(x, y, ctoi(c1flat[RED]), ctoi(c1flat[GREEN]), ctoi(c1flat[BLUE]), 1, interpolatedZ);
					}
					// Gourad interpolation
					else if (interp_mode == 1) {
						float interpolatedR = interpolate(v1, v2, v3, c1[RED], c2[RED], c3[RED], x, y);
						float interpolatedG = interpolate(v1, v2, v3, c1[GREEN], c2[GREEN], c3[GREEN], x, y);
						float interpolatedB = interpolate(v1, v2, v3, c1[BLUE], c2[BLUE], c3[BLUE], x, y);

						GzPut(x, y, ctoi(texColor[RED] * interpolatedR), ctoi(texColor[BLUE] * interpolatedG), ctoi(texColor[GREEN] * interpolatedB), 1, interpolatedZ);
					}
					// Phong interpolation
					else if (interp_mode == 2) {
						float normal[] = {
							interpolate(v1, v2, v3, n1[X], n2[X], n3[X], x, y),
							interpolate(v1, v2, v3, n1[Y], n2[Y], n3[Y], x, y),
							interpolate(v1, v2, v3, n1[Z], n2[Z], n3[Z], x, y)
						};
						normalize(normal);

						float c[] = { 0.0, 0.0, 0.0 };
						if (tex_fun != nullptr) {
							getColor(normal, c, texColor, Ks, texColor);
						}
						else {
							getColor(normal, c, Kd, Ks, Ka);
						}

						GzPut(x, y, ctoi(c[RED]), ctoi(c[GREEN]), ctoi(c[BLUE]), 1, interpolatedZ);
					}
					else if (interp_mode == 3) {
						float normal[] = {
							interpolate(v1, v2, v3, n1[X], n2[X], n3[X], x, y),
							interpolate(v1, v2, v3, n1[Y], n2[Y], n3[Y], x, y),
							interpolate(v1, v2, v3, n1[Z], n2[Z], n3[Z], x, y)
						};
						normalize(normal);

						float c[] = { 0.0, 0.0, 0.0 };
						if (tex_fun != nullptr) {
							getColor(normal, c, texColor, Ks, texColor);
						}
						else {
							getColor(normal, c, Kd, Ks, Ka);
						}

						GzPut(x, y, ctoi(c[RED]), ctoi(c[GREEN]), ctoi(c[BLUE]), 1, interpolatedZ);
					}
					// No interpolation
					else {
						GzPut(x, y, ctoi(flatcolor[RED]), ctoi(flatcolor[GREEN]), ctoi(flatcolor[BLUE]), 1, interpolatedZ);
					}

				}
			}
		}
	}

	return GZ_SUCCESS;
}

// Project functions

bool GzRender::isValidPixel(int x, int y)
{
	return (
		(x > 0) &&
		(x < this->xres - 1) &&
		(y > 0) &&
		(y < this->yres - 1));
}

float GzRender::getGradient(int x, int y, GzPixel* pixelbuffer, int param)
{
	

	float g = 0.0;
	float l = 0.0;
	float g_contour = 0.0;
	if (isValidPixel(x, y)) {
		float A = pixelbuffer[ARRAY(x - 1, y - 1)].z;
		float B = pixelbuffer[ARRAY(x, y - 1)].z;
		float C = pixelbuffer[ARRAY(x + 1, y - 1)].z;
		float D = pixelbuffer[ARRAY(x - 1, y)].z;
		float curr = pixelbuffer[ARRAY(x, y)].z;
		float E = pixelbuffer[ARRAY(x + 1, y)].z;
		float F = pixelbuffer[ARRAY(x - 1, y + 1)].z;
		float G = pixelbuffer[ARRAY(x, y + 1)].z;
		float H = pixelbuffer[ARRAY(x + 1, y + 1)].z;
		
		g = (abs(A + 2 * B + C - F - 2 * G - H) + abs(C + 2 * E + H - A - 2 * D - F)) / 8;
		l = (8 * curr - A - B - C - D - E - F - G - H) / 3;
	}

	if (param == 1) return g;
	else if (param == 2) return l;

}

int GzRender::GzSobelEdgeDetection()
{
	float thresholdZ = 3000000;
	float thresholdZ2 = 70000;
	float thresholdN = 50;
	float thresholdDark = 0.25;
	std::random_device rd;
	std::mt19937 gen(rd()); 
	std::uniform_int_distribution<> dis(1, 100);

	// for every pixel
	for (int i = 0; i < this->xres; i++)
	{
		for (int j = 0; j < this->yres; j++)
		{

			if (!pixelbuffer[ARRAY(i, j)].visible) continue;
			GzIntensity r, g, b, a;
			GzDepth z;
			GzGet(i, j, &r, &g, &b, &a, &z);

			int u = (int)(pixelbuffer[ARRAY(i, j)].u * 4000);
			int v = (int)(pixelbuffer[ARRAY(i, j)].v * 4000);
			float dark = (0.2126f * itoc(r) + 0.7152f * itoc(g) + 0.0722f * itoc(b));
			int k = floor((20 - 6) * dark + 6);
			
			
			float zGradient = getGradient(i, j, pixelbuffer, 1);
			float zGradient2 = getGradient(i, j, pixelbuffer, 2);
			float contourGradient = getGradient(i, j, pixelbuffer, 3);
			// first order edge detection
			if (zGradient > thresholdZ)
			{
				GzPut(i, j, ctoi(0.0), ctoi(0.0), ctoi(0.0), a, z);				
			}
			// second order edge detection
			else if (zGradient2 > thresholdZ2)
			{
				GzPut(i, j, ctoi(0.0), ctoi(0.0), ctoi(0.0), a, z);
			}
			
			// cross-hatch shading
			else if ((u % k == 0 || v % k == 0) && dark < 0.5) {
				//float col = std::pow(4 * dark, 2);
				GzPut(i, j, ctoi(0.0), ctoi(0.0), ctoi(0.0), a, z);
			}

			// stipple
			else if (dis(gen) > (dark * 350) && (int)dis(gen) % k == 0) {
				//float col = std::pow(4 * dark, 2);
				GzPut(i, j, ctoi(0.0), ctoi(0.0), ctoi(0.0), a, z);
			}
			// else make white
			else {
				GzPut(i, j, ctoi(1.0f), ctoi(1.0f), ctoi(1.0f), a, z);
			}
		}
	}
	return GZ_SUCCESS;
}