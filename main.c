#include <stdio.h>
#include <math.h>
#include <stdlib.h>
#include <limits.h>
#include <unistd.h>
#include "tga.h"
#include "model.h"

#define WIDTH 1400
#define HEIGHT 1400
#define DEPTH 255

void swap(int *a, int *b);
int iabs(int a);
int sign(double x);


/*
* Using Bresenham's algorithm
* to draw interval connecting (x0, y0) with (x1, y1)
* on image using color
*/
void BRLine (tgaImage *image, 
		   int x0, int y0,
		   int x1, int y1,
		   tgaColor color);

Model *scaleModel(Model *model, double scale) {
	for(unsigned i = 0; i < model->nvert; i++) {
		for(unsigned j = 0; j < 3; j++) {
			(model->vertices[i])[j] = (model->vertices[i])[j] * scale;
		}
	}

	return model;
}

Model *offsetModel(Model *model, double x, double y, double z) {
	for(unsigned i = 0; i < model->nvert; i++) {
		(model->vertices[i])[0] += x;
		(model->vertices[i])[1] += y;
		(model->vertices[i])[2] += z;
	}

	return model;
}

double getZCoord(int x0, int y0, int z0,
				int x1, int y1, int z1,
				int x2, int y2, int z2,
				int x, int y) {
					int a = y0*(z1 - z2) + y1*(z2 - z0) + y2*(z0 - z1);
					int b = z0*(x1 - x2) + z1*(x2 - x0) + z2*(x0 - x1);
					int c = x0*(y1 - y2) + x1*(y2 - y0) + x2*(y0 - y1);
					
					#ifdef DEBUG
					printf("getZCoord: a = %d; b = %d; c = %d\n", a, b, c);
					#endif

					if(c == 0) {
						return 0;
					}

					double z = (a*(x - x0) + b*(y - y0) - c*z0)/(double)c;

					#ifdef DEBUG
					printf("getZCoord: z = %f\n", z);
					#endif

					return z;
				}

void triangle(tgaImage *image, int x0, int y0, int z0, int x1, int y1, int z1,
								int x2, int y2, int z2, tgaColor color, int *zBuffer) {
	if((y0 == y1) && (y0 == y2)) {
		return;
	}

	if(y0 > y1) {
		swap(&y0, &y1);
		swap(&x0, &x1);
		swap(&z0, &z1);
	}

	if(y0 > y2) {
		swap(&y0, &y2);
		swap(&x0, &x2);
		swap(&z0, &z2);
	}
	
	if(y1 > y2) {
		swap(&y1, &y2);
		swap(&x1, &x2);
		swap(&z1, &z2);
	}

	#ifdef DEBUG
	printf("y0 y1 y2 %d %d %d\n", y0, y1, y2);
	printf("DEBUG x0 x1 x2 %d %d %d\n", x0, x1, x2);
	#endif

	int xa, xb, y;

	for(y = y0; y <= y1; y++) {
		

		if((y1 - y0) != 0) {
			xa = round(x0 + (x1-x0) * (double)(y - y0)/(y1 - y0));
		} else {
			xa = x0;
		}
		if((y2 - y0) != 0) {
			xb = round(x0 + (x2-x0) * (double)(y - y0)/(y2 - y0));
		} else {
			xb = x0;
		}

		if(xa > xb) {
			swap(&xa, &xb);
		}
		#ifdef DEBUG
		printf("DEBUG xa xb %d %d\n", xa, xb);
		#endif

		for(int x = xa; x <= xb; x++) {
			int idx = x + y * WIDTH;
			int z = round(getZCoord(x0, y0, z0, x1, y1, z1, x2, y2, z2, x, y));

			if(zBuffer[idx] < z) {
				zBuffer[idx] = z;
				#ifdef DEBUG
				printf("SET PIXEL X = %d; Y = %d; Z = %d\n", x, y, z);
				#endif

				tgaSetPixel(image, x, y, color);
			}
		}
		
	}
	/*
	* Redefinition upper vertices (x, y)
	*/

	y0 = y1 = y;
	x0 = xa;
	x1 = xb;

	#ifdef DEBUG
	printf("ADDITION x0 y0 %d %d; x1 y1 %d %d; x2 y2 %d %d\n", x0, y0, x1, y1, x2, y2);
	#endif

	for(int y = y0; y <= y2; y++) {
		if((y2 - y0) != 0) {
			xa = round(x0 + (x2 - x0)*(double)(y - y0)/(y2 - y0));

			#ifdef DEBUG
			printf("xb => y2 - y0 != 0 => %d\n", xb);
			#endif
		} else {
			xa = x0;
		}
		if((y2 - y1) != 0) {
			xb = round(x1 + (x2 - x1)*(double)(y - y1)/(y2 - y1));

			#ifdef DEBUG
			printf("xb => y2 - y1 != 0 => %d\n", xb);
			#endif

		} else {
			xb = x0;
		}

		if(xa > xb) {
			swap(&xa, &xb);
		}

		#ifdef DEBUG
		printf("DEBUG ADDITION xa xb %d %d\n", xa, xb);
		#endif

		for(int x = xa; x <= xb; x++) {

			int idx = x + y * WIDTH;
			int z = round(getZCoord(x0, y0, z0, x1, y1, z1, x2, y2, z2, x, y));

			if(zBuffer[idx] < z) {
				zBuffer[idx] = z;
				#ifdef DEBUG
				printf("SET PIXEL X = %d; Y = %d; Z = %d\n", x, y, z);
				#endif

				tgaSetPixel(image, x, y, color);
			}
		}
	}
}

double newgetAngleNormal(Vec3 lightDirection, double x0, double y0, double z0,
		double x1, double y1, double z1,
		double x2, double y2, double z2) {
			Vec3 wc1 = {x0 - x1, y0 - y1, z0 - z1};
			Vec3 wc2 = {x1 - x2, y1 - y2, z1 - z2};

			// ^

			Vec3 result = {
				wc1[1]*wc2[2] - wc1[2]*wc2[1],
				wc1[2]*wc2[0] - wc1[0]*wc2[2],
				wc1[0]*wc2[1] - wc1[1]*wc2[0]};

			// normalize
			//for(int i = 0; i < 3; i++) {
			//	result[i] /= sqrt(result[0]*result[0] + result[1]*result[1] + result[2]*result[2]);
			//}
			double fracA = sqrt(result[0] * result[0] + result[1] * result[1] + result[2] * result[2]);
			double fracB = sqrt(lightDirection[0] * lightDirection[0] + lightDirection[1] * lightDirection[1] + lightDirection[2] * lightDirection[2]);
			double res = (result[0]*lightDirection[0] + result[1]*lightDirection[1] + result[2]*lightDirection[2])/(fracA*fracB);
			
			#ifdef DEBUG
			printf("getAngleNormal: angle: %f\n", res);
			#endif

			return res;
		}

double getAngleNormal(Vec3 lightDirection, double x0, double y0, double z0,
		double x1, double y1, double z1,
		double x2, double y2, double z2) {
			double a, b, c;
			a = y0*(z1 - z2) + y1*(z2 - z0) + y2*(z0 - z1);
			b = z0*(x1 - x2) + z1*(x2 - x0) + z2*(x0 - x1);
			c = x0*(y1 - y2) + x1*(y2 - y0) + x2*(y0 - y1);

			#ifdef DEBUG
			printf("getAngleNormal: Surface: %f, %f, %f\n", a, b, c);
			#endif

			if((a == b) && (b == c)) {
				return 1;
			}

			// normalize normal vector
			double lenFrac = sqrt(a*a + b*b + c*c);
			Vec3 normal = {a/lenFrac, b/lenFrac, c/lenFrac};

			#ifdef DEBUG
			printf("getAngleNormal: Normal: %f, %f, %f\n", normal[0], normal[1], normal[2]);
			#endif

			// cos(normal ^ lightDirection)
			double fracA = sqrt(normal[0] * normal[0] + normal[1] * normal[1] + normal[2] * normal[2]);
			double fracB = sqrt(lightDirection[0] * lightDirection[0] + lightDirection[1] * lightDirection[1] + lightDirection[2] * lightDirection[2]);
			double angle = (normal[0]*lightDirection[0] + normal[1]*lightDirection[1] + normal[2]*lightDirection[2]);

			#ifdef DEBUG
			printf("getAngleNormal: angle: %f\n", angle);
			#endif

			return angle;
		}

void meshgrid(tgaImage *image, Model *model, char *argv) {
	double lightIntensity = 1;
	double color = 255;
	Vec3 lightDirection = {0, 0, 1};
	Vec3 *vertices[3];
	int zBuffer[HEIGHT * WIDTH];

	for(int i = 0; i < WIDTH*HEIGHT; i++) {
		zBuffer[i] = INT_MIN;
	}

	for (unsigned i = 0; i < model->nface; ++i) {
		int screen_coords[3][3];
		for (unsigned j = 0; j < 3; ++j) {
			Vec3 *v = &(model->vertices[model->faces[i][3*j]]);
			vertices[j] = &(model->vertices[model->faces[i][3*j]]);
			screen_coords[j][0] = round(((*v)[0] + 1) * image->width / 2);
			screen_coords[j][1] = round((1 - (*v)[1]) * image->height / 2);
			screen_coords[j][2] = round((1 - (*v)[2]) * DEPTH/2);
		}
		
		double nCosAngle = getAngleNormal(lightDirection,
			(*vertices[0])[0], (*vertices[0])[1], (*vertices[0])[2],
			(*vertices[1])[0], (*vertices[1])[1], (*vertices[1])[2],
			(*vertices[2])[0], (*vertices[2])[1], (*vertices[2])[2]);

		double colorIntensity = lightIntensity * nCosAngle;

		#ifdef DEBUG
		printf("lightIntensity = %f; colorIntensity = %f\n", lightIntensity, colorIntensity);
		#endif

		/*
		* Note: If call triangle when colorIntensity <= 0, then we have lost poly???
		*/

		//if(colorIntensity <= 0) {
			int j = 0;
			double colorCode = ceil(fabs(colorIntensity) * color);

			#ifdef DEBUG
			printf("tgaColor: %f\n", colorCode);
			#endif
			tgaColor randColor = tgaRGB(colorCode, colorCode, colorCode);

			triangle(image, screen_coords[j][0], screen_coords[j][1], screen_coords[j][2],
				screen_coords[(j+1)%3][0], screen_coords[(j+1)%3][1], screen_coords[(j+1)%3][2],
				screen_coords[(j+2)%3][0], screen_coords[(j+2)%3][1], screen_coords[(j+2)%3][2], randColor, zBuffer); 
		//}

	}
}

int main(int argc, char **argv)
{
	int rv = 0;
	if (argc < 3) {
		fprintf(stderr, "Usage: %s infile outfile\n", argv[0]);
		return -1;
	}

	tgaImage *image = tgaNewImage(HEIGHT, WIDTH, RGB);
	Model *model = loadFromObj(argv[1]);

	/*
	* ./main obj/cat.obj cat.tga <SCALE> <XOFFSET> <YOFFSET> <ZOFFSET>
	* SCALE, XOFFSET, YOFFSET, ZOFFSET - fractional values.
	*/
	if(argc > 6) {
		scaleModel(model, strtod(argv[3], NULL));
		model = offsetModel(model, strtod(argv[4], NULL), strtod(argv[5], NULL), strtod(argv[6], NULL));
	}

	meshgrid(image, model, argv[2]);

	if (-1 == tgaSaveToFile(image, argv[2])) {
		perror("tgaSateToFile");
		rv = -1;
	}

	tgaFreeImage(image);    
	return rv;
}

void BRLine(tgaImage *image, 
		   int x0, int y0,
		   int x1, int y1,
		   tgaColor color)
{
	char steep = 0;
	if (iabs(y1-y0) > iabs(x1-x0)) {
		steep = 1;
		swap(&x0, &y0);
		swap(&x1, &y1);
	}

	if (x0 > x1) {
		swap(&x0, &x1);
		swap(&y0, &y1);
	}

	int dx = x1 - x0;
	int dy = y1 - y0;
	int dError = 2*iabs(dy);
	int error = 0;
	unsigned x, y;
	for (x = x0, y = y0; x <= x1; ++x) {
		if(steep) {
			tgaSetPixel(image, (unsigned int)y, (unsigned int)x, color);
		} else {
			tgaSetPixel(image, (unsigned int)x, (unsigned int)y, color);
		}
		error += dError;
		if(error > dx) {
			y += sign(dy);
			error -= 2*dx;
		}
	}
}

void swap(int *a, int *b) {
	int t = *a;
	*a = *b;
	*b = t;
}

int iabs(int a) {
	return (a >= 0) ? a : -a;
}

int sign(double x) {
	return (x >= 0) ? 1: -1;
}
