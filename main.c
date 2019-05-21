#include <stdio.h>
#include <math.h>
#include <stdlib.h>
#include <limits.h>
#include <string.h>
#include <unistd.h>
#include "tga.h"
#include "model.h"

#define WIDTH 1400
#define HEIGHT 1400
#define DEPTH 255

typedef int Vec2i[2];
typedef int Vec3i[3];

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

typedef struct {
	double m[4][4];
} Matrix44;

typedef struct {
	double m[4][1];
} Matrix41;

typedef struct {
	double m[4][2];
} Matrix42;

void m442v(Vec3 *vector, Matrix44 *m) {
	//Vec3 vector = {m.m[0][0]/m.m[3][0], m.m[1][0]/m.m[3][0], m.m[2][0]/m.m[3][0]};
	*vector[0] = m->m[0][0]/m->m[3][0];
	*vector[1] = m->m[1][0]/m->m[3][0];
	*vector[2] = m->m[2][0]/m->m[3][0];
}

void m412v(Vec3 *vector, Matrix41 *m) {
	//Vec3 vector = {m.m[0][0]/m.m[3][0], m.m[1][0]/m.m[3][0], m.m[2][0]/m.m[3][0]};
	*vector[0] = m->m[0][0]/m->m[3][0];
	*vector[1] = m->m[1][0]/m->m[3][0];
	*vector[2] = m->m[2][0]/m->m[3][0];
}


void initMatrixIdentity(Matrix44 *matrix) {
	for(int i = 0; i < 4; i++) {
    	for(int j = 0; j < 4; j++) {
    		if(i == j) {
    			matrix->m[i][j] = 1;
    		} else {
    			matrix->m[i][j] = 0;
    		}
    	}
    }
}

void v2m(Vec3 *vector, Matrix41 *m) {
    m->m[0][0] = *vector[0];
    m->m[1][0] = *vector[1];
    m->m[2][0] = *vector[2];
    m->m[3][0] = 1.f;
}

void viewport(int x, int y, int w, int h, Matrix44 *matrix) {
    for(int i = 0; i < 4; i++) {
    	for(int j = 0; j < 4; j++) {
    		if(i == j) {
    			matrix->m[i][j] = 1;
    		} else {
    			matrix->m[i][j] = 0;
    		}
    	}
    }
    matrix->m[0][3] = x+w/2.f;
    matrix->m[1][3] = y+h/2.f;
    matrix->m[2][3] = DEPTH/2.f;

    matrix->m[0][0] = w/2.f;
    matrix->m[1][1] = h/2.f;
    matrix->m[2][2] = DEPTH/2.f;
}

void mulMatrix44x44(Matrix44 *distMatrix, Matrix44 *firstMatrix, Matrix44 *secondMatrix) {
    for (int i=0; i<4; i++) {
        for (int j=0; j<4; j++) {
            distMatrix->m[i][j] = 0.f;
            for (int k=0; k<4; k++) {
                distMatrix->m[i][j] += firstMatrix->m[i][k]*secondMatrix->m[k][j];
            }
        }
    }
}

void mulMatrix44x42(Matrix42 *distMatrix, Matrix44 *firstMatrix, Matrix42 *secondMatrix) {
    for (int i=0; i<4; i++) {
        for (int j=0; j<2; j++) {
            distMatrix->m[i][j] = 0.f;
            for (int k=0; k<4; k++) {
                distMatrix->m[i][j] += firstMatrix->m[i][k]*secondMatrix->m[k][j];
            }
        }
    }
}

void mulMatrix42x41(Matrix42 *distMatrix, Matrix44 *firstMatrix, Matrix42 *secondMatrix) {
	for (int i=0; i<4; i++) {
        for (int j=0; j<1; j++) {
            distMatrix->m[i][j] = 0.f;
            for (int k=0; k<2; k++) {
                distMatrix->m[i][j] += firstMatrix->m[i][k]*secondMatrix->m[k][j];
            }
        }
    }
}

void mulMatrix44x41(Matrix41 *distMatrix, Matrix44 *firstMatrix, Matrix41 *secondMatrix) {
	for (int i=0; i<4; i++) {
        for (int j=0; j<1; j++) {
            distMatrix->m[i][j] = 0.f;
            for (int k=0; k<4; k++) {
                distMatrix->m[i][j] += firstMatrix->m[i][k]*secondMatrix->m[k][j];
            }
        }
    }
}

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

void swapf(double *a, double *b) {
	double t = *a;
	*a = *b;
	*b = t;
}

void newTriangle(tgaImage *image, int x0, int y0, int z0, int x1, int y1, int z1,
								int x2, int y2, int z2,
								tgaColor color, Vec3 uv0, Vec3 uv1, Vec3 uv2, Model *model, int *zBuffer) {
	if((y0 == y1) && (y0 == y2)) {
		return;
	}

	if(y0 > y1) {
		swap(&y0, &y1);
		swap(&x0, &x1);
		swap(&z0, &z1);

		for(int i = 0; i < 2; i++) {
			swapf(&uv0[i], &uv1[i]);
		}
	}

	if(y0 > y2) {
		swap(&y0, &y2);
		swap(&x0, &x2);
		swap(&z0, &z2);

		for(int i = 0; i < 2; i++) {
			swapf(&uv0[i], &uv2[i]);
		}
	}
	
	if(y1 > y2) {
		swap(&y1, &y2);
		swap(&x1, &x2);
		swap(&z1, &z2);

		for(int i = 0; i < 2; i++) {
			swapf(&uv1[i], &uv2[i]);
		}
	}

	int totalHeight = y2-y0;
	for(int i = 0; i < totalHeight; i++) {
		int secondHalf = i > y1-y0 || y1 == y0;
		int segmentHeight = secondHalf ? y2-y1 : y1-y0;
		double alpha = (double)i/totalHeight;
		double beta = (double)(i-(secondHalf ? y1-y0 : 0))/segmentHeight;
		
		Vec3i A = {x0 + ((x2 - x0)*alpha), y0 + ((y2 - y0)*alpha), z0 + ((z2 - z0)*alpha)};
		Vec3i B = {x1 + ((x2 - x1)*beta), y1 + ((y2 - y1)*beta), z1 + ((z2 - z1)*beta)};
		if(!secondHalf) {
			B[0] = x0 + ((x1 - x0)*beta);
			B[1] = y0 + ((y1 - y0)*beta);
			B[2] = z0 + ((z1 - z0)*beta);
		}

        Vec3 uvA = {uv0[0] + (uv2[0] - uv0[0])*alpha, uv0[1] + (uv2[1] - uv0[1])*alpha};
        Vec3 uvB = {uv0[0] + (uv1[0] - uv0[0])*beta, uv0[1] + (uv1[1] - uv0[1])*beta};

        if(secondHalf) {
        	for(int i = 0; i < 2; i++) {
        		uvB[i] = uv1[i] + (uv2[i] - uv1[i])*beta;
        	}
        }
		if(A[0] > B[0]) {
			for(int i = 0; i < 3; i++) {
				swap(&A[i], &B[i]);
				swapf(&uvA[i], &uvB[i]);
			}
		}

		for(int j = A[0]; j <=B[0]; j++) {
			double phi = B[0]==A[0] ? 1. : (double)(j-A[0])/(double)(B[0]-A[0]);
			Vec3 P = {A[0] + ((B[0] - A[0])*phi), A[1] + ((B[1] - A[1])*phi), A[2] + ((B[2] - A[2])*phi)};
			int idx = P[0] + P[1] * WIDTH;

			if(zBuffer[idx] < P[2]) {
				zBuffer[idx] = P[2];
				#ifdef DEBUG
				printf("SET PIXEL X = %d; Y = %d; Z = %f\n", j, i, P[2]);
				#endif

				Vec3 uvP = {uvA[0] + (uvB[0] - uvA[0])*phi, uvA[1] + (uvB[1] - uvA[1])*phi};
				color = getDiffuseColor(model, &uvP);

				tgaSetPixel(image, P[0], P[1], color);
			}
		}
	}
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
			xa = x1;
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

			#ifdef DEBUG
			printf("INDEX ZBUFFER: %d\n", idx);
			#endif

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
			xb = x1;
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

			#ifdef DEBUG
			printf("INDEX ZBUFFER: %d\n", idx);
			#endif

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

void normalizeVec3(Vec3 vector) {
	double fraction = 1/sqrt((vector)[0] * (vector)[0] + (vector)[1] * (vector)[1] + (vector)[2] * (vector)[2]);

	for(int i = 0; i < 3; i++) {
		vector[i] *= fraction;
	}
}

double norm(Vec3 *vector) {
	return sqrt((*vector[0]) * (*vector[0]) + (*vector[1]) * (*vector[1]) + (*vector[2]) * (*vector[2]));
}

void lookat(Matrix44 *ModelView, Vec3 eye, Vec3 center, Vec3 up) {
    Vec3 z = {eye[0]-center[0], eye[1]-center[1], eye[2]-center[2]};
    normalizeVec3(z);

    Vec3 x = {up[1] * z[2] - up[2] * z[1], up[2] * z[0] - up[0] * z[2], up[0] * z[1] - up[1] * z[0]}; //Cross (up ^ z)
    normalizeVec3(x);

    Vec3 y = {z[1] * x[2] - z[2] * x[1], z[2]  *x[0] - z[0] * x[2], z[0] * x[1] - z[1] * x[0]}; //Cross (z ^ x)
    normalizeVec3(y);

    Matrix44 Minv;
    initMatrixIdentity(&Minv);

    Matrix44 Tr;
    initMatrixIdentity(&Tr);
    for (int i=0; i<3; i++) {
        Minv.m[0][i] = x[i];
        Minv.m[1][i] = y[i];
        Minv.m[2][i] = z[i];
        Tr.m[i][3] = -center[i];
    }

    mulMatrix44x44(ModelView, &Minv, &Tr);

}
void meshgrid(tgaImage *image, Model *model, char *argv) {
	double lightIntensity = 1;
	double color = 255;
	Vec3 lightDirection = {0, 0, 1};
	//Vec3 eye = {1, 0, 2};
	//Vec3 center = {1.5, 0, 2.5}; //rotated
	Vec3 eye = {0.5, 0.2, 2};
	Vec3 center = {0.5, 0.2, 2.5};
	Vec3 up = {0, 4, 1};
	Vec3 eye_minus_center = {eye[0]-center[0], eye[1]-center[1], eye[2]-center[2]};
	normalizeVec3(eye_minus_center);
	Matrix44 ModelView;
	initMatrixIdentity(&ModelView);
	lookat(&ModelView, eye, center, up);
	Vec3 *vertices[3];
	int zBuffer[HEIGHT * WIDTH];
	Vec3 *uv[3];

	Matrix44 Projection;
	initMatrixIdentity(&Projection);

	Matrix44 ViewPort;
	viewport(WIDTH/8, HEIGHT/8, WIDTH*3/4, HEIGHT*3/4, &ViewPort);
	#ifdef DEBUG
	printf("ViewPort:\n");
			for(int i = 0; i < 4; i++) {
				printf("%f %f %f %f\n", ViewPort.m[i][0], ViewPort.m[i][1], ViewPort.m[i][2], ViewPort.m[i][3]);
			}
	#endif
	Projection.m[3][2] = -1.f/norm(&eye_minus_center);

	for(int i = 0; i < WIDTH*HEIGHT; i++) {
		zBuffer[i] = INT_MIN;
	}

	for (unsigned i = 0; i < model->nface; ++i) {
		int screen_coords[3][3];
		for (unsigned j = 0; j < 3; ++j) {
			Vec3 *v = &(model->vertices[model->faces[i][3*j]]);
			vertices[j] = &(model->vertices[model->faces[i][3*j]]);
			/*screen_coords[j][0] = round(((*v)[0] + 1) * image->width / 2);
			screen_coords[j][1] = round((1 - (*v)[1]) * image->height / 2);
			screen_coords[j][2] = round((1 - (*v)[2]) * DEPTH/2);*/
			Vec3 vertices_vec = {(*v)[0], (*v)[1], (*v)[2]};
			#ifdef DEBUG
			printf("vertices_vec %f %f %f\n", vertices_vec[0], vertices_vec[1], vertices_vec[2]);
			#endif




			Matrix44 finally_screen_coords_pre;
			mulMatrix44x44(&finally_screen_coords_pre, &ViewPort, &Projection);
			#ifdef DEBUG
			printf("finally_screen_coords_4x4matrix:\n");
			for(int i = 0; i < 4; i++) {
				printf("%f %f %f %f\n", finally_screen_coords_pre.m[i][0], finally_screen_coords_pre.m[i][1], finally_screen_coords_pre.m[i][2], finally_screen_coords_pre.m[i][3]);
			}
			#endif

			Matrix44 finally_screen_coords;
			mulMatrix44x44(&finally_screen_coords, &finally_screen_coords_pre, &ModelView);



			Matrix41 finally_screen_coords_final;
			Matrix41 vertices_matrix;
			v2m(&vertices_vec, &vertices_matrix);
			for(int i = 0; i < 4; i++) {
				if(i < 3) {
					vertices_matrix.m[i][0] = vertices_vec[i];
				} else {
					vertices_matrix.m[i][0] = 1.f;
				}
			}
			#ifdef DEBUG
			printf("vertices_matrix:\n");
			for(int i = 0; i < 4; i++) {
				printf("%f\n", vertices_matrix.m[i][0]);
			}
			#endif
			mulMatrix44x41(&finally_screen_coords_final, &finally_screen_coords, &vertices_matrix);
			#ifdef DEBUG
			printf("finally_screen_coords_final4x1:\n");
			for(int i = 0; i < 4; i++) {
				printf("%f\n", finally_screen_coords_final.m[i][0]);
			}
			#endif

			Vec3 finalVecEnd;
			m412v(&finalVecEnd, &finally_screen_coords_final);
			finalVecEnd[0] = finally_screen_coords_final.m[0][0]/finally_screen_coords_final.m[3][0];
			finalVecEnd[1] = finally_screen_coords_final.m[1][0]/finally_screen_coords_final.m[3][0];
			finalVecEnd[2] = finally_screen_coords_final.m[2][0]/finally_screen_coords_final.m[3][0];

			#ifdef DEBUG
			printf("finalVecEnd:\n");
			for(int i = 0; i < 3; i++) {
				printf("%f\n", finalVecEnd[i]);
			}
			#endif

			Matrix44 test;
			Matrix44 testres;
			initMatrixIdentity(&test);
			mulMatrix44x44(&testres, &test, &test);

			screen_coords[j][0] = round(finalVecEnd[0]);
			screen_coords[j][1] = round(finalVecEnd[1]);
			screen_coords[j][2] = round(finalVecEnd[2]);
			#ifdef DEBUG
			printf("COORDS %f %f %f\n", finalVecEnd[0], finalVecEnd[1], finalVecEnd[2]);
			printf("COORDS SCREEN %d %d %d\n", screen_coords[j][0], screen_coords[j][1], screen_coords[j][2]);
			#endif

			uv[j] = getDiffuseUV(model, i, j);
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

		if(colorIntensity >= 0) {
			double colorCode = round(fabs(colorIntensity) * color);

			#ifdef DEBUG
			printf("tgaColor: %f\n", colorCode);
			#endif
			tgaColor randColor = tgaRGB(colorCode, colorCode, colorCode);

			#ifdef DEBUG
			for(int j = 0; j < 3; j++) {
				printf("COORDS SCREEN TRIANGLE j:%d: %d %d %d\n", j, screen_coords[j][0], screen_coords[j][1], screen_coords[j][2]);
			}
			#endif

			triangle(image, screen_coords[0][0], screen_coords[0][1], screen_coords[0][2],
				screen_coords[1][0], screen_coords[1][1], screen_coords[1][2],
				screen_coords[2][0], screen_coords[2][1], screen_coords[2][2],
				randColor, zBuffer);
			/*newTriangle(image, screen_coords[0][0], screen_coords[0][1], screen_coords[0][2],
				screen_coords[1][0], screen_coords[1][1], screen_coords[1][2],
				screen_coords[2][0], screen_coords[2][1], screen_coords[2][2],
				randColor, *uv[0], *uv[1], *uv[2], model, zBuffer);*/
			
		}
		#ifdef DEBUG
		printf("debug END!!!!\n");
		#endif

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

	loadDiffuseMap(model, "obj/cat_diff.tga");

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
