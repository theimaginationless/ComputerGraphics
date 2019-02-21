#include <stdio.h>
#include <math.h>
#include "model.h"
#include "tga.h"

#define WIDTH 800
#define HEIGHT 800

void swap(int *a, int *b);
int iabs(int a);

/*
* Using Digital Differential Analyzer algorihm
* to draw interval connecting (x0, y0) with (x1, y1)
* on image using color
*/
void line (tgaImage *image, 
           int x0, int y0,
           int x1, int y1,
           tgaColor color);

void meshgrid(tgaImage *image, Model *model) {
    tgaColor white = tgaRGB(255, 255, 255);

    for (unsigned i = 0; i < model->nface; ++i) {
        int screen_coords[3][2];
        for (unsigned j = 0; j < 3; ++j) {
            Vec3 *v = &(model->vertices[model->faces[i][3*j]]);
            screen_coords[j][0] = ((*v)[0] + 1) * image->width / 2;
            screen_coords[j][1] = (1 - (*v)[1]) * image->height / 2;
        }

        for (unsigned j = 0; j < 3; ++j) {
            line(image, screen_coords[j][0],screen_coords[j][1],
                screen_coords[(j+1)%3][0], screen_coords[(j+1)%3][1],white);
        }
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

    meshgrid(image, model);

    if (-1 == tgaSaveToFile(image, argv[2])) {
        perror("tgaSateToFile");
        rv = -1;
    }

    tgaFreeImage(image);    
    return rv;
}

void line (tgaImage *image, 
           int x0, int y0,
           int x1, int y1,
           tgaColor color)
{
    int steep = 0;
    if (iabs(y1 - y0) > iabs(x1 - x0)) {
        steep = 1;
        swap(&x0, &y0);
        swap(&x1, &y1);
    }

    if (x0 > x1) {
        swap(&x0, &x1);
        swap(&y0, &y1);
    }

    int x;
    double y;
    double k = ((double)(y1 - y0))/(x1 - x0);
    for (x = x0, y = y0; x <= x1; ++x, y += k) {
        if (steep) {
            tgaSetPixel(image, (unsigned int)y, (unsigned int)x, color);
        } else {
            tgaSetPixel(image, (unsigned int)x, (unsigned int)y, color);
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
