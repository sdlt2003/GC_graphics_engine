
#ifndef TRANSFORMATIONS_H
#define TRANSFORMATIONS_H

#include <stdio.h>
#include <string.h>
#include <math.h>
#include <GL/glut.h>
#include "utils.h"

void scale(int sign);
void translate(double x, double y, double z);
void delete_object(triobj *optr);
void rotate_around_object(char axis, int sign);
void rotate(double x, double y, double z);
void x_trfm(int sign);
void y_trfm(int sign);
void z_trfm(int sign);

#endif
