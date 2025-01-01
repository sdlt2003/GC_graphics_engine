
#ifndef CAMERA_H
#define CAMERA_H

#include "utils.h"
#include <GL/glut.h>
#include <stdio.h>
#include <string.h>
#include <math.h>

void calc_cam_matrix(triobj *first_cam_ptr);
void calc_projection_matrix();
void obtener_CSR_partiendo_de_M(GLdouble *M, GLdouble *MCSR);
void orientate_cam_to_obj();
void update_cam_params();
void init_camera();

#endif