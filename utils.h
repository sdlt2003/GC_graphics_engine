
#ifndef UTILS_H
#define UTILS_H

#include <stdio.h>
#include <string.h>
#include <math.h>
#include <GL/glut.h>

/////////////////
// ESTRUCTURAS //
/////////////////

// Almacena las coordenadas de un punto junto a sus coordenadas de textura.
typedef struct punto
{
    float x, y, z, u, v;
} punto;

// Almacena las coordenadas de un vector.
typedef struct vector
{
    float x, y, z;
} vector;

// Almacena las coordenadas de los vértices de un triángulo y su vector normal.
typedef struct triangulo
{
    punto p1, p2, p3;
    vector N;
} triangulo;

// Almacena matrices de transformación que se aplicarán a los objetos.
typedef struct mlist
{
    double m[16];
    struct mlist *next;
} mlist;

// Contiene la información de los triángulos a dibujar y una referencia a su matriz de transformación.
typedef struct triobj
{
    triangulo *triangulos;
    int num_triangles;
    mlist *mptr;
    struct triobj *next;
    unsigned char *rgb;
} triobj;


extern triobj *first_cam_ptr; //first_cam_pointer

// información de textura
extern int load_ppm(char *file, unsigned char **bufferptr, int *dimxptr, int *dimyptr);
extern unsigned char *buff;
extern int dimx, dimy, dim;

extern int indexx; // me dice que triangulo del objeto que estoy dibujando
extern triobj *first_object_pointer;   //?
extern triobj *sel_ptr; //?
extern int all;
extern int lines;
extern int objects;
extern int perspective;
extern int cam_val; // si estamos aplicando transformaciones a la camara o  al objeto
extern char trfm;
extern int local_trfm;
extern int view_mode; // 0 es la camara, 1 es el objeto
extern int analisis; // if 1, then we are in analisis mode, if 0, we are in flight mode
extern int backface;

//matrices globales
extern double modelview_matrix[16];
extern double projection_matrix[16];

extern double l, r, b, t, n, f;
// punto E
extern punto at;
extern vector look_dir, up, right;
extern double dist_to_obj;
extern char filename[100];

void do_matrix();
void print_matrix(char *str);
void print_matrix_from_param(double *m);
void dibujar_normal(float x, float y, float z, float nx, float ny, float nz);
void mxp(punto *pptr, double m[16], punto p);
void mPxP(punto *pptr, double m[16], punto p);
void mxv(vector *vptr, double m[16], vector v);
void mxm(double *m1, double *m2, double *res);
double dot_product(vector v1, vector v2);
void cross_product(vector v1, vector v2, vector *result);
void ordenar(punto *p1, punto *p2, punto *p3, punto **superior, punto **intermediate, punto **inferior);
float distancia(punto p1, punto p2);
float area_triangulo(float lado1, float lado2, float lado3);
void calcular_alturas(float lado1, float lado2, float lado3, float *h1, float *h2, float *h3);
void intercambiar_puntos(punto *p1, punto *p2);
void normalize(vector *v);
int tipo_triangulo(punto p1, punto p2, punto p3, float *lado1, float *lado2, float *lado3);
void calc_obj_centre();

#endif