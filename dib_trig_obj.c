//	Program developed by
//
//	Informatika Fakultatea
//	Euskal Herriko Unibertsitatea
//	http://www.ehu.eus/if
//
// to compile it: gcc *.c -lGL -lGLU -lglut -lm -o ejecutable
// to execute it: ./ejecutable
// to do both at the same time: gcc *.c -lGL -lGLU -lglut -lm -o ejecutable; ./ejecutable
//
//
//

#include <GL/glut.h>
#include <stdio.h>
#include <string.h>
#include <math.h>
#include "cargar-triangulo.h"

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

triobj *first_cam_ptr; //first_cam_pointer

// información de textura
extern int load_ppm(char *file, unsigned char **bufferptr, int *dimxptr, int *dimyptr);
unsigned char *buff;
int dimx, dimy, dim;

int indexx; // me dice que triangulo del objeto que estoy dibujando
triobj *first_object_pointer;   //?
triobj *sel_ptr; //?
int all;
int lines;
int objects;
int perspective;
int cam_val; // si estamos aplicando transformaciones a la camara o  al objeto
char trfm;
int local_trfm;
int view_mode; // 0 es la camara, 1 es el objeto
int analisis; // if 1, then we are in analisis mode, if 0, we are in flight mode
int backface;

//matrices globales
double modelview_matrix[16];
double projection_matrix[16];

double l, r, b, t, n, f;
// punto E
punto at;
vector look_dir, up, right;
double dist_to_obj;
char filename[100];

// TODO
// debe devolver el pointer correspondiente a las coordenadas u y v
unsigned char *color_textura(float u, float v)
{
    int displacement, x, y;
    char *lag;
    // printf("texturan...%x\n",buffer);

    x = ((int)(u * (float)(dimx - 1))) % dimx;       // me aseguro siempre que estos valores estén entre 0 y dimx -1
    y = ((int)((1 - v) * (float)(dimy - 1))) % dimy; // formato ppm tiene el (0,0) arriba a la izquierda
    displacement = y * dimx + x;

    lag = (unsigned char *)buff;
    return (lag + 3 * displacement);
}

void dibujar_linea_z(float linea, float c1x, float c1z, float c1u, float c1v, float c2x, float c2z, float c2u, float c2v)
{
    float x_coord, z_coord;
    float cx, cz, cu, cv;
    float u, v;
    int i, point_count;
    unsigned char r, g, b;
    unsigned char *colorv;

    // TODO x balioak -1 eta 1 artekoak direla ziurtatu eta ondorioz z, u i = ord(str1[index]) if eta v egokitu
    // x balioak -1 eta 1 artekoa izan behar du, ezkerretik eskuinera 2 unitateko trfm dauka x balioak, eta
    // dim adina pixel behar dira ezkerretik eskuinera joateko, beraz, 2-ko trfmri "dim" pixel dagozkio.
    // ondorioz c1x-tik c2x-ra doan tartean behar ditudan pixelak = (c2x-c1x)dim/2 pixel behar ditut
    // ondorioz, Z-ren u-ren eta v-ren trfmn zenbaki hori erabili behar dut.
    // para un cambio de -1 a 1 en x hay que dibujar "dim" pixels, para un cambio que va de c1x a c2x, cuantos?
    point_count = (c2x - c1x) * (float)dim / 2.0;

    if (point_count > 0)
    {
        // pixel batetik bestera dagoen distantzia x koordenatuan: cx = (c2x-c1x)/(float)point_count;
        // diustancia en x de un pixel al siguiente: cx
        cx = (c2x - c1x) / (float)point_count;
        // TODO kalkulatu cz, cu eta cv
        //      calcula los valores adecuados
        cz = 0;
        cu = 0;
        cv = 0;
    }
    else
    {
        cx = 0;
        cz = 0;
        cu = 0;
        cv = 0;
    }

    glBegin(GL_POINTS);
    for (i = 0, x_coord = c1x, z_coord = c1z, u = c1u, v = c1v;
         i < point_count /*x_coord <= c2x*/;
         i++, x_coord += cx /* 500 puntu -1 eta 1 artean */)
    {
        colorv = color_textura(u, v); // si esta función es correcta se ve la foto en la ventana
        r = colorv[0];
        g = colorv[1];
        b = colorv[2];
        glColor3ub(r, g, b);
        glVertex3f(x_coord, linea, z_coord);
        z_coord += cz;
        u += cu;
        v += cv;
    }
    glEnd();
}

void dibujar_normal(float x, float y, float z, float nx, float ny, float nz)
{
    glBegin(GL_LINES);
    glVertex3f(x, y, z);
    glVertex3f(x + nx, y + ny, z + nz);
    glEnd();
} 

// TODO: Revisar por si está bien
void do_matrix ()
{
    if (sel_ptr->mptr == NULL)
    {
        sel_ptr->mptr = (mlist *)malloc(sizeof(mlist));
        if (sel_ptr->mptr == NULL)
        {
            printf("Error alocando mptr\n");
            exit(1);
        }
        glLoadIdentity();
        glGetDoublev(GL_MODELVIEW_MATRIX, sel_ptr->mptr->m);
    }
    else
    {
        glLoadMatrixd(sel_ptr->mptr->m);
    }


}

void print_matrix(char *str)
{
    int i;

    printf("%s\n", str);
    for (i = 0; i < 4; i++)
        printf("%lf, %lf, %lf, %lf\n", sel_ptr->mptr->m[i * 4], sel_ptr->mptr->m[i * 4 + 1], sel_ptr->mptr->m[i * 4 + 2],
               sel_ptr->mptr->m[i * 4 + 3]);
}

void print_matrix_from_param(double *m)
{
    int i;

    for (i = 0; i < 4; i++)
        printf("%lf, %lf, %lf, %lf\n", m[i * 4], m[i * 4 + 1], m[i * 4 + 2], m[i * 4 + 3]);
    printf("\n");
}

void mxp(punto *pptr, double m[16], punto p)
{
    // print_matrix("objektuaren matrix\n");
    // printf("puntua = %lf,%lf,%lf\n ",p.x,p.y,p.z);
    pptr->x = m[0] * p.x + m[1] * p.y + m[2] * p.z + m[3];
    pptr->y = m[4] * p.x + m[5] * p.y + m[6] * p.z + m[7];
    pptr->z = m[8] * p.x + m[9] * p.y + m[10] * p.z + m[11];

    pptr->u = p.u;
    pptr->v = p.v;
}

//matriz de proyeccion por punto, teninedo en cuenta ahora que la 4ta componente de la matriz de proyeccion puede ser 0
void mPxP(punto *pptr, double m[16], punto p)
{
    double w;
    pptr->x = m[0] * p.x + m[1] * p.y + m[2] * p.z + m[3];
    pptr->y = m[4] * p.x + m[5] * p.y + m[6] * p.z + m[7];
    pptr->z = m[8] * p.x + m[9] * p.y + m[10] * p.z + m[11];

    w = m[12] * p.x + m[13] * p.y + m[14] * p.z + m[15];

    pptr->x /= w;
    pptr->y /= w;
    pptr->z /= w;

    pptr->u = p.u;
    pptr->v = p.v;
}

void mxv(vector *vptr, double m[16], vector v)
{
    vptr->x = m[0] * v.x + m[1] * v.y + m[2] * v.z + m[3]*0;
    vptr->y = m[4] * v.x + m[5] * v.y + m[6] * v.z + m[7]*0;
    vptr->z = m[8] * v.x + m[9] * v.y + m[10] * v.z + m[11]*0;
}

void mxm(double *m1, double *m2, double *result)
{
    for (int i = 0; i < 4; i++)
    {
        for (int j = 0; j < 4; j++)
        {
            result[i*4 + j] = 0;
            for (int k = 0; k < 4; k++)
            {
                result[i*4 + j] += m1[i*4 + k] * m2[k*4 + j];
            }
        }
    }
}

double dot_product(vector v1, vector v2)
{
    return v1.x * v2.x + v1.y * v2.y + v1.z * v2.z;
}

// el vector c se obtiene de la multiplicación de los vectores a ^ b
void cross_product (vector v1, vector v2, vector *result)
{
    result->x = v1.y * v2.z - v1.z * v2.y;
    result->y = v1.z * v2.x - v1.x * v2.z;
    result->z = v1.x * v2.y - v1.y * v2.x;
}

// inicializacion de la camara
void calc_cam_matrix(triobj *first_cam_ptr){

    if (first_cam_ptr->mptr == NULL) {
        first_cam_ptr->mptr = (mlist *)malloc(sizeof(mlist));
        if (first_cam_ptr->mptr == NULL) {
            fprintf(stderr, "Error: No se pudo asignar memoria para la matriz de la cámara\n");
            exit(1);
        }
    }
    
    // Inicializa la matriz de la cámara
    double *m = first_cam_ptr->mptr->m;

    first_cam_ptr->mptr->m[0] = 1.0;
    first_cam_ptr->mptr->m[4] = 0.0;
    first_cam_ptr->mptr->m[8] = 0.0;
    first_cam_ptr->mptr->m[12] = 0.0;

    first_cam_ptr->mptr->m[1] = 0.0;
    first_cam_ptr->mptr->m[5] = 1.0;
    first_cam_ptr->mptr->m[9] = 0.0;
    first_cam_ptr->mptr->m[13] = 0.0;
    
    first_cam_ptr->mptr->m[2] = 0.0;
    first_cam_ptr->mptr->m[6] = 0.0;
    first_cam_ptr->mptr->m[10] = 1.0;
    first_cam_ptr->mptr->m[14] = 0.0;

    first_cam_ptr->mptr->m[3] = 0.0;
    first_cam_ptr->mptr->m[7] = 0.0;
    first_cam_ptr->mptr->m[11] = 2.0;
    first_cam_ptr->mptr->m[15] = 1.0;
}


void calc_projection_matrix()
{
    if (perspective == 1)
    {
        // (P)erspectiva
        // n es positiva
        // tengo que invertir el signo de la linea de la z
        l = -0.1, r = 0.1, b = -0.1, t = 0.1, n = 0.1, f = 100.0; 

        projection_matrix[0] = 2 * n / (r - l);
        projection_matrix[1] = 0;
        projection_matrix[2] = (r + l) / (r - l);
        projection_matrix[3] = 0;

        projection_matrix[4] = 0;
        projection_matrix[5] = 2 * n / (t - b);
        projection_matrix[6] = (t + b) / (t - b);
        projection_matrix[7] = 0;

        projection_matrix[8] = 0;
        projection_matrix[9] = 0;
        projection_matrix[10] = -(f + n) / (n - f);
        projection_matrix[11] = -(2 * f * n) / (n - f);

        projection_matrix[12] = 0;
        projection_matrix[13] = 0;
        projection_matrix[14] = -1;
        projection_matrix[15] = 0;
    }
    else
    {
        // (p)aralelo

        l = -1.0, r = 1.0, b = -1.0, t = 1.0, n = 0.1, f = -100.0;

        projection_matrix[0] = 2 / (r - l);
        projection_matrix[1] = 0;
        projection_matrix[2] = 0;
        projection_matrix[3] = -(r + l) / (r - l);

        projection_matrix[4] = 0;
        projection_matrix[5] = 2 / (t - b);
        projection_matrix[6] = 0;
        projection_matrix[7] = -(t + b) / (t - b);

        projection_matrix[8] = 0;
        projection_matrix[9] = 0;
        projection_matrix[10] = 2 / (n - f);
        projection_matrix[11] = -(n + f) / (n - f);

        projection_matrix[12] = 0;
        projection_matrix[13] = 0;
        projection_matrix[14] = 0;
        projection_matrix[15] = 1;
    }
    // printf("La matriz de proyección es la siguiente: \n");
    // for (int i = 0; i < 4; i++) {
    //     for (int j = 0; j < 4; j++) {
    //         printf("%lf ", projection_matrix[i * 4 + j]);
    //     }
    //     printf("\n");
    // }

}

// TODO (transform): obtener la matriz que pasa del sistema de referencia del objeto al sistema de referencia del mundo
void obtener_CSR_partiendo_de_M(GLdouble *M, GLdouble *MCSR)
{
    MCSR[0] = M[0];
    MCSR[4] = M[1];
    MCSR[8] = M[2];
    MCSR[12] = 0;

    MCSR[1] = M[4];
    MCSR[5] = M[5];
    MCSR[9] = M[6];
    MCSR[13] = 0;

    MCSR[2] = M[8];
    MCSR[6] = M[9];
    MCSR[10] = M[10];
    MCSR[14] = 0;

    MCSR[3] = -(M[0] * M[3] + M[4] * M[7] + M[8] * M[11]);
    MCSR[7] = -(M[1] * M[3] + M[5] * M[7] + M[9] * M[11]);
    MCSR[11] = -(M[2] * M[3] + M[6] * M[7] + M[10] * M[11]);
    MCSR[15] = 1;
}



void ordenar(punto *p1, punto *p2, punto *p3, punto **superior, punto **intermediate, punto **inferior)
{
    // los casos especiales me dan igual. este codigo cubre bien estas posibilidades.
    if (p1->y > p2->y)
    {
        *superior = p1;
        *intermediate = p2;
    }
    else
    {
        *superior = p2;
        *intermediate = p1;
    }

    if (p3->y > (*superior)->y)
    {
        *inferior = *intermediate;
        *intermediate = *superior;
        *superior = p3;
    }
    else if (p3->y > (*intermediate)->y)
    {
        *inferior = *intermediate;
        *intermediate = p3;
    }
    else
    {
        *inferior = p3;
    }
}

float distancia(punto p1, punto p2)
{
    return sqrt((p1.x - p2.x) * (p1.x - p2.x) + (p1.y - p2.y) * (p1.y - p2.y) + (p1.z - p2.z) * (p1.z - p2.z));
}

float area_triangulo(float lado1, float lado2, float lado3)
{
    float s = (lado1 + lado2 + lado3) / 2;
    return sqrt(s * (s - lado1) * (s - lado2) * (s - lado3));
}

void calcular_alturas(float lado1, float lado2, float lado3, float *h1, float *h2, float *h3)
{
    float A = area_triangulo(lado1, lado2, lado3);

    *h1 = (2 * A) / lado1; // altura desde p1 hacia el lado opuesto (p2, p3)
    *h2 = (2 * A) / lado2; // altura desde p2 hacia el lado opuesto (p1, p3)
    *h3 = (2 * A) / lado3; // altura desde p3 hacia el lado opuesto (p1, p2)
}

void intercambiar_puntos(punto *p1, punto *p2)
{
    punto temp = *p1;
    *p1 = *p2;
    *p2 = temp;
}

void normalize(vector *v) {
    float norma = sqrt(v->x * v->x + v->y * v->y + v->z * v->z);
    v->x /= norma;
    v->y /= norma;
    v->z /= norma;
}

int tipo_triangulo(punto p1, punto p2, punto p3, float *lado1, float *lado2, float *lado3)
{
    /**
    * Devuelve el tipo de triangulo: acutángulo, rectángulo u obtusángulo
    */

    // Calcular las longitudes de los lados del triángulo
    *lado1 = distancia(p2, p3); // lado opuesto a p1
    *lado2 = distancia(p1, p3); // lado opuesto a p2
    *lado3 = distancia(p1, p2); // lado opuesto a p3

    // Calcular los cuadrados de las longitudes de los lados
    float a2 = *lado1 * *lado1;
    float b2 = *lado2 * *lado2;
    float c2 = *lado3 * *lado3;

    // Calcular el coseno de los ángulos
    float cos_a = (b2 + c2 - a2) / (2 * *lado2 * *lado3);
    float cos_b = (a2 + c2 - b2) / (2 * *lado1 * *lado3);
    float cos_c = (a2 + b2 - c2) / (2 * *lado1 * *lado2);

    // printf("cos_a: %f, cos_b: %f, cos_c: %f\n", cos_a, cos_b, cos_c);

    // Calcular el tipo de triángulo
    if (cos_a > 0 && cos_b > 0 && cos_c > 0)
    {
        return 0; // acutángulo
    }
    else if (cos_a == 0 || cos_b == 0 || cos_c == 0)
    {
        return 1; // rectángulo
    }
    else
    {
        return 2; // obtusángulo
    }
}

void dibujar_triangulo(triobj *optr, int ti)
{
    if (ti >= optr->num_triangles)
        return;

    triangulo *tptr = optr->triangulos + ti;
    
    punto p;          // punto actual
    punto p1, p2, p3; // puntos del triangulo
    // vector n;         // normal del triangulo

    punto p1aux, p2aux, p3aux;
    vector nAux, dir_cam, pos_cam;
    float c;    

    float alpha, beta, gamma;
    unsigned char *colorv;
    float incremento_alpha, incremento_beta;
    float lado1, lado2, lado3; // Longitudes de los lados del triángulo

    int puntos_totales = 0;

    mxp(&p1aux, modelview_matrix, tptr->p1);
    mxp(&p2aux, modelview_matrix, tptr->p2);
    mxp(&p3aux, modelview_matrix, tptr->p3);

    // Verificar si alguno de los puntos tiene z <= -0.1
    if (p1aux.z >= -n || p2aux.z >= -n || p3aux.z >= -n) {
        return; // No dibujamos el triángulo
    }

    mxv(&nAux, modelview_matrix, tptr->N);
    // printf("Normal: %f, %f, %f\n", nAux.x, nAux.y, nAux.z);

    // tengo que decidir si dibujar o no el triangulo dependiendo de si el punto está por detrás de mi plano de proyección
    if (!perspective){
        if (nAux.z < 0){    
            if (!backface)
                return;
            glColor3ub(0, 0, 0);
        }else{
            glColor3ub(255, 255, 255);
        }
    }else{
        
        c = -nAux.x * p1aux.x - nAux.y * p1aux.y - nAux.z * p1aux.z; 

        if (c < 0)
        {
            if (!backface)
                return;
            glColor3ub(0, 0, 0);
        }
        else
        {
            glColor3ub(255, 255, 255);
        }
    }

    // como la 4ta componenete de la matriz de proyeccion en perspectiva puede ser 0, ya no podemos usar mxp
    // creamos una nueva función, mPxP

    mPxP(&p1, projection_matrix, p1aux);
    mPxP(&p2, projection_matrix, p2aux);
    mPxP(&p3, projection_matrix, p3aux);
    
    if (lines == 1)
    {
        glBegin(GL_POLYGON);
        glVertex3f(p1.x, p1.y, p1.z);
        glVertex3f(p2.x, p2.y, p2.z);
        glVertex3f(p3.x, p3.y, p3.z);
        glEnd();

        // glBegin(GL_LINES);
        // glVertex3f(p1.x, p1.y, p1.z);
        // glVertex3f(p1.x + n.x, p1.y + n.y, p1.z + n.z);        
        // glEnd();

        return;
    }

    // Tipo de triangulo
    int tipo; // 0: acutángulo, 1: rectángulo, 2: obtusángulo

    // Calcular el tipo de triángulo
    tipo = tipo_triangulo(p1, p2, p3, &lado1, &lado2, &lado3);
    // printf("Tipo de triángulo: %d\n", tipo);

    float h1, h2, h3; // Alturas del triángulo

    // Si es obtusángulo, el incremento depende del lado mayor
    if (tipo == 2)
    {
        h1 = lado2;
        h2 = lado3;
    }
    else
    {
        // Calcular las alturas del triángulo
        calcular_alturas(lado1, lado2, lado3, &h1, &h2, &h3);

        float temp;
        // Ordenar los puntos de mayor a menor altura
        if (h1 < h2)
        {
            intercambiar_puntos(&p1, &p2);
            temp = h1;
            h1 = h2;
            h2 = temp;
        }
        if (h1 < h3)
        {
            intercambiar_puntos(&p1, &p3);
            temp = h1;
            h1 = h3;
            h3 = temp;
        }
        if (h2 < h3)
        {
            intercambiar_puntos(&p2, &p3);
            temp = h2;
            h2 = h3;
            h3 = temp;
        }
        // printf("h1: %f, h2: %f, h3: %f\n", h1, h2, h3);
    }
    
    float pixeldist = 2.0 / (dim - 1);

    incremento_alpha = pixeldist / h1;
    incremento_beta = pixeldist / h2;

    int puntos_dibujados = 0;

    glBegin(GL_POINTS);
    for (alpha = 0.0; alpha <= 1.0; alpha += incremento_alpha)
    {
        for (beta = 0.0; beta <= 1.0 - alpha; beta += incremento_beta)
        {
            puntos_dibujados++;

            gamma = 1.0 - alpha - beta;

            p.x = p1.x * alpha + p2.x * beta + p3.x * gamma;
            p.y = p1.y * alpha + p2.y * beta + p3.y * gamma;
            p.z = p1.z * alpha + p2.z * beta + p3.z * gamma;

            p.u = p1.u * alpha + p2.u * beta + p3.u * gamma;
            p.v = p1.v * alpha + p2.v * beta + p3.v * gamma;

            colorv = color_textura(p.u, p.v);
            glColor3ub(colorv[0], colorv[1], colorv[2]);

            // glColor3ub(255*alpha, 255*beta, 255*gamma);

            glVertex3f(p.x, p.y, p.z);
        }
    }
    glEnd();

    // printf("Puntos dibujados: %d\n", puntos_dibujados);
}



static void draw(void)
{
    float u, v;
    int i, j;
    double M_csr[16];
    triobj *auxptr;

    if (first_cam_ptr == NULL || first_cam_ptr->mptr == NULL) {
        fprintf(stderr, "Error: La cámara no está inicializada correctamente a la hora de dibujar.\n");
        return;
    }

    if (view_mode == 0) {
        // Vista de la cámara
        obtener_CSR_partiendo_de_M(first_cam_ptr->mptr->m, M_csr);
    } else {
        // Vista del objeto seleccionado
        if (sel_ptr && sel_ptr->mptr) {
            obtener_CSR_partiendo_de_M(sel_ptr->mptr->m, M_csr);
        } else {
            // Si no hay objeto seleccionado, usar la vista de la cámara
            obtener_CSR_partiendo_de_M(first_cam_ptr->mptr->m, M_csr);
        }
    }

    calc_projection_matrix();

    // no se puede dibujar sin objetos
    if (first_object_pointer == 0)
        return;

    // clear viewport...
    if (objects == 1)
        glClear(GL_DEPTH_BUFFER_BIT | GL_COLOR_BUFFER_BIT);
    else
    {
        if (all == 0)
            glClear(GL_DEPTH_BUFFER_BIT | GL_COLOR_BUFFER_BIT);
    }
    
    if (objects == 1) // si tengo que dibujar objetos completos
    {
        if (all == 1) // si sin todos los objetos
        {
            for (auxptr = first_object_pointer; auxptr != 0; auxptr = auxptr->next) // con un auxiliar recorro toda la lista
            {
                if (view_mode == 0 || auxptr != sel_ptr) {
                    // Vista de cámara o objeto no seleccionado
                    mxm(M_csr, auxptr->mptr->m, modelview_matrix);
                } else {
                    // Vista del objeto seleccionado
                    mxm(M_csr, sel_ptr->mptr->m, modelview_matrix);
                }
                for (i = 0; i < auxptr->num_triangles; i++)
                {
                    dibujar_triangulo(auxptr, i); // dibujo el auxiliar
                }
            }
        }
        else
        {
            mxm(M_csr, sel_ptr->mptr->m, modelview_matrix);
            for (i = 0; i < sel_ptr->num_triangles; i++)
            {
                mxm(M_csr, sel_ptr->mptr->m, modelview_matrix);
                dibujar_triangulo(sel_ptr, i);
            }
        }
    }
    else
    {
        mxm(M_csr, sel_ptr->mptr->m, modelview_matrix); 
        dibujar_triangulo(sel_ptr, indexx);
    }

    glFlush();
    // print_matrix("matriz de transformación");

}

void read_from_file(char *fitx)
{
    int i, retval;
    triobj *optr;
    vector v1, v2;
    double norma;
    triangulo *tri;

    printf("%s Cogiendo datos de fichero \n", fitx);
    optr = (triobj *)malloc(sizeof(triobj));
    optr->rgb = 0;
    retval = cargar_triangulos_color(fitx, &(optr->num_triangles), &(optr->triangulos), &(optr->rgb));
    // TODO (transform...)
    // int cargar_triangulos_color(char *filename, int *hkopptr, triangulo **nextptr,unsigned char **rgbptr);
    // retval = cargar_triangulos_color(...)
    if ((retval != 15) && (retval != 9))
    {
        printf("%s Error al procesar datos del fichero\n", filename);
        free(optr);
    }
    else
    {

        // calcular la normal de cada triangulo
        for (i = 0; i < optr->num_triangles; i++)
        {
            tri = &optr->triangulos[i];
            //calcular la normal del triangulo
            v1.x = tri->p2.x - tri->p1.x;
            v1.y = tri->p2.y - tri->p1.y;
            v1.z = tri->p2.z - tri->p1.z;

            v2.x = tri->p3.x - tri->p1.x;
            v2.y = tri->p3.y - tri->p1.y;
            v2.z = tri->p3.z - tri->p1.z;

            //calcular la normal
            tri->N.x = v1.y * v2.z - v1.z * v2.y;
            tri->N.y = v1.z * v2.x - v1.x * v2.z;
            tri->N.z = v1.x * v2.y - v1.y * v2.x;

            //normalizar la normal
            norma = sqrt(tri->N.x * tri->N.x + tri->N.y * tri->N.y + tri->N.z * tri->N.z);
            tri->N.x /= norma;
            tri->N.y /= norma;
            tri->N.z /= norma;
        }

        optr->mptr = (mlist *)malloc(sizeof(mlist));
        for (i = 0; i < 16; i++)
            optr->mptr->m[i] = 0;
        optr->mptr->m[0] = 1.0;
        optr->mptr->m[5] = 1.0;
        optr->mptr->m[10] = 1.0;
        optr->mptr->m[15] = 1.0;
        optr->mptr->next = 0;
        // printf("objektu zerrendara doa informazioa...\n");
        optr->next = first_object_pointer;
        first_object_pointer = optr;
        sel_ptr = optr;
    }
    printf("datuak irakurrita\n");
}

void calc_obj_centre() {
    // if (!sel_ptr) return;
    
    // Reiniciar at a 0 antes de calcular el nuevo centro
    at.x = 0.0;
    at.y = 0.0;
    at.z = 0.0;
    
    punto p_transformed;
    
    for (int i = 0; i < sel_ptr->num_triangles; i++) {
        // Transformar cada vértice usando la matriz del objeto
        mxp(&p_transformed, sel_ptr->mptr->m, sel_ptr->triangulos[i].p1);
        at.x += p_transformed.x;
        at.y += p_transformed.y;
        at.z += p_transformed.z;
        
        mxp(&p_transformed, sel_ptr->mptr->m, sel_ptr->triangulos[i].p2);
        at.x += p_transformed.x;
        at.y += p_transformed.y;
        at.z += p_transformed.z;
        
        mxp(&p_transformed, sel_ptr->mptr->m, sel_ptr->triangulos[i].p3);
        at.x += p_transformed.x;
        at.y += p_transformed.y;
        at.z += p_transformed.z;
    }

    // Dividir por el número total de vértices (3 por triángulo)
    at.x /= (sel_ptr->num_triangles * 3);
    at.y /= (sel_ptr->num_triangles * 3);
    at.z /= (sel_ptr->num_triangles * 3);
}

void orientate_cam_to_obj() {
    // Actualizar la matriz de la cámara para que mire al objeto seleccionado
    first_cam_ptr->mptr->m[0] = right.x;
    first_cam_ptr->mptr->m[4] = right.y;
    first_cam_ptr->mptr->m[8] = right.z;

    first_cam_ptr->mptr->m[1] = up.x;
    first_cam_ptr->mptr->m[5] = up.y;
    first_cam_ptr->mptr->m[9] = up.z;
    
    first_cam_ptr->mptr->m[2] = -look_dir.x;
    first_cam_ptr->mptr->m[6] = -look_dir.y;
    first_cam_ptr->mptr->m[10] = -look_dir.z;
}


void update_cam_params() {

    if (analisis) {
        // Calcular el centro del objeto seleccionado
        calc_obj_centre();

        // En modo análisis, calcular look_dir hacia el centro del objeto
        look_dir.x = at.x - first_cam_ptr->mptr->m[3];
        look_dir.y = at.y - first_cam_ptr->mptr->m[7];
        look_dir.z = at.z - first_cam_ptr->mptr->m[11];
        normalize(&look_dir);

        // Mantener up constante
        up.x = 0.0;
        up.y = 1.0;
        up.z = 0.0;

        // Calcular right
        cross_product(up, look_dir, &right);
        normalize(&right);

        // Recalcular up para asegurar ortogonalidad
        cross_product(look_dir, right, &up);
        normalize(&up);

        orientate_cam_to_obj();

        print_matrix_from_param(first_cam_ptr->mptr->m);

    } else {
        // En modo vuelo, obtener look_dir de la matriz de la cámara
        look_dir.x = -first_cam_ptr->mptr->m[2];
        look_dir.y = -first_cam_ptr->mptr->m[6];
        look_dir.z = -first_cam_ptr->mptr->m[10];
        normalize(&look_dir);
    }
}


void init_camera() {
    first_cam_ptr = (triobj *)malloc(sizeof(triobj));
    if (first_cam_ptr == NULL) {
        fprintf(stderr, "Error: No se pudo asignar memoria para la cámara\n");
        exit(1);
    }
    first_cam_ptr->num_triangles = 0;
    first_cam_ptr->triangulos = NULL;
    first_cam_ptr->next = NULL;
    first_cam_ptr->rgb = NULL;
    
    calc_cam_matrix(first_cam_ptr);

    // Inicializar E y at
    update_cam_params();
}

void scale(int sign)
{
    // do_matrix();

    // printf("sign: %d\n", sign);
    double factor = sign * 0.05;
    mlist *target = (cam_val == 1) ? first_cam_ptr->mptr : sel_ptr->mptr;
    
    target->m[0] += factor;
    target->m[5] += factor;
    target->m[10] += factor;
}

void translate(double x, double y, double z) {
    double factor = 0.05;

    if (cam_val == 1) {
        vector direction;
        if (analisis == 0) { // Modo vuelo
            // Obtener la dirección de vista actual de la cámara
            direction.x = -first_cam_ptr->mptr->m[2];
            direction.y = -first_cam_ptr->mptr->m[6];
            direction.z = -first_cam_ptr->mptr->m[10];
            normalize(&direction);

            // Actualizar la posición usando la dirección actual
            first_cam_ptr->mptr->m[3] += direction.x * z * factor;
            first_cam_ptr->mptr->m[7] += direction.y * z * factor;
            first_cam_ptr->mptr->m[11] += direction.z * z * factor;

        } else {  // modo analisis
            vector new_pos;
            // Calcular nueva posición propuesta
            new_pos.x = first_cam_ptr->mptr->m[3] + look_dir.x * z * 0.05;
            new_pos.y = first_cam_ptr->mptr->m[7] + look_dir.y * z * 0.05;
            new_pos.z = first_cam_ptr->mptr->m[11] + look_dir.z * z * 0.05;

            // Calcular distancia al objeto usando at (centro del objeto)
            vector to_obj;
            to_obj.x = at.x - new_pos.x;
            to_obj.y = at.y - new_pos.y;
            to_obj.z = at.z - new_pos.z;
            double new_distance = sqrt(to_obj.x * to_obj.x + 
                                    to_obj.y * to_obj.y + 
                                    to_obj.z * to_obj.z);

            // Solo actualizar si la distancia es válida
            double min_distance = 0.2;
            double max_distance = 10.0;
            if (new_distance >= min_distance && new_distance <= max_distance) {
                first_cam_ptr->mptr->m[3] = new_pos.x;
                first_cam_ptr->mptr->m[7] = new_pos.y;
                first_cam_ptr->mptr->m[11] = new_pos.z;
            }
        }
        update_cam_params();
    } else {
        // Para objetos mantener el comportamiento actual
        mlist *target = sel_ptr->mptr;
        target->m[3] += x * factor;
        target->m[7] += y * factor;
        target->m[11] += z * factor;
    }
}

void delete_object(triobj *optr)
{
    if (first_object_pointer != 0)
    {
        triobj *temp = sel_ptr;
        sel_ptr = sel_ptr->next;
        if (sel_ptr == 0)
            sel_ptr = first_object_pointer;
        
        // Eliminar el objeto actual
        if (temp == first_object_pointer)
        {
            first_object_pointer = first_object_pointer->next;
        }
        else
        {
            triobj *prev = first_object_pointer;
            while (prev->next != temp)
            {
                prev = prev->next;
            }
            prev->next = temp->next;
        }

        // Liberar memoria
        free(temp->triangulos);

        // Liberar matrices
        mlist *mtemp = temp->mptr;
        while (mtemp != 0)
        {
            mlist *mnext = mtemp->next;
            free(mtemp);
            mtemp = mnext;
        }

        free(temp);

        indexx = 0;

        //gl flush
        if (first_object_pointer == 0)
        {
            glClear(GL_DEPTH_BUFFER_BIT | GL_COLOR_BUFFER_BIT);
            glFlush();
        }   
    }
}

void rotate_around_object(char axis, int sign) {
    if (!sel_ptr) return;
    
    double angle = sign * 0.05;
    double c = cos(angle);
    double s = sin(angle);
    
    // Matriz de traslación al punto de atención
    double M1[16] = {
        1, 0, 0, -at.x,
        0, 1, 0, -at.y,
        0, 0, 1, -at.z,
        0, 0, 0, 1
    };
    
    // Matriz de traslación inversa
    double M2[16] = {
        1, 0, 0, at.x,
        0, 1, 0, at.y,
        0, 0, 1, at.z,
        0, 0, 0, 1
    };
    
    // Matriz de rotación
    double R[16] = {0};
    for (int i = 0; i < 16; i++) 
        R[i] = (i % 5 == 0) ? 1.0 : 0.0;
    
    if (axis == 'x') {
        // Rotación alrededor del eje X
        R[5] = c;
        R[6] = -s;
        R[9] = s;
        R[10] = c;
    } else { // axis == 'y'
        // Rotación alrededor del eje Y
        R[0] = c;
        R[2] = s;
        R[8] = -s;
        R[10] = c;
    }
    
    // Calcular Ml = M2 * R * M1
    double temp[16], result[16];
    mxm(R, M1, temp);       // temp = R * M1
    mxm(M2, temp, result);  // result = M2 * temp
    
    // M_csr' = Ml * M_csr
    mxm(result, first_cam_ptr->mptr->m, temp);
    
    // Actualizar la matriz de la cámara
    memcpy(first_cam_ptr->mptr->m, temp, sizeof(double) * 16);
    
    // Actualizar los parámetros de la cámara
    // update_cam_params();
}




void rotate(double x, double y, double z)
{
    double rotation_matrix[16] = {0};
    double factor = 0.05;
    double angle = factor * (x != 0 ? x : (y != 0 ? y : z));
    double c = cos(angle);
    double s = sin(angle);

    // Seleccionar la matriz a transformar, cámara u objeto
    mlist *target = (cam_val == 1) ? first_cam_ptr->mptr : sel_ptr->mptr;

    // Matriz identidad
    rotation_matrix[0] = rotation_matrix[5] = rotation_matrix[10] = rotation_matrix[15] = 1.0;

    if (x != 0)
    {
        // Rotación alrededor del eje X
        rotation_matrix[5] = c;
        rotation_matrix[6] = -s;
        rotation_matrix[9] = s;
        rotation_matrix[10] = c;
    }
    else if (y != 0)
    {
        // Rotación alrededor del eje Y
        rotation_matrix[0] = c;
        rotation_matrix[2] = s;
        rotation_matrix[8] = -s;
        rotation_matrix[10] = c;
    }
    else if (z != 0)
    {
        // Rotación alrededor del eje Z
        rotation_matrix[0] = c;
        rotation_matrix[1] = -s;
        rotation_matrix[4] = s;
        rotation_matrix[5] = c;
    }

    double result[16] = {0};

    if (local_trfm == 1) {
        mxm(target->m, rotation_matrix, result);
    } else {
        mxm(rotation_matrix, target->m, result);
    }
    
    for (int i = 0; i < 16; i++) {
        target->m[i] = result[i];
    }

    if (cam_val == 1) {
        update_cam_params();
    }
}

void x_trfm(int sign)
{
    if (cam_val == 1 && analisis == 0) { // modo vuelo
        rotate(sign, 0, 0);
    } else
    {
        if (trfm == 't')
            translate(sign, 0, 0);
        else
            rotate(sign, 0, 0);
    }
}

void y_trfm(int sign)
{
    if (cam_val == 1 && analisis == 0) { // modo vuelo
        rotate(0, sign, 0);
    } else
    {
        if (trfm == 't')
            translate(0, sign, 0);
        else
            rotate(0, sign, 0);
    }
}

void z_trfm(int sign)
{
    if (trfm == 't')
        translate(0, 0, sign);
    else
        rotate(0, 0, sign);

}

static void keyboard(unsigned char key, int x, int y)
{
    int retval;
    int i;
    char selected;
    FILE *obj_file;

    switch (key)
    {
    case 13:            // intro
        if (first_object_pointer != 0) // si no hay objeto que no haga nada
        {
            indexx++; // pero si es el último? hay que controlarlo! // se refiere al indexxesimo triangulo del objeto seleccionado
            if (indexx == sel_ptr->num_triangles)
            {
                indexx = 0;
                // "esto es una chorrada, da igual"
                if ((all == 1) && (objects == 0))
                {
                    glClear(GL_DEPTH_BUFFER_BIT | GL_COLOR_BUFFER_BIT);
                    glFlush();
                }
            }
        }
        break;
    case 'd':
        all = 1 - all;
        break;
    case 'b':
    case 'B':
        backface = 1 - backface;
        printf("Backface culling %s\n", backface ? "activado" : "desactivado");
        break;
    case 'o':
        if (objects == 1)
            objects = 0;
        else
            objects = 1;
        break;
    case 's':
        scale(1);
        break;
    case 'S':
        scale(-1);
        break;
    case 'l':
        if (lines == 1)
            lines = 0;
        else
            lines = 1;
        break;
    case 't':
        trfm = 't';
        printf("=== \n");
        printf("- Traslación activada\n");
        break;
    case 'r':
        trfm = 'r';
        printf("=== \n");
        printf("- Rotación activada\n");
        break;
    case 'G':
    case 'g':
        printf("=== \n");      
        if (cam_val == 1) { // si estamos aplicando transformaciones a la camara
            analisis = 1 - analisis; // cambiamos entre modo analisis y modo vuelo
            printf("La camara se encuentra en modo %s\n", analisis ? "Analisis" : "Vuelo");
            if (analisis){
                update_cam_params();
            }
        } else {            // si estamos aplicando transformaciones a un objeto
            local_trfm = 1 - local_trfm; // cambiamos entre transformaciones locales y globales 
            printf("Transformaciones locales %s\n", local_trfm ? "activadas" : "desactivadas");
        }
        break;

    case 'C':
            //TODO: Si le damos a esta tecla, pasaremos de visualizar lo que ve la camara a lo que ve el objeto seleccionado (entiendo que el objeto seleccionado por el modo analisis?)
        view_mode = 1 - view_mode;
        printf("Visualizamos lo que ve %s\n", view_mode ? "el objeto seleccionado" : "la cámara");
        break;
    case 'c':
        cam_val = 1 - cam_val;
        printf("=== \n");
        printf("Transformaciones aplicadas a: %s\n", cam_val ? "Cámara" : "Objeto");
        if (cam_val == 1) {
            printf("- La camara se encuentra en modo %s\n", analisis ? "Analisis" : "Vuelo");
        }
        printf("- %s activada\n", trfm == 't' ? "Traslación" : "Rotación");
        break;
    case 'p':
    case 'P':
        perspective = 1 - perspective;
        printf("=== \n");
        printf("Visualización de vista en %s\n", perspective ? "Perspectiva" : "Paralelo");
        break;
    case 'x':
        if (cam_val == 0){ // transformaciones en el objeto
            x_trfm(1);
        } else if (cam_val == 1){
            if (analisis == 0){
                y_trfm(1);
            } else {
                rotate_around_object('x', 1);
            }
            update_cam_params();
        }
        break;
    case 'y':
        if (cam_val == 0){ //
            y_trfm(1);
        } else if (cam_val == 1){
            if (analisis == 0){
                x_trfm(1);
            } else {
                selected = 'y';
                rotate_around_object('y', 1);
            }
            update_cam_params();
        }
        break;
    case 'z':
        z_trfm(1);
        break;
    case 'X':
        if (cam_val == 0){
            x_trfm(-1);
        } else if (cam_val == 1){
            if (analisis == 0){
                y_trfm(-1);
            } else {
                rotate_around_object('x', -1);
            }
            update_cam_params();
        }
        break;
    case 'Y':
        if (cam_val == 0){
            y_trfm(-1);
        } else if (cam_val == 1){
            if (analisis == 0){
                x_trfm(-1);
            } else {
                selected = 'y';
                rotate_around_object('y', -1);
            }
            update_cam_params();
        }
        break;
    case 'Z':
        z_trfm(-1);
        break;
    case 'u':
        // undo();
        break;
    case 'f':
        /*Ask for file*/
        printf("=== \n");
        printf("Escribe el nombre del fichero: \n");
        scanf("%s", &(filename[0]));
        read_from_file(filename);
        indexx = 0;
        break;

    case 8:
        if (first_object_pointer != 0)
        {
            delete_object(sel_ptr);
            if (first_object_pointer != 0)
                sel_ptr = first_object_pointer;
            else
                sel_ptr = 0;
        }
        break;
    case 9:             /* <TAB> */
        if (first_object_pointer != 0) // si no hay objeto no hace nada
        {
            sel_ptr = sel_ptr->next;
            
            /*The selection is circular, thus if we move out of the list we go back to the first element*/
            if (sel_ptr == 0)
                sel_ptr = first_object_pointer;

            if (analisis == 1) {
            
            update_cam_params();
            // printf("Se ha dado al tab y cambiamos de objeto. \n");
            // calc_obj_centre();
            // printf("Centro del objeto: %f, %f, %f\n", at.x, at.y, at.z);
                
            //     printf("orientamos la camara al objeto\n");
            //     orientate_cam_to_obj();
            //     printf("camara orientada al objeto\n");

            //     printf("actualizamos los parametros de la camara\n");
            //     update_cam_params();
            //     printf("parametros de la camara actualizados\n");
            }
            indexx = 0; // the selected polygon is the first one
        }
        break;
    case 27: // <ESC>
        exit(0);
        break;
    default:
        printf("%d %c\n", key, key);
    }

    // The screen must be drawn to show the new triangle
    glutPostRedisplay();
}

void newViewPort(int zabal, int garai)
{
    if (zabal < garai)
        dim = zabal;
    else
        dim = garai;
    glViewport(0, 0, dim, dim);
    printf("Nueva cantidad de líneas = %d\n", dim);
}

int main(int argc, char **argv)
{
    int retval;

    // inicializamos la matriz de la cámara
    init_camera();

    printf(" Triangeluak: barneko puntuak eta testura\n Triángulos con puntos internos y textura \n");
    printf("Press <ESC> to finish\n");
    glutInit(&argc, argv);
    glutInitDisplayMode(GLUT_RGB | GLUT_DEPTH);
    dim = 500;
    glutInitWindowSize(dim, dim);
    glutInitWindowPosition(100, 100);
    glutCreateWindow("KBG/GO praktika");

    glutDisplayFunc(draw);
    glutKeyboardFunc(keyboard);
    glutReshapeFunc(newViewPort);
    /* we put the information of the texture in the buffer pointed by buff. The dimensions of the texture are loaded into dimx and dimy */
    retval = load_ppm("testura.ppm", &buff, &dimx, &dimy);
    if (retval != 1)
    {
        printf("No hay archivo de textura (testura.ppm)\n");
        exit(-1);
    }

    glClearColor(0.0f, 0.0f, 0.2f, 0.2f);
    glPolygonMode(GL_FRONT_AND_BACK, GL_LINE);
    glEnable(GL_DEPTH_TEST); // activar el test de profundidad (Z-buffer)
    glDepthFunc(GL_GREATER);
    glClearDepth(0.0);
    all = 1;                 // dibujame todos los triangulos a la vez
    lines = 1;               // dibujame los polignos mediante lineas (0 es que no)
    objects = 1;             // todos los objetos a la vez
    first_object_pointer = 0;
    backface = 0;
    view_mode = 0;          
    perspective = 1;
    cam_val = 0;             // 1 = las transformaciones se aplican a la camara
    sel_ptr = 0;
    trfm = 'r';
    local_trfm = 1;

    if (argc > 1)
        read_from_file(argv[1]);
    else
        read_from_file("abioia-1+1.txt");
        if (sel_ptr != 0) 
            {
            sel_ptr->mptr->m[3] = -1.0;
            sel_ptr = sel_ptr->next;
            }   
        read_from_file("abioia-1+1.txt");
        if (sel_ptr != 0) 
            {
            sel_ptr->mptr->m[3] = 1.0;
            sel_ptr = sel_ptr->next;
            }   
        read_from_file("abioia-1+1.txt");
        if (sel_ptr != 0) 
            { 
            sel_ptr->mptr->m[7] = -0.4;
            sel_ptr->mptr->m[11] = 0.5;
            sel_ptr = sel_ptr->next;
            }   
        read_from_file("abioia-1+1.txt");
        if (sel_ptr != 0) 
            {
            sel_ptr->mptr->m[11] = -0.7;
            sel_ptr = sel_ptr->next;
            }        
        read_from_file("abioia-1+1.txt");
        if (sel_ptr != 0) 
            {
            sel_ptr->mptr->m[7] = 0.7;
            sel_ptr = sel_ptr->next;
            }                        
        read_from_file("abioia-1+1.txt");
        // read_from_file("z-1+1.txt");
        // read_from_file("z-1+1.txt");
                
        printf("Empezamos con el backface culling %s\n", backface ? "desactivado" : "activado");
        printf("Empezamos con la visualización de vista en %s\n", perspective ? "Perspectiva" : "Paralelo");
        printf("Empezamos con las transformaciones aplicadas a %s\n", cam_val ? "Cámara" : "Objeto");
        printf("Empezamos con el modo de visualización de %s\n", view_mode ? "Objeto" : "Cámara");
    glutMainLoop();


    return 0;
}
