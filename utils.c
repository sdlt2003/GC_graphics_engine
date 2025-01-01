
#include "utils.h"
triobj *first_cam_ptr; //first_cam_pointer

int load_ppm(char *file, unsigned char **bufferptr, int *dimxptr, int *dimyptr);
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

double modelview_matrix[16];
double projection_matrix[16];
double l, r, b, t, n, f;

punto at;
vector look_dir, up, right;
double dist_to_obj;
char filename[100];

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

void dibujar_normal(float x, float y, float z, float nx, float ny, float nz)
{
    glBegin(GL_LINES);
    glVertex3f(x, y, z);
    glVertex3f(x + nx, y + ny, z + nz);
    glEnd();
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
