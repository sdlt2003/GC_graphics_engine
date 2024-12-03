
typedef struct punto
{
    float x, y, z, u, v;
} punto;

typedef struct vector
{
    float x, y, z;
} vector;

typedef struct triangulo
{
    punto p1, p2, p3;
    vector N;
} triangulo;

// @filename: nombre del fichero a cargar
// @trig_total: puntero a la variable donde se guardará el número de triángulos del archivo leído
// @hptrptr: puntero a la variable donde se guardará la dirección de memoria del primer triángulo leído
int cargar_triangulos(char *filename, int *trig_total, triangulo **first_trig);
int cargar_triangulos_color(char *filename, int *trig_total, triangulo **first_trig, unsigned char **rgbptr);
