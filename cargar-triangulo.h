#ifndef CARGAR_TRIANGULO_H
#define CARGAR_TRIANGULO_H

#include "utils.h"

// @filename: nombre del fichero a cargar
// @trig_total: puntero a la variable donde se guardará el número de triángulos del archivo leído
// @hptrptr: puntero a la variable donde se guardará la dirección de memoria del primer triángulo leído
int cargar_triangulos(char *filename, int *trig_total, triangulo **first_trig);
int cargar_triangulos_color(char *filename, int *trig_total, triangulo **first_trig, unsigned char **rgbptr);

#endif