#include <stdio.h>
#include <malloc.h>
#define MAXLINE 200

/*
typedef struct punto
{
float x, y, z, u,v;
} punto;

typedef struct triangulo
{
punto p1,p2,p3;
} triangulo;
*/

#include "cargar-triangulo.h"

int cargar_triangulos(char *filename, int *trig_total, triangulo **first_trig)
{
	FILE *obj_file;
	char line[MAXLINE];
	int i, num_triangles;

	if ((obj_file = fopen(filename, "r")) == NULL)
	{
		*trig_total = 0;
		return (-1);
	}
	num_triangles = 0;
	while (fscanf(obj_file, "\n%[^\n]", line) > 0)
	{
		if (line[0] == 't') // triangulo!
		{
			num_triangles++;
		}
	}
	fclose(obj_file);
	*trig_total = num_triangles;
	
	// Asignación de memoria dinámica
	*first_trig = (triangulo *)malloc(num_triangles * sizeof(triangulo));

	obj_file = fopen(filename, "r");

	i = 0;
	while (fscanf(obj_file, "\n%[^\n]", line) > 0)
	{
		if (line[0] == 't') // triangulo!
		{
			sscanf(line + 2, "%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f", 
				   &((*first_trig)[i].p1.x), &((*first_trig)[i].p1.y), &((*first_trig)[i].p1.z),
				   &((*first_trig)[i].p1.u), &((*first_trig)[i].p1.v),
				   &((*first_trig)[i].p2.x), &((*first_trig)[i].p2.y), &((*first_trig)[i].p2.z),
				   &((*first_trig)[i].p2.u), &((*first_trig)[i].p2.v),
				   &((*first_trig)[i].p3.x), &((*first_trig)[i].p3.y), &((*first_trig)[i].p3.z),
				   &((*first_trig)[i].p3.u), &((*first_trig)[i].p3.v));
			i++;
		}
	}
	fclose(obj_file);
	return (1);
}

// zkop: número de valores exitosamente leídos por sscanf desde una línea del archivo, 
// utilizado para verificar si los datos de un triángulo están completos (15 valores) 
// o si solo incluyen coordenadas de vértices (9 valores).

int cargar_triangulos_color(char *filename, int *trig_total, triangulo **first_trig, unsigned char **rgbptr)
{
	FILE *obj_file;
	char line[MAXLINE];
	int i, num_triangles;
	int zkop;
	int r, g, b, color;
	float v1, v2, v3, v4, v5, v6, v7, v8, v9, v10, v11, v12, v13, v14, v15;

	if ((obj_file = fopen(filename, "r")) == NULL)
	{
		*trig_total = 0;
		printf("ezin ireki fitxategia\n");
		return (-1);
	}
	num_triangles = 0;
	color = 0;
	r = 255;
	g = 255;
	b = 255;
	while (fscanf(obj_file, "\n%[^\n]", line) > 0)
	{
		if (line[0] == 't') // triangulo!
		{
			num_triangles++;
		}
		if (line[0] == 'c') // color!
		{
			color = 1;
			sscanf(line + 2, "%d%d%d", &r, &g, &b);
		}
	}
	fclose(obj_file);
	*trig_total = num_triangles;
	*first_trig = (triangulo *)malloc(num_triangles * sizeof(triangulo));

	obj_file = fopen(filename, "r");

	i = 0;
	while (fscanf(obj_file, "\n%[^\n]", line) > 0)
	{
		if (line[0] == 't') // triangulo!
		{
			zkop = sscanf(line + 2, "%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f", &v1, &v2, &v3, &v4, &v5, &v6, &v7, &v8, &v9,
						  &v10, &v11, &v12, &v13, &v14, &v15);
			printf("lerroko balio kopurua = %d\n", zkop);
			if (zkop == 15)
			{
				(*first_trig)[i].p1.x = v1;
				(*first_trig)[i].p1.y = v2;
				(*first_trig)[i].p1.z = v3;
				(*first_trig)[i].p1.u = v4;
				(*first_trig)[i].p1.v = v5;
				(*first_trig)[i].p2.x = v6;
				(*first_trig)[i].p2.y = v7;
				(*first_trig)[i].p2.z = v8;
				(*first_trig)[i].p2.u = v9;
				(*first_trig)[i].p2.v = v10;
				(*first_trig)[i].p3.x = v11;
				(*first_trig)[i].p3.y = v12;
				(*first_trig)[i].p3.z = v13;
				(*first_trig)[i].p3.u = v14;
				(*first_trig)[i].p3.v = v15;
			}
			else
			{
				if (zkop == 9)
				{
					color = 1;
					(*first_trig)[i].p1.x = v1;
					(*first_trig)[i].p1.y = v2;
					(*first_trig)[i].p1.z = v3;
					(*first_trig)[i].p2.x = v4;
					(*first_trig)[i].p2.y = v5;
					(*first_trig)[i].p2.z = v6;
					(*first_trig)[i].p3.x = v7;
					(*first_trig)[i].p3.y = v8;
					(*first_trig)[i].p3.z = v9;
					(*first_trig)[i].p1.u = 0.0;
					(*first_trig)[i].p1.v = 0.0;
					(*first_trig)[i].p2.u = 0.0;
					(*first_trig)[i].p2.v = 0.0;
					(*first_trig)[i].p3.u = 0.0;
					(*first_trig)[i].p3.v = 0.0;
				}
				else
				{
					free(first_trig);
					fclose(obj_file);
					printf("errorea: %s...\n", filename);
					return (-1);
				}
			}
			i++;
		}
	}

	fclose(obj_file);
	if (color)
	{
		for (i = 0; i < num_triangles; i++)
		{
			(*first_trig)[i].p1.u = 0.0;
			(*first_trig)[i].p1.v = 0.0;
			(*first_trig)[i].p2.u = 0.0;
			(*first_trig)[i].p2.v = 0.0;
			(*first_trig)[i].p3.u = 0.0;
			(*first_trig)[i].p3.v = 0.0;
		}
		*rgbptr = (unsigned char *)malloc(3 * sizeof(unsigned char));
		(*rgbptr)[0] = r;
		(*rgbptr)[1] = g;
		(*rgbptr)[2] = b;
		return (9);
	}
	else
		return (15);
}

/* * /
void main(int argc, char*argv[])
{
int num_triangles,i;
triangulo *tptr;
unsigned char kolore[3];

printf("%s fitxategitik triangeluak kargatzera\n",argv[1]);
cargar_triangulos_color(argv[1], &num_triangles, &tptr,&(kolore[0]));
for (i=0; i<num_triangles; i++)
	{
	printf("t: (%.1f, %.1f, %.1f) (%.1f, %.1f)", tptr[i].p1.x, tptr[i].p1.y, tptr[i].p1.z, tptr[i].p1.u, tptr[i].p1.v);
	printf("   (%.1f, %.1f, %.1f) (%.1f, %.1f)", tptr[i].p2.x, tptr[i].p2.y, tptr[i].p2.z, tptr[i].p2.u, tptr[i].p2.v);
	printf("   (%.1f, %.1f, %.1f) (%.1f, %.1f)\n", tptr[i].p3.x, tptr[i].p3.y, tptr[i].p3.z, tptr[i].p3.u, tptr[i].p3.v);
	}
printf("color: %d, %d,%d\n",kolore[0],kolore[1],kolore[2]);
}


 // */
