#include <stdio.h>
#include <unistd.h>
#include <malloc.h>

int load_ppm(char *file, unsigned char **bufferptr, int *dimxptr, int *dimyptr)
{
  char line[40];
  int length, num;
  FILE *obj_file;
  int fd;

  if ((obj_file = fopen(file, "r")) == NULL)
  {
    *dimxptr = 0;
    *dimyptr = 0;
    *bufferptr = (unsigned char *)0;
    return (-1);
  }

  /* Lectura de formato de fichero */
  length = fscanf(obj_file, "%[^\n]\n", line);
  if (length > 1)
  {
    if ((line[0] == 'P') && (line[1] == '6'))
      printf("formato correcto\n");
    else
    {
      *dimxptr = 0;
      *dimyptr = 0;
      *bufferptr = (unsigned char *)0;
      printf("el formato tiene que ser P6\n");
      return (-1);
    }
  }

  /* fitxategi neurria irakurtzera */
  length = fscanf(obj_file, "%[^\n]\n", line);
  if (length > 0)
  {
    length = sscanf(line, "%d %d", dimxptr, dimyptr);
    if (length == 2)
      printf("dimensiones leídas: %d,%d\n", *dimxptr, *dimyptr);
    else
    {
      *dimxptr = 0;
      *dimyptr = 0;
      *bufferptr = (unsigned char *)0;
      printf("Fallos al leer las dimensiones\n");
      return (-1);
    }
  }
  /* Lectura de colores del fichero */
  length = fscanf(obj_file, "%[^\n]\n", line);
  if (length > 0)
  {
    length = sscanf(line, "%d", &num);
    if (length == 1)
      printf("Valor máximo de color leído: %d\n", num);
    else
    {
      *dimxptr = 0;
      *dimyptr = 0;
      *bufferptr = (unsigned char *)0;
      printf("Fallo en lectura de color máximo\n");
      return (-1);
    }
  }
  length = (*dimxptr) * (*dimyptr) * 3;
  *bufferptr = (unsigned char *)malloc(length);
  num = fread(*bufferptr, 1, length, obj_file);
  /*
  fd = fileno(obj_file);
  num = read(fd, bufferptr, length);
  */
  if (num != length)
  {
    *dimxptr = 0;
    *dimyptr = 0;
    free(*bufferptr);
    *bufferptr = (void *)0;
    printf("Errores al rellenarel buffer...num = %d, length =%d\n", num, length);
    return (-1);
  }
  else
  {
    printf("Se ha leído el buffer correctamente\n");
    return (1);
  }
}

