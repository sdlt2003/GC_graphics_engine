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
#include "utils.h"
#include "transformations.h"
#include "camera.h"

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
