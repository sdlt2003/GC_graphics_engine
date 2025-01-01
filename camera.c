
#include "camera.h"

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
