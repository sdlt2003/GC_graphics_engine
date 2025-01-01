
#include "transformations.h"
#include "utils.h"
#include "camera.h"

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