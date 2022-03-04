#include <stdio.h>
#include <stdbool.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>


void gravity_evaluator(int n, double *array_in, double *m_vec, double *array_out)
{
    int num_bodies;
    num_bodies = n / 6;
    int ioffset;
    int joffset;
    double c_cgs = -6.673 * pow(10, -8);
    long double dx, dy, dz, ax, ay, az, r;

    for (int i = 0; i < n; i++)
    {
        array_out[i] = 0;
    }
    for (int i = 0; i < num_bodies; i++)
    {
        ioffset = i * 6;
        for (int j = 0; j < num_bodies; j++)
        {
            joffset = j * 6;
            array_out[ioffset] = array_in[ioffset+3];
            array_out[ioffset+1] = array_in[ioffset+4];
            array_out[ioffset+2] = array_in[ioffset+5];
            if (i != j) {
                dx = array_in[ioffset] - array_in[joffset];
                dy = array_in[ioffset+1] - array_in[joffset+1];
                dz = array_in[ioffset+2] - array_in[joffset+2];
                r = pow(dx, 2) + pow(dy,2) + pow(dz,2);
                r = pow(r, .5);
                ax = (c_cgs * m_vec[j] / pow(r, 3)) * dx;
                ay = (c_cgs * m_vec[j] / pow(r, 3)) * dy;
                az = (c_cgs * m_vec[j] / pow(r, 3)) * dz;
                array_out[ioffset+3] += ax;
                array_out[ioffset+4] += ay;
                array_out[ioffset+5] += az;
            } 
        }
    }
}
