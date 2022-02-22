#include <stdio.h>
#include <stdbool.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>


void free_array(double **array, int rows, int cols)
{
    for (int j = 0; j < rows; j++){
        free(array[j]);
    }
    free(array);
}

