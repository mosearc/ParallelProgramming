#include <stdio.h>
#include <stdint.h> // Per uintptr_t
#include <stdlib.h>


void checkAlignment(void *ptr, const char *name) {
    uintptr_t address = (uintptr_t)ptr;
    int alignment = 1;

    // Trova il massimo allineamento possibile
    while (address % (alignment * 2) == 0) {
    	alignment *= 2;
    }

    printf("%s è allineato a %d byte.\n", name, alignment);
}

int main() {
    for (int n = 4; n<=4096; n=n*2){
    	// Allocazione normale
    	float *array = malloc(n * sizeof(float));
    	
	printf("%d : ", n);
    	// Controlla l'allineamento
    	checkAlignment(array, "array");

    	// Libera memoria
    	free(array);
    }

    return 0;
}

