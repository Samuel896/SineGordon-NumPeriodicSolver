#include "condiciones_iniciales.h"
#include "parametros.h"
#include <cmath>

// Inicializa el campo phi en t=0 (por ejemplo, todo en cero o con una perturbación)
void inicializar_phi(std::vector<std::vector<double>> &phi) {
    for (int i = 0; i < Nx; ++i) {
        for (int j = 0; j < Ny; ++j) {
            // Ejemplo: campo en cero
            phi[i][j] = 0.0;
            // Ejemplo alternativo: pequeña perturbación gaussiana en el centro
            // double x = i * dx;
            // double y = j * dy;
            // phi[i][j] = exp(-((x-Lx/2)*(x-Lx/2) + (y-Ly/2)*(y-Ly/2)));
        }
    }
}

// Inicializa la derivada temporal de phi en t=0 (por ejemplo, todo en cero)
void inicializar_phi_t(std::vector<std::vector<double>> &phi_t) {
    for (int i = 0; i < Nx; ++i) {
        for (int j = 0; j < Ny; ++j) {
            phi_t[i][j] = 0.0;
        }
    }
}