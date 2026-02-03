#include "condiciones_iniciales.h"
#include "parametros.h"
#include <cmath>

// Parámetros del kink viajero (deben coincidir con solucion_exacta.cpp)
const double v = 0.5;        // Velocidad del solitón
const double x0_kink = 0.0;  // Posición inicial del kink
const double gamma_lorentz = 1.0 / std::sqrt(1.0 - v*v);  // Factor de Lorentz

// Inicializa el campo phi en t=0 con la solución exacta del kink viajero
void inicializar_phi(std::vector<std::vector<double>> &phi) {
    for (int i = 0; i < Nx; ++i) {
        double x = i*dx - Lx/2.0;
        for (int j = 0; j < Ny; ++j) {
            // Solución exacta del kink viajero en t=0
            phi[i][j] = 4.0 * atan(exp(gamma_lorentz * (x - x0_kink)));
        }
    }
}

// Inicializa la derivada temporal de phi en t=0 según la solución exacta
void inicializar_phi_t(std::vector<std::vector<double>> &phi_t) {
    for (int i = 0; i < Nx; ++i) {
        double x = i*dx - Lx/2.0;
        for (int j = 0; j < Ny; ++j) {
            // Derivada temporal de la solución exacta en t=0
            double arg = gamma_lorentz * (x - x0_kink);
            double sech_arg = 1.0 / cosh(arg);
            phi_t[i][j] = -4.0 * gamma_lorentz * v * sech_arg * sech_arg / (1.0 + exp(-2.0*arg));
        }
    }
}