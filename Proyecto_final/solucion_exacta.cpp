#include <iostream>
#include <vector>
#include <cmath>
#include "parametros.h"

using namespace std;

// Parámetros del kink viajero
const double v = 0.5;        // Velocidad del solitón
const double x0_kink = 0.0;  // Posición inicial del kink
const double gamma_lorentz = 1.0 / sqrt(1.0 - v*v);  // Factor de Lorentz

// Solución analítica del kink viajero
double kink_exacto(double x, double t) {
    return 4.0 * atan(exp(gamma_lorentz * (x - v*t - x0_kink)));
}

int main() {
    cerr << "Generando solución exacta del kink viajero..." << endl;
    cerr << "Velocidad v = " << v << endl;
    cerr << "Factor gamma = " << gamma_lorentz << endl;
    
    // Crear malla para la solución
    vector<vector<double>> phi_exacta(Nx, vector<double>(Ny));
    
    // Imprimir datos a stdout para captura en Python
    // Formato: tiempo valor1 valor2 valor3 ...
    
    // Generar soluciones para diferentes tiempos
    for (int t_idx = 0; t_idx <= Nt; t_idx += 50) {
        double t_actual = t_idx * dt;
        
        // Calcular solución en toda la malla
        for (int i = 0; i < Nx; ++i) {
            double x = i*dx - Lx/2.0;
            for (int j = 0; j < Ny; ++j) {
                phi_exacta[i][j] = kink_exacto(x, t_actual);
            }
        }
        
        // Imprimir perfil en x (línea central en y)
        int j_mid = Ny / 2;
        cout << t_actual;
        for (int i = 0; i < Nx; ++i) {
            cout << " " << phi_exacta[i][j_mid];
        }
        cout << "\n";
    }
    
    cerr << "Solución generada." << endl;
    
    return 0;
}
