#include<iostream>
#include<vector>
#include<cmath>
#include "parametros.h"
#include "condiciones_iniciales.h"

using namespace std;

// Parámetros del kink viajero (deben coincidir con condiciones_iniciales.cpp)
const double v = 0.5;
const double x0_kink = 0.0;
const double gamma_lorentz = 1.0 / sqrt(1.0 - v*v);

// Solución exacta del kink viajero
double kink_exacto(double x, double t) {
    return 4.0 * atan(exp(gamma_lorentz * (x - v*t - x0_kink)));
}

//Función para calcular la energía total
double calcular_energia(const vector<vector<double>>& curr, const vector<vector<double>>& old) {
    double energia = 0.0;
    for (int i = 1; i < Nx - 1; ++i) {
        for (int j = 1; j < Ny - 1; ++j) {
            double phi_t = (curr[i][j] - old[i][j]) / dt;
            double grad_x = (curr[i+1][j] - curr[i-1][j]) / (2.0 * dx);
            double grad_y = (curr[i][j+1] - curr[i][j-1]) / (2.0 * dy);
            
            energia += 0.5 * (phi_t * phi_t + grad_x * grad_x + grad_y * grad_y) + (1.0 - cos(curr[i][j]));
        }
    }
    return energia * dx * dy;
}

//Función para calcular la carga topológica

double calcular_carga(const vector<vector<double>>& curr) {
    double carga = 0.0;
    // Para un Kink que varía en X, integramos d_phi/dx a lo largo de una línea media en Y
    int j_mid = Ny / 2;
    for (int i = 1; i < Nx - 1; ++i) {
        carga += (curr[i+1][j_mid] - curr[i-1][j_mid]) / 2.0;
    }
    return carga / (2.0 * M_PI); // M_PI requiere #include <cmath>
}

int main() { 
    //Se reserva memoria para los tres pasos de tiempo
    vector<vector<double>> phi_old(Nx,vector<double>(Ny));
    vector<vector<double>> phi_curr(Nx, vector<double>(Ny));
    vector<vector<double>> phi_next(Nx, vector<double>(Ny));
    vector<vector<double>> phi_t(Nx, vector<double>(Ny));

    //inicializar
    inicializar_phi(phi_curr);
    inicializar_phi_t(phi_t);

    //Primer paso de tiempo 

    for (int i = 1; i < Nx - 1; ++i) {
        for (int j = 1; j < Ny - 1; ++j) {
            double lap = (phi_curr[i+1][j] - 2*phi_curr[i][j] + phi_curr[i-1][j])/(dx*dx) +
                         (phi_curr[i][j+1] - 2*phi_curr[i][j] + phi_curr[i][j-1])/(dy*dy);
            phi_old[i][j] = phi_curr[i][j] - dt * phi_t[i][j] + 0.5 * dt * dt * (lap - sin(phi_curr[i][j]));
        }
    }
    
    // Aplicar condiciones de frontera para phi_old (en t = -dt)
    double t_old = -dt;
    for (int i = 0; i < Nx; ++i) {
        double x = i*dx - Lx/2.0;
        for (int j = 0; j < Ny; ++j) {
            // Aplicar en todos los bordes
            if (i == 0 || i == Nx-1 || j == 0 || j == Ny-1) {
                phi_old[i][j] = kink_exacto(x, t_old);
            }
        }
    }

    //Bucle principal de tiempo 

    for (int t=0; t<Nt ; ++t){
        for (int i = 1; i < Nx - 1; ++i) {
            for (int j = 1; j < Ny - 1; ++j) {
                double lap = (phi_curr[i+1][j] - 2*phi_curr[i][j] + phi_curr[i-1][j])/(dx*dx) +
                             (phi_curr[i][j+1] - 2*phi_curr[i][j] + phi_curr[i][j-1])/(dy*dy);
                
                // Ecuación de evolución
                phi_next[i][j] = 2.0 * phi_curr[i][j] - phi_old[i][j] + 
                                 (dt * dt) * (lap - sin(phi_curr[i][j]));
            }
        }

        //condiciones de frontera: usar solución exacta en todos los bordes
        double t_actual = (t+1) * dt;
        
        for (int i = 0; i < Nx; ++i) {
            double x = i*dx - Lx/2.0;
            for (int j = 0; j < Ny; ++j) {
                // Aplicar en todos los bordes
                if (i == 0 || i == Nx-1 || j == 0 || j == Ny-1) {
                    phi_next[i][j] = kink_exacto(x, t_actual);
                }
            }
        }
        if (t%50 == 0){
            double E = calcular_energia(phi_curr,phi_old);
            double Q = calcular_carga(phi_curr);
            cerr << "Paso : " << t << "\t| Energía : " << E << "\t| Carga topológica: " << Q << endl;
            
            // Imprimir perfil en línea central a stdout
            int j_mid = Ny / 2;
            cout << t*dt;
            for (int i = 0; i < Nx; ++i) {
                cout << " " << phi_curr[i][j_mid];
            }
            cout << "\n";
        }
        //Actualización de punteros
        phi_old = phi_curr;
        phi_curr = phi_next;
    }

    cerr << "Simulación completa" << endl;
    return 0;
}  