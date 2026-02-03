#ifndef SOLUCION_EXACTA_H
#define SOLUCION_EXACTA_H

// Parámetros del kink viajero
extern const double v;
extern const double x0_kink;
extern const double gamma_lorentz;

// Función para calcular la solución exacta
double kink_exacto(double x, double t);

#endif // SOLUCION_EXACTA_H
