#include <stdlib.h>
#include <stdio.h>
#include <math.h>

#define Nx 32 // интервалы по x
#define Ny 32 // интервалы по y
#define L 1.0 // длина кюветы
#define H 1.0 // выстока кюветы
#define A_tau 500  // Интенсивность тангенсальные течений
#define Mid 0.5 // центр открытой области
#define Length 1.0 // длина открытой области
#define x_in Mid - Length/2.0 // Крайняя левая точка открытой полости
#define x_out Mid + Length/2.0 // Крайняя правая точка открытой полости



#define _Pr 17.0
#define _Gr 1e5
#define _Re 1.0
#define _Ra 74000
#define _Ma 7800  
#define max_psi_time 500 // макс фиктивное время по пси
#define psi_eps 1e-8 
#define Nu_eps 1e-4
#define Theta0 20
#define Theta1 21

// Глобальные массивы
double T_old[Nx+1][Ny+1];
double T_half[Nx+1][Ny+1];
double T_new[Nx+1][Ny+1];

double psi_old[Nx+1][Ny+1];
double psi_half[Nx+1][Ny+1];
double psi_new[Nx+1][Ny+1];

double omega_old[Nx+1][Ny+1];
double omega_half[Nx+1][Ny+1];
double omega_new[Nx+1][Ny+1];

double u[Nx+1][Ny+1];
double v[Nx+1][Ny+1];


void save_isotherms(const char *filename, double data[Nx+1][Ny+1], double hx, double hy) {
    FILE *f = fopen(filename, "w");
    if (f == NULL) {
        printf("Ошибка открытия файла!\n");
        return;
    }
    // Заголовок (необязательно, но полезно)
    fprintf(f, "x\ty\tT\n");

    for (int i = 0; i <= Nx; i++) {
        for (int j = 0; j <= Ny; j++) {
            fprintf(f, "%lf\t%lf\t%lf\n", i * hx, j * hy, data[i][j]);
        }
        fprintf(f, "\n"); // Пустая строка для разделения сканов по X
    }
    fclose(f);
    printf("Данные для изотерм сохранены в %s\n", filename);
}

void save_vectors(const char *filename, double u[Nx+1][Ny+1], double v[Nx+1][Ny+1], double hx, double hy) {
    FILE *f = fopen(filename, "w");
    for (int i = 0; i <= Nx; i++) {
        for (int j = 0; j <= Ny; j++) {
            // Формат: X  Y  U  V
            fprintf(f, "%.20lf\t%.20lf\t%.20lf\t%.20lf\n", i * hx, j * hy, u[i][j], v[i][j]);
        }
    }
    fclose(f);
}

void save_surfer(const char *filename, double data[Nx+1][Ny+1], double hx, double hy){
    FILE *f = fopen(filename, "w");
    if (f == NULL){
        printf("Error opening\n");
        return;
    }

    for (int i =0; i <=Nx; i++){
        for (int j = 0; j <= Ny; j++){
            double x = i * hx;
            double y = j * hy;
            fprintf(f, "%lf\t%lf\t%lf\n", x, y, data[i][j] );
        }
        fprintf(f, "\n");
    }

    fclose(f);
    printf("Data saved!\n");

}

double calculate_Nu_bott(double T[Nx+1][Ny+1], double hx, double hy){
    double f[Nx + 1];
    for (int i = 0; i <= Nx; i++){
        f[i]= (-3.0 * T[i][0] + 4.0 * T[i][1] - T[i][2])/(2.0 * hy);
    }

    double sum = f[0] + f[Nx];

    for (int i = 1; i < Nx; i++){
        if (i % 2 != 0){
            sum += 4.0 * f[i]; 
        } else {
            sum += 2.0 * f[i];
        }

    }

    return -sum * (hx/3.0);

}

double calculate_Nu_mid(double T[Nx+1][Ny+1], double u[Nx+1][Ny+1], double hx, double hy){
    double f[Ny + 1];
    int mid = Nx / 2; 

    for (int j = 0; j <= Ny; j++) {

        double dTdx = (T[mid + 1][j] - T[mid - 1][j]) / (2.0 * hx);
        
        f[j] = dTdx - u[mid][j] * T[mid][j];
    }

    double sum = f[0] + f[Ny];

    for (int j = 1; j < Ny; j++) {
        if (j % 2 != 0) {
            sum += 4.0 * f[j]; 
        } else {
            sum += 2.0 * f[j];
        }
    }

    return sum * (hy / 3.0);
}

double calculate_Nu(double T[Nx+1][Ny+1], double hx, double hy){
    double f[Ny + 1];

    // производные на правой стенке (2 порядка )
    for (int j = 0; j <= Ny; j++) {
        f[j] = (3.0 * T[Nx][j] - 4.0 * T[Nx-1][j] + T[Nx-2][j]) / (2.0 * hx);
    }

    double sum = f[0] + f[Ny];

    for (int j = 1; j < Ny; j++) {
        if (j % 2 != 0) {
            sum += 4.0 * f[j]; // 2n+1
        } else {
            sum += 2.0 * f[j]; // 2n
        }
    }

    return sum * (hy / 3.0);
}
void adi_solve_psi(double psi_old[Nx+1][Ny+1], double psi_half[Nx+1][Ny+1], double psi_new[Nx+1][Ny+1], double omega_old[Nx+1][Ny+1], double tau_f, double hx, double hy  ){

    for (int n = 0; n < max_psi_time; n++)
    {
        double max_diff = 0.0;
            // массивы для коэффициентов 
        double ax[Nx+1];
        double bx[Nx+1];
        double cx[Nx+1];
        double dx[Nx+1];

        double alpha_x[Nx+1];
        double betta_x[Nx+1];

        double ay[Ny+1];
        double by[Ny+1];
        double cy[Ny+1];
        double dy[Ny+1];

        double alpha_y[Ny+1];
        double betta_y[Ny+1];

        // X-sweep:
        for (int j =1 ; j < Ny; j++){

            alpha_x[1] = 0;
            betta_x[1] = 0;

            for (int i = 1; i < Nx; i++ ){

                ax[i] = -tau_f / (2.0 * hx*hx);
                bx[i] = 1.0 + tau_f/(hx*hx);
                cx[i] = -tau_f / (2.0 * hx*hx);
                dx[i] = tau_f/2.0 * omega_old[i][j] + psi_old[i][j] + tau_f/(2.0* hy*hy) * (psi_old[i][j-1] - 2.0 * psi_old[i][j] + psi_old[i][j+1]);

                alpha_x[i+1] =  -cx[i] / (ax[i]*alpha_x[i] + bx[i]);
                betta_x[i+1] = (dx[i] - ax[i]*betta_x[i]) / (ax[i]*alpha_x[i] + bx[i]);
            }

            psi_half[Nx][j] = 0.0;

            for (int i = Nx-1; i >= 1; i--){
                psi_half[i][j] = alpha_x[i+1] * psi_half[i+1][j] + betta_x[i+1]; 
            }
            psi_half[0][j] = 0.0;
        }

        for (int i = 0; i <= Nx; i++){
            psi_half[i][0] = 0.0;
            psi_half[i][Ny] = 0.0;
        }
        for (int i = 1; i < Nx; i++){

            alpha_y[1] = betta_y[1] = 0.0;

            for (int j = 1; j < Ny; j++){
                ay[j] = -tau_f/(2.0*hy*hy);
                by[j] =  1.0 + tau_f/(hy*hy);
                cy[j] = -tau_f/(2.0*hy*hy);
                dy[j] = tau_f*omega_old[i][j] / 2.0 + psi_half[i][j] + tau_f/(2.0*hx*hx) * (psi_half[i-1][j] - 2.0 * psi_half[i][j] + psi_half[i+1][j]);
            
                alpha_y[j+1] = -cy[j] / (ay[j]*alpha_y[j] + by[j]);
                betta_y[j+1] = (dy[j] - ay[j]*betta_y[j])/(ay[j]*alpha_y[j] + by[j]);
            }

            psi_new[i][Ny] = 0.0; 

            for (int j = Ny-1; j >= 1; j--){
                psi_new[i][j] = alpha_y[j+1] * psi_new[i][j+1] + betta_y[j+1];
            }

            psi_new[i][0] = 0.0;
        }

        for (int j = 0; j <= Ny; j++){
            psi_new[0][j] = 0.0;
            psi_new[Nx][j] = 0.0;
        }



        for (int i = 1; i < Nx; i++) {
            for (int j = 1; j < Ny; j++) {
                double diff = fabs(psi_new[i][j] - psi_old[i][j]); // До переприсваивания!
                if (diff > max_diff) max_diff = diff;
            }
        }

        for (int i = 0; i <= Nx; i++){
            for (int j = 0; j <= Ny; j++){
                psi_old[i][j] = psi_new[i][j];
            }
        }
        
        // условие остановки по времени:
        if (max_diff < psi_eps){
            break;
        }
    }
}

void adi_solve_omega(double omega_old[Nx+1][Ny+1],double omega_half[Nx+1][Ny+1],double omega_new[Nx+1][Ny+1], double T_old[Nx+1][Ny+1], double psi_old[Nx+1][Ny+1], double u[Nx+1][Ny+1], double v[Nx+1][Ny+1], double hx, double hy, double tau, double Re){

    
    // массивы для коэффициентов 
    double ax[Nx+1];
    double bx[Nx+1];
    double cx[Nx+1];
    double dx[Nx+1];

    double alpha_x[Nx+1];
    double betta_x[Nx+1];

    double ay[Ny+1];
    double by[Ny+1];
    double cy[Ny+1];
    double dy[Ny+1];

    double alpha_y[Ny+1];
    double betta_y[Ny+1];

    // X-sweep:

    for (int j = 1; j < Ny; j++){
        alpha_x[1] = 0;
        betta_x[1] = -2.0 * ( psi_old[1][j])/(hx*hx);

        for (int i = 1; i < Nx; i++){
            ax[i] = -u[i][j]*tau/(2.0*hx*2.0) - tau/(2.0*hx*hx);
            bx[i] = 1.0 + tau/(hx*hx);
            cx[i] = u[i][j]*tau/(2.0*hx*2.0) - tau/(2.0*hx*hx);
            dx[i] = omega_old[i][j] + (tau/(2.0*hy*hy))*(omega_old[i][j-1] - 2.0*omega_old[i][j] + omega_old[i][j+1])  - (v[i][j] * tau / (4.0 * hy)) * (omega_old[i][j+1] - omega_old[i][j-1]) + tau/2.0 * _Ra/_Pr * (T_old[i+1][j] - T_old[i-1][j]) / (2.0 * hx);


            alpha_x[i+1] = -cx[i] / (ax[i]*alpha_x[i] + bx[i]);
            betta_x[i+1] = (dx[i] - ax[i]*betta_x[i])/(ax[i]*alpha_x[i] + bx[i]);
        }

        omega_half[Nx][j] = -2.0 * ( psi_old[Nx-1][j])/(hx*hx);

        for (int i = Nx-1 ; i>= 1; i--){
            omega_half[i][j] = alpha_x[i+1] * omega_half[i+1][j] + betta_x[i+1];
        }
        omega_half[0][j] = - 2.0 * psi_old[1][j] / (hx*hx);
    }
    for (int i = 0; i <= Nx; i++){
            double A = A_tau ;//* (hx * i);
            omega_half[i][0] = - 2.0 * psi_old[i][1] / (hy*hy);
            if (i == 0) omega_half[i][Ny] = _Ma/_Pr * ((T_old[i+1][Ny] - T_old[i][Ny]) / (hx)) - A;
            if (i == Nx) omega_half[i][Ny] = _Ma/_Pr * ((T_old[i][Ny] - T_old[i-1][Ny]) / (hx)) - A;
            else omega_half[i][Ny] = _Ma/_Pr * ((T_old[i+1][Ny] - T_old[i-1][Ny]) / (2.0 * hx)) - A;
    }
    // Y-sweep:
    for (int i = 1; i < Nx; i++){

        double xi = i * hx;

        alpha_y[1] = 0.0;
        betta_y[1] = -2.0*psi_old[i][1]/(hy*hy);

        for (int j = 1; j < Ny; j++){
            ay[j] = -v[i][j] * tau / (2.0*hy*2.0)  - tau/(2.0*hy*hy);
            by[j] = 1.0 + tau/(hy*hy);
            cy[j] = v[i][j] * tau / (2.0*hy*2.0)  - tau/(2.0*hy*hy);
            dy[j] = omega_half[i][j] + tau/(2.0*hx*hx) * (omega_half[i-1][j] - 2.0*omega_half[i][j] + omega_half[i+1][j])  - (u[i][j] * tau / (4.0 * hx)) * (omega_half[i+1][j] - omega_half[i-1][j]);

            alpha_y[j+1] = -cy[j] / (ay[j]*alpha_y[j] + by[j]);
            betta_y[j+1] = (dy[j] - ay[j]*betta_y[j])/(ay[j]*alpha_y[j] + by[j]);
        }
        double A = A_tau ; // * (hx * i);

        if (i == 0)omega_new[i][Ny] = _Ma/_Pr * ((T_old[i+1][Ny] - T_old[i][Ny]) / (hx)) - A;
        if (i == Nx) omega_new[i][Ny] = _Ma/_Pr * ((T_old[i][Ny] - T_old[i-1][Ny]) / (hx)) - A;
        else omega_new[i][Ny] = _Ma/_Pr * ((T_old[i+1][Ny] - T_old[i-1][Ny]) / (2.0 * hx)) - A;

        for (int j = Ny-1; j >= 1; j--){
            omega_new[i][j] = alpha_y[j+1] * omega_new[i][j+1] + betta_y[j+1];

        }
        omega_new[i][0] = -2.0 * psi_old[i][1] / (hy*hy);
    }

    for (int j = 0; j <= Ny; j++){
        omega_new[0][j] = -2.0 * psi_old[1][j] / (hx*hx);
        omega_new[Nx][j] = -2.0 * psi_old[Nx-1][j] / (hx*hx);
    }

}

void adi_solve_T(double T_new[Nx+1][Ny+1], double T_half[Nx+1][Ny+1], double T_old[Nx+1][Ny+1], double u[Nx+1][Ny+1], double v[Nx+1][Ny+1], double hx, double hy, double tau ){

    // массивы для коэффициентов 
    double ax[Nx+1];
    double bx[Nx+1];
    double cx[Nx+1];
    double dx[Nx+1];

    double alpha_x[Nx+1];
    double betta_x[Nx+1];

    double ay[Ny+1];
    double by[Ny+1];
    double cy[Ny+1];
    double dy[Ny+1];

    double alpha_y[Ny+1];
    double betta_y[Ny+1];

    // X-sweep:
    for (int j = 1; j < Ny; j++){

        alpha_x[1] = 1.0;
        betta_x[1] = 0.0;

        for (int i = 1; i < Nx; i++){

            ax[i] = -tau * u[i][j]/(4.0 * hx) - tau/(2.0*_Pr*hx*hx);
            bx[i] =  1.0 + tau/(hx*_Pr*hx);
            cx[i] = tau * u[i][j]/(4.0 * hx) - tau/(2.0*_Pr*hx*hx);
            dx[i] = T_old[i][j] + tau/(2.0 * _Pr * hy * hy) * (T_old[i][j+1] -  2.0 * T_old[i][j] + T_old[i][j-1]) - v[i][j] * tau / (4.0*hy) * (T_old[i][j+1] - T_old[i][j-1]);

            alpha_x[i+1] = -cx[i] / (ax[i]*alpha_x[i] + bx[i]);
            betta_x[i+1] = (dx[i] - ax[i]*betta_x[i])/(ax[i]*alpha_x[i] + bx[i]);
        }

        T_half[Nx][j] = T_half[Nx-1][j];
        for (int i = Nx-1; i >= 1; i--){
            T_half[i][j] = alpha_x[i+1] * T_half[i+1][j] + betta_x[i+1];
        }
        T_half[0][j] = T_half[1][j];
    }

    for (int i = 0; i <= Nx; i++){
        T_half[i][0] = Theta0;
        T_half[i][Ny] = Theta0 + (Theta1 - Theta0) * (i * hx) / L;;
    }

    // Y-sweep:

    for (int i = 1; i < Nx; i++){
        double xi = i * hx;

        alpha_y[1] = 0.0;
        betta_y[1] = Theta0;

        for (int j = 1; j < Ny; j++){

            ay[j] = -tau * v[i][j] / (4.0*hy) - tau/(2.0 * _Pr* hy*hy);
            by[j] = 1.0 + tau/(hy*_Pr*hy);
            cy[j] = tau * v[i][j] / (4.0*hy) - tau/(2.0 * _Pr * hy*hy);

            dy[j] = T_half[i][j] + tau/(2.0 * _Pr * hx * hx) * (T_half[i+1][j] - 2.0 * T_half[i][j] + T_half[i-1][j]) - tau * u[i][j]/(4.0 * hx) * (T_half[i+1][j] - T_half[i-1][j]);


            alpha_y[j+1] = -cy[j] / (ay[j]*alpha_y[j] + by[j]);
            betta_y[j+1] = (dy[j] - ay[j]*betta_y[j])/(ay[j]*alpha_y[j] + by[j]);

        }
        T_new[i][Ny] = Theta0 + (Theta1 - Theta0) * (i * hx) / L;

        for (int j = Ny-1; j >= 1; j--){
            T_new[i][j] = alpha_y[j+1] * T_new[i][j+1] + betta_y[j+1];
        }
        T_new[i][0] = Theta0;
    }

    for (int j = 0; j <= Ny; j++){

        T_new[0][j] = T_new[1][j];
        T_new[Nx][j] = T_new[Nx-1][j];

    }

}

int main(){

    double hx = (double)L/Nx;
    double hy = (double)H/Ny;
    double tau = (hx*hx)/4.0;
    int max_n = 500000;
    double tau_f = 0.001;
    double omega_eps = 1e-8;

    //init
    for (int i = 0; i <= Nx; i++){
        double xi = (double) hx * i;
        for (int j = 0; j <= Ny; j++){
            T_old[i][j] = 0.0;
            omega_old[i][j] = 0.0;
            psi_old[i][j] = 0.0;
            u[i][j] = 0.0;
            v[i][j] = 0.0;
        }
    }

    printf("Starting calculation with Ra = %.d, Ma = %.d\n", _Ra, _Ma);
    printf("Theta0 = %d, Theta1 = %d\n", Theta0, Theta1);

    for (int n = 1 ; n < max_n ; n++){
        // уравнение температуры
        adi_solve_T(T_new,T_half,T_old,u,v,hx,hy,tau);
        // уравнение вихря
        adi_solve_omega(omega_old,omega_half,omega_new, T_old, psi_old, u, v, hx, hy, tau, _Re);
        // уравнения тока
        adi_solve_psi(psi_old,psi_half,psi_new,omega_old,tau_f,hx,hy);

        // обновляем скорости 
        for (int i = 1; i < Nx; i++) {
            for (int j = 1; j < Ny; j++) {
                u[i][j] = (psi_new[i][j+1] - psi_new[i][j-1]) / (2.0 * hy);
                v[i][j] = -(psi_new[i+1][j] - psi_new[i-1][j]) / (2.0 * hx);
            }
        }


        double max_diff_omega = 0.0;
        for (int i = 0; i <= Nx; i++) {
            for (int j = 0; j <= Ny; j++) {
                double diff = fabs(omega_new[i][j] - omega_old[i][j]);
                if (diff > max_diff_omega) max_diff_omega = diff;
            }
        }
        
        // обновим массивы
        for (int i = 0; i <= Nx; i++){
            for (int j = 0; j <= Ny; j++){
                T_old[i][j] = T_new[i][j];
                omega_old[i][j] = omega_new[i][j];
                psi_old[i][j] = psi_new[i][j];
            }
        }

        if (n % 500 == 0) {
            printf("Step %d: max_diff_omega = %.2e\n", n, max_diff_omega);
        }


        if (n > 100 && max_diff_omega < omega_eps) {
            printf("\n--- CONVERGED at step %d ---\n", n);
            printf("Final max_diff_omega: %.2e\n", max_diff_omega);
            
            save_surfer("streamlines.dat", psi_new, hx, hy);
            save_surfer("vorticity.dat", omega_new, hx, hy);    
            save_vectors("vector.dat", u, v, hx, hy);
            save_isotherms("isoterm.dat", T_new, hx,  hy );
            break;
        }
    }

    return 0;
}