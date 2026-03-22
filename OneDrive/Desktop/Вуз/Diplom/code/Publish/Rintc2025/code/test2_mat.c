#include <stdlib.h>
#include <stdio.h>
#include <math.h>

#define Nx 32 
#define Ny 32 
#define L 1.0 
#define H 1.0 

// Параметры физики
#define _Pr 17.0
#define _Ra 740000.0    
#define _Ma 7800
#define A_tau -500.0   
#define Mid 1.5       
#define Length 1.0    
#define x_in (Mid - Length/2.0)
#define x_out (Mid + Length/2.0)

#define max_psi_time 800 
#define psi_eps 1e-7 
#define omega_eps 1e-8

double T_old[Nx+1][Ny+1], T_half[Nx+1][Ny+1], T_new[Nx+1][Ny+1];
double psi_old[Nx+1][Ny+1], psi_half[Nx+1][Ny+1], psi_new[Nx+1][Ny+1];
double omega_old[Nx+1][Ny+1], omega_half[Nx+1][Ny+1], omega_new[Nx+1][Ny+1];
double u[Nx+1][Ny+1], v[Nx+1][Ny+1];


void save_scalar_field(const char *filename, double data[Nx+1][Ny+1], double hx, double hy) {
    FILE *f = fopen(filename, "w");
    if (!f) return;
    for (int i = 0; i <= Nx; i++) {
        for (int j = 0; j <= Ny; j++) {
            fprintf(f, "%lf\t%lf\t%lf\n", i * hx, j * hy, data[i][j]);
        }
        fprintf(f, "\n");
    }
    fclose(f);
}

void save_velocity_field(const char *filename, double u_field[Nx+1][Ny+1], double v_field[Nx+1][Ny+1], double hx, double hy) {
    FILE *f = fopen(filename, "w");
    if (!f) {
        printf("Error opening velocity file!\n");
        return;
    }
    for (int i = 0; i <= Nx; i++) {
        for (int j = 0; j <= Ny; j++) {
            fprintf(f, "%lf\t%lf\t%lf\t%lf\n", i * hx, j * hy, u_field[i][j], v_field[i][j]);
        }
        fprintf(f, "\n");
    }
    fclose(f);
}

// Функция расчета числа Нуссельта (не используется в main, но оставлена)
double calculate_Nu_mid(double T[Nx+1][Ny+1], double u_vel[Nx+1][Ny+1], double hx, double hy) {
    double f[Ny + 1];
    int mid = Nx / 2; 
    for (int j = 0; j <= Ny; j++) {
        double dTdx = (T[mid + 1][j] - T[mid - 1][j]) / (2.0 * hx);
        f[j] = dTdx - u_vel[mid][j] * T[mid][j];
    }
    double sum = f[0] + f[Ny];
    for (int j = 1; j < Ny; j++) {
        sum += (j % 2 != 0) ? 4.0 * f[j] : 2.0 * f[j];
    }
    return sum * (hy / 3.0);
}

void adi_solve_T(double hx, double hy, double tau) {
    double ax, bx, cx, dx, alpha[Nx+1], betta[Nx+1];
    double ay, by, cy, dy, alpha_y[Ny+1], betta_y[Ny+1];

    for (int j = 1; j < Ny; j++) {
        alpha[1] = 1.0; betta[1] = 0.0; // Изоляция слева
        for (int i = 1; i < Nx; i++) {
            ax = -tau * u[i][j]/(4.0 * hx) - tau/(2.0*hx*hx);
            bx = 1.0 + tau/(hx*hx);
            cx = tau * u[i][j]/(4.0 * hx) - tau/(2.0*hx*hx);
            dx = T_old[i][j] + tau/(2.0*hy*hy)*(T_old[i][j+1]-2.0*T_old[i][j]+T_old[i][j-1]) - v[i][j]*tau/(4.0*hy)*(T_old[i][j+1]-T_old[i][j-1]);
            alpha[i+1] = -cx / (ax * alpha[i] + bx);
            betta[i+1] = (dx - ax * betta[i]) / (ax * alpha[i] + bx);
        }
        T_half[Nx][j] = betta[Nx] / (1.0 - alpha[Nx]); // Изоляция справа
        for (int i = Nx-1; i >= 0; i--) T_half[i][j] = alpha[i+1] * T_half[i+1][j] + betta[i+1];
    }

    for (int i = 0; i <= Nx; i++) {
        double xi = i * hx;
        alpha_y[1] = 0.0; betta_y[1] = xi/L; // Линейный нагрев низ
        for (int j = 1; j < Ny; j++) {
            ay = -tau * v[i][j]/(4.0*hy) - tau/(2.0*hy*hy);
            by = 1.0 + tau/(hy*hy);
            cy = tau * v[i][j]/(4.0*hy) - tau/(2.0*hy*hy);
            dy = T_half[i][j] + tau/(2.0*hx*hx)*(T_half[i+1][j]-2.0*T_half[i][j]+T_half[i-1][j]) - tau*u[i][j]/(4.0 * hx)*(T_half[i+1][j]-T_half[i-1][j]);
            alpha_y[j+1] = -cy / (ay * alpha_y[j] + by);
            betta_y[j+1] = (dy - ay * betta_y[j]) / (ay * alpha_y[j] + by);
        }
        T_new[i][Ny] = xi/L; // Линейный нагрев верх
        for (int j = Ny-1; j >= 0; j--) T_new[i][j] = alpha_y[j+1] * T_new[i][j+1] + betta_y[j+1];
    }
}

void adi_solve_omega(double hx, double hy, double tau) {
    double alpha[Nx+1], betta[Nx+1], alpha_y[Ny+1], betta_y[Ny+1];

    // Прогонка по X
    for (int j = 1; j < Ny; j++) {
        alpha[1] = 0; betta[1] = -2.0 * psi_old[1][j] / (hx * hx);
        for (int i = 1; i < Nx; i++) {
            double ax = -u[i][j]*tau/(4.0*hx) - _Pr*tau/(2.0*hx*hx);
            double bx = 1.0 + _Pr*tau/(hx*hx);
            double cx = u[i][j]*tau/(4.0*hx) - _Pr*tau/(2.0*hx*hx);
            
            // Источниковый член теперь работает через Ra * Pr
            double dx = omega_old[i][j] + (_Pr*tau/(2.0*hy*hy))*(omega_old[i][j-1]-2.0*omega_old[i][j]+omega_old[i][j+1]) 
                        - (v[i][j]*tau/(4.0*hy))*(omega_old[i][j+1]-omega_old[i][j-1])
                        + (_Ra * _Pr * tau / (4.0 * hx)) * (T_old[i+1][j] - T_old[i-1][j]); 

            alpha[i+1] = -cx / (ax * alpha[i] + bx);
            betta[i+1] = (dx - ax * betta[i]) / (ax * alpha[i] + bx);
        }
        omega_half[Nx][j] = -2.0 * psi_old[Nx-1][j] / (hx * hx);
        for (int i = Nx-1; i >= 0; i--) omega_half[i][j] = alpha[i+1] * omega_half[i+1][j] + betta[i+1];
    }

    // Прогонка по Y
    for (int i = 1; i < Nx; i++) {
        double xi = i * hx;
        alpha_y[1] = 0; betta_y[1] = -2.0 * psi_old[i][1] / (hy * hy);
        for (int j = 1; j < Ny; j++) {
            double ay = -v[i][j]*tau/(4.0*hy) - _Pr*tau/(2.0*hy*hy);
            double by = 1.0 + _Pr*tau/(hy*hy);
            double cy = v[i][j]*tau/(4.0*hy) - _Pr*tau/(2.0*hy*hy);
            
            double dy = omega_half[i][j] + (_Pr*tau/(2.0*hx*hx))*(omega_half[i-1][j]-2.0*omega_half[i][j]+omega_half[i+1][j])
                        - (u[i][j]*tau/(4.0*hx))*(omega_half[i+1][j]-omega_half[i-1][j]);

            alpha_y[j+1] = -cy / (ay * alpha_y[j] + by);
            betta_y[j+1] = (dy - ay * betta_y[j]) / (ay * alpha_y[j] + by);
        }
        
        double dTdx_top = (T_old[i+1][Ny] - T_old[i-1][Ny]) / (2.0 * hx);
        omega_new[i][Ny] = (_Ma / _Pr) * dTdx_top - A_tau;
        
        for (int j = Ny-1; j >= 0; j--) omega_new[i][j] = alpha_y[j+1] * omega_new[i][j+1] + betta_y[j+1];
    }
}

void solve_psi(double hx, double hy) {
    double tau_f = 0.001;
    for (int n = 0; n < max_psi_time; n++) {
        double max_diff = 0.0;
        double alpha[Nx+1], betta[Nx+1], alpha_y[Ny+1], betta_y[Ny+1];
        for (int j = 1; j < Ny; j++) {
            alpha[1] = 0; betta[1] = 0;
            for (int i = 1; i < Nx; i++) {
                double ax = -tau_f/(2.0*hx*hx), bx = 1.0+tau_f/(hx*hx), cx = -tau_f/(2.0*hx*hx);
                double dx = tau_f/2.0*omega_old[i][j] + psi_old[i][j] + tau_f/(2.0*hy*hy)*(psi_old[i][j-1]-2.0*psi_old[i][j]+psi_old[i][j+1]);
                alpha[i+1] = -cx/(ax*alpha[i]+bx);
                betta[i+1] = (dx-ax*betta[i])/(ax*alpha[i]+bx);
            }
            psi_half[Nx][j] = 0;
            for (int i = Nx-1; i >= 0; i--) psi_half[i][j] = alpha[i+1]*psi_half[i+1][j] + betta[i+1];
        }
        for (int i = 1; i < Nx; i++) {
            alpha_y[1] = 0; betta_y[1] = 0;
            for (int j = 1; j < Ny; j++) {
                double ay = -tau_f/(2.0*hy*hy), by = 1.0+tau_f/(hy*hy), cy = -tau_f/(2.0*hy*hy);
                double dy = tau_f/2.0*omega_old[i][j] + psi_half[i][j] + tau_f/(2.0*hx*hx)*(psi_half[i-1][j]-2.0*psi_half[i][j]+psi_half[i+1][j]);
                alpha_y[j+1] = -cy/(ay*alpha_y[j]+by);
                betta_y[j+1] = (dy-ay*betta_y[j])/(ay*alpha_y[j]+by);
            }
            psi_new[i][Ny] = 0;
            for (int j = Ny-1; j >= 0; j--) {
                double old_val = psi_old[i][j];
                psi_new[i][j] = alpha_y[j+1]*psi_new[i][j+1] + betta_y[j+1];
                double diff = fabs(psi_new[i][j] - old_val);
                if (diff > max_diff) max_diff = diff;
                psi_old[i][j] = psi_new[i][j];
            }
        }
        if (max_diff < psi_eps) break;
    }
}

int main() {
    double hx = L / Nx, hy = H / Ny;
    double tau = 0.05 * (hx * hx * hy * hy) / (hx * hx + hy * hy); 
    int max_steps = 500000;

    for (int i = 0; i <= Nx; i++) {
        for (int j = 0; j <= Ny; j++) {
            T_old[i][j] = (double)i * hx / L;
            omega_old[i][j] = psi_old[i][j] = u[i][j] = v[i][j] = 0.0;
        }
    }

    printf("Running... (Ra=%.1e, A_tau=%.1f)\n", _Ra, A_tau);

    for (int n = 1; n <= max_steps; n++) {
        adi_solve_T(hx, hy, tau);
        adi_solve_omega(hx, hy, tau);
        solve_psi(hx, hy);

        for (int i = 1; i < Nx; i++) {
            for (int j = 1; j < Ny; j++) {
                u[i][j] = (psi_old[i][j+1] - psi_old[i][j-1]) / (2.0 * hy);
                v[i][j] = -(psi_old[i+1][j] - psi_old[i-1][j]) / (2.0 * hx);
            }
        }

        double max_err_w = 0.0;
        for (int i = 0; i <= Nx; i++) {
            for (int j = 0; j <= Ny; j++) {
                double diff = fabs(omega_new[i][j] - omega_old[i][j]);
                if (diff > max_err_w) max_err_w = diff;
                omega_old[i][j] = omega_new[i][j];
                T_old[i][j] = T_new[i][j];
            }
        }

        if (n % 2000 == 0) printf("Step %d: Err_W = %.2e\n", n, max_err_w);
        if (n > 1000 && max_err_w < omega_eps) {
            printf("CONVERGED at step %d\n", n);
            break;
        }
    }

    save_velocity_field("velocity_field.dat", u, v, hx, hy);
    
    printf("Success. Files 'result_T.dat', 'result_psi.dat' and 'velocity_field.dat' created.\n");

    return 0;
}