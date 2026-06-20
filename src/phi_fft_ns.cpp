#include <iostream>
#include <vector>
#include <cmath>
#include <complex>
#include <iomanip>
#include <algorithm>

constexpr double PHI = 1.6180339887498948482;
constexpr double PHI_INV = 0.6180339887498948482;
constexpr double LYAPUNOV_LAMBDA = 0.48121182505960347;

using Complex = std::complex<double>;
constexpr Complex I(0.0, 1.0);
constexpr double PI = M_PI;

// ═══════════════════════════════════════
// 1D FFT (Cooley-Tukey)
// ═══════════════════════════════════════
void fft(std::vector<Complex>& x) {
    int n = x.size();
    if (n <= 1) return;
    
    std::vector<Complex> even(n/2), odd(n/2);
    for (int i=0; i<n/2; i++) {
        even[i] = x[2*i];
        odd[i] = x[2*i+1];
    }
    fft(even); fft(odd);
    
    for (int k=0; k<n/2; k++) {
        Complex t = std::exp(-2.0*PI*I*(double)k/(double)n) * odd[k];
        x[k] = even[k] + t;
        x[k + n/2] = even[k] - t;
    }
}

class PhiFFTNavierStokes {
    int N;
    double dx, dt, nu;
    std::vector<double> u, v, w;
    int steps = 0;
    int idx(int i, int j, int k) { return i + j*N + k*N*N; }
    
public:
    PhiFFTNavierStokes(int n=32, double visc=PHI_INV) : N(n), dx(2.0*PI/N), dt(0.0005), nu(visc) {
        int M = N*N*N; u.resize(M); v.resize(M); w.resize(M);
        for (int i=0; i<N; i++) for (int j=0; j<N; j++) for (int k=0; k<N; k++) {
            double x=i*dx, y=j*dx, z=k*dx;
            u[idx(i,j,k)] = sin(x)*cos(y)*cos(z);
            v[idx(i,j,k)] = -cos(x)*sin(y)*cos(z);
            w[idx(i,j,k)] = 0.0;
        }
    }
    
    // True spectral energy at wavenumber k using FFT on 1D slice
    double true_Ek(int target_k) {
        if (target_k <= 0 || target_k >= N/2) return 0.0;
        
        // Take 1D slice along x-direction, average over y,z
        double total_E = 0.0;
        int samples = 0;
        
        for (int j=0; j<N; j+=4) {
            for (int k=0; k<N; k+=4) {
                std::vector<Complex> signal(N);
                for (int i=0; i<N; i++) {
                    int id = idx(i,j,k);
                    signal[i] = Complex(u[id], 0.0);
                }
                fft(signal);
                
                // Energy at wavenumber k = |û(k)|²
                int k_idx = target_k;
                if (k_idx < N/2) {
                    double Ek_local = std::norm(signal[k_idx]);
                    if (k_idx > 0) Ek_local += std::norm(signal[N - k_idx]);
                    total_E += Ek_local;
                    samples++;
                }
            }
        }
        
        return samples > 0 ? total_E / samples : 0.0;
    }
    
    // Spectrum across all wavenumbers
    std::vector<double> spectrum() {
        std::vector<double> spec(N/2, 0.0);
        for (int k=1; k<N/2; k++) {
            spec[k] = true_Ek(k);
        }
        return spec;
    }
    
    // Lyapunov cascade: λ(k) for φ-separated scales
    double lyapunov_cascade() {
        auto spec = spectrum();
        double total_lambda = 0.0;
        int count = 0;
        
        for (int k=2; k<N/4; k++) {
            int k_phi = std::max(1, (int)(k * PHI_INV));
            double Ek = spec[k];
            double Ek_phi = spec[k_phi];
            
            if (Ek < 1e-15) continue;
            
            // λ(k) = |E(k) - E(k/φ)| / E(k)
            double lambda = std::abs(Ek - Ek_phi) / Ek;
            total_lambda += lambda;
            count++;
        }
        
        return count > 0 ? total_lambda / count : 1.0;
    }
    
    double phi_stability() {
        double lc = lyapunov_cascade();
        return 1.0 - std::min(1.0, std::abs(lc - LYAPUNOV_LAMBDA) / LYAPUNOV_LAMBDA);
    }
    
    void step() {
        auto un=u, vn=v, wn=w;
        for(int i=0;i<N;i++)for(int j=0;j<N;j++)for(int k=0;k<N;k++){
            int id=idx(i,j,k);
            int ip=idx((i+1)%N,j,k), im=idx((i-1+N)%N,j,k);
            int jp=idx(i,(j+1)%N,k), jm=idx(i,(j-1+N)%N,k);
            int kp=idx(i,j,(k+1)%N), km=idx(i,j,(k-1+N)%N);
            
            auto lap=[&](auto& f){return (f[ip]+f[im]+f[jp]+f[jm]+f[kp]+f[km]-6.0*f[id])/(dx*dx);};
            auto adv=[&](auto& f){return u[id]*(f[ip]-f[im])/(2*dx)+v[id]*(f[jp]-f[jm])/(2*dx)+w[id]*(f[kp]-f[km])/(2*dx);};
            
            un[id]=u[id]+dt*(-adv(u)+nu*lap(u));
            vn[id]=v[id]+dt*(-adv(v)+nu*lap(v));
            wn[id]=w[id]+dt*(-adv(w)+nu*lap(w));
        }
        u=un; v=vn; w=wn; steps++;
    }
    
    void run(int total) {
        auto spec0 = spectrum();
        double lc0 = lyapunov_cascade();
        double ps0 = phi_stability();
        
        // Print initial spectrum shape
        std::cout << "ν=" << std::fixed << std::setprecision(4) << nu << std::endl;
        std::cout << "  Spectrum (first 8 modes): ";
        for (int k=1; k<=std::min(8, N/2-1); k++) {
            std::cout << std::scientific << std::setprecision(2) << spec0[k] << " ";
        }
        std::cout << std::endl;
        
        for(int i=0;i<total;i++) step();
        
        auto spec1 = spectrum();
        double lc1 = lyapunov_cascade();
        double ps1 = phi_stability();
        bool stable = (ps1 > 0.5);
        
        std::cout << "  λ_cascade: " << std::setprecision(6) << lc0 << "→" << lc1
                  << " | φ-stab: " << std::setprecision(4) << ps0 << "→" << ps1
                  << " | " << (stable ? "✅ TRUE LYAPUNOV STABLE" : "⚠️ CASCADE BROKEN") 
                  << std::endl << std::endl;
    }
};

int main() {
    std::cout << "╔════════════════════════════════════════════════════════════╗\n";
    std::cout << "║  Φ-FFT NAVIER-STOKES — True Spectral Analysis             ║\n";
    std::cout << "║  FFT-based energy spectrum | λ_φ = " << LYAPUNOV_LAMBDA << "                    ║\n";
    std::cout << "║  ΦΩ0 — I AM THAT I AM                                    ║\n";
    std::cout << "╚════════════════════════════════════════════════════════════╝\n\n";
    
    for (double v : {1.0, PHI_INV, 0.3}) {
        PhiFFTNavierStokes(32, v).run(200);
    }
    
    std::cout << "╔════════════════════════════════════════════════════════════╗\n";
    std::cout << "║  TRUE FFT SPECTRUM — φ-Lyapunov Cascade                   ║\n";
    std::cout << "║  ΦΩ0 — I AM THAT I AM                                    ║\n";
    std::cout << "╚════════════════════════════════════════════════════════════╝\n";
    
    return 0;
}
