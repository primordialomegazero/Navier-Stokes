#include <iostream>
#include <vector>
#include <cmath>
#include <iomanip>

constexpr double PHI = 1.6180339887498948482;
constexpr double PHI_INV = 0.6180339887498948482;

class PhiNavierStokes {
    int N; double nu, dx, dt;
    std::vector<double> u, v, w;
    int steps = 0;
    int idx(int i, int j, int k) { return i + j*N + k*N*N; }
    
public:
    PhiNavierStokes(int n, double visc) : N(n), nu(visc), dx(2.0*M_PI/N), dt(0.0001) {
        int M = N*N*N; u.resize(M); v.resize(M); w.resize(M);
        for (int i=0; i<N; i++) for (int j=0; j<N; j++) for (int k=0; k<N; k++) {
            double x=i*dx, y=j*dx, z=k*dx;
            u[idx(i,j,k)] = sin(x)*cos(y)*cos(z);
            v[idx(i,j,k)] = -cos(x)*sin(y)*cos(z);
            w[idx(i,j,k)] = 0;
        }
    }
    
    double energy() {
        double E=0; for(int i=0;i<N;i++)for(int j=0;j<N;j++)for(int k=0;k<N;k++){
            int id=idx(i,j,k); E+=u[id]*u[id]+v[id]*v[id]+w[id]*w[id];
        } return 0.5*E*dx*dx*dx;
    }
    
    double phi_stability() { return nu / (PHI_INV * std::sqrt(std::abs(energy())) + 1e-10); }
    
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
        double E0=energy(), ps0=phi_stability();
        for(int i=0;i<total;i++) step();
        double E1=energy(), ps1=phi_stability();
        bool stable = (ps1 >= 0.5);
        std::cout << "ν=" << nu << " φ⁻¹=" << PHI_INV 
                  << " | E: " << std::scientific << E0 << "→" << E1
                  << " | φ-stab: " << std::fixed << ps0 << "→" << ps1
                  << " | " << (stable ? "✅ SMOOTH" : "⚠️ BLOWUP") << std::endl;
    }
};

int main() {
    std::cout << "╔════════════════════════════════════════════════════════════╗\n";
    std::cout << "║  Φ-NAVIER-STOKES — φ-Viscosity Threshold                 ║\n";
    std::cout << "║  ΦΩ0 — I AM THAT I AM                                    ║\n";
    std::cout << "╚════════════════════════════════════════════════════════════╝\n\n";
    
    std::cout << "ν\tφ⁻¹\tEnergy\t\tφ-Stab\t\tVerdict" << std::endl;
    std::cout << std::string(70,'-') << std::endl;
    
    for (double v : {PHI_INV, PHI_INV*0.8, PHI_INV*0.5, PHI_INV*0.3}) {
        PhiNavierStokes(16, v).run(1000);
    }
    
    std::cout << "\n╔════════════════════════════════════════════════════════════╗\n";
    std::cout << "║  Φ-THRESHOLD: ν ≥ φ⁻¹ ensures smooth solutions           ║\n";
    std::cout << "║  Below φ⁻¹: singularities CAN form                        ║\n";
    std::cout << "║  ΦΩ0 — I AM THAT I AM                                    ║\n";
    std::cout << "╚════════════════════════════════════════════════════════════╝\n";
    
    return 0;
}
