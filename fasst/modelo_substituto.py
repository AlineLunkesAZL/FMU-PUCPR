"""
MODELO SUBSTITUTO (placeholder) -- NAO E O FASST.

Serve so para testar o pipeline da FMU de ponta a ponta. Lucas e Allan devem
trocar este arquivo pelo FASST transcrito (ou apontar o adaptador para ele).

Balanco simplificado: folhagem (1 no) + substrato (N camadas) + base acoplada a laje.
"""
import numpy as np

SIGMA = 5.67e-8


class ModeloSubstituto:
    def __init__(self, LAI=2.0, h_sub=0.10, n_camadas=5, theta0=0.30, T0=293.15):
        self.LAI = LAI
        self.n = n_camadas
        self.dz = h_sub / n_camadas
        self.T = np.full(n_camadas, T0)      # temperatura do substrato por camada [K]
        self.theta = theta0                   # umidade volumetrica media [m3/m3]
        self.T_folha = T0
        self.k = 0.4                          # condutividade do substrato [W/mK]
        self.rho_c = 1.5e6                    # capacidade termica volumetrica [J/m3K]

    def passo(self, dt, T_ar, UR, G_sol, L_ceu, vento, chuva, T_laje):
        sig_f = 1.0 - np.exp(-0.5 * self.LAI)                 # fracao de cobertura
        h_conv = 5.7 + 3.8 * vento
        # evapotranspiracao simplificada (limitada pela umidade)
        LE = sig_f * 0.6 * G_sol * min(1.0, self.theta / 0.30) * (1.0 - UR * 0.5)
        # balanco da folhagem (quase-estacionario)
        R_liq = (1 - 0.2) * G_sol + L_ceu - 0.95 * SIGMA * self.T_folha**4
        self.T_folha = T_ar + (R_liq - LE) / (h_conv + 4 * 0.95 * SIGMA * T_ar**3)
        # fluxo no topo do substrato
        q_topo = (1 - sig_f) * (0.8 * G_sol) + sig_f * h_conv * (self.T_folha - self.T[0]) * 0.3 \
            + (1 - sig_f) * h_conv * (T_ar - self.T[0])
        # conducao 1D explicita com subpassos para estabilidade
        dt_max = 0.4 * self.rho_c * self.dz**2 / self.k
        nsub = max(1, int(np.ceil(dt / dt_max)))
        h = dt / nsub
        for _ in range(nsub):
            q = np.zeros(self.n + 1)
            q[0] = q_topo
            q[1:-1] = -self.k * np.diff(self.T) / self.dz
            q[-1] = -self.k * (T_laje - self.T[-1]) / (self.dz / 2)
            self.T += h * (q[:-1] - q[1:]) / (self.rho_c * self.dz)
        q_laje = q[-1]                                         # >0 = calor entrando na laje
        # balanco hidrico simplificado
        evap = LE / 2.45e6                                     # kg/m2s
        self.theta += dt * (chuva - evap) / (1000.0 * self.n * self.dz)
        escoamento = max(0.0, (self.theta - 0.45)) * 1000.0 * self.n * self.dz / dt
        self.theta = float(np.clip(self.theta, 0.05, 0.45))
        return dict(q_laje=float(q_laje), T_folha=float(self.T_folha),
                    T_sub=float(self.T.mean()), theta=self.theta, LE=float(LE),
                    escoamento=float(escoamento))
