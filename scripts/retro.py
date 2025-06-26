
####################################################################
# Minimal Retroreflective Microfacet Model
####################################################################

import math
import numpy as np
import matplotlib.pyplot as plt
import colorsys

DENOM_TOLERANCE       = 1.0e-20
FLT_EPSILON           = 1.1920929e-7

def sqr(X): return X*X
def sqrt(X): return X**0.5

def D_ggx(M, r):
    ax = max(sqr(r), DENOM_TOLERANCE)
    ay = max(sqr(r), DENOM_TOLERANCE)
    Ddenom = math.pi * ax * ay * sqr(sqr(M[0]/ax) + sqr(M[1]/ay) + sqr(M[2]))
    return 1.0 / max(Ddenom, DENOM_TOLERANCE)


def ggx_lambda(V, r):
    if abs(V[2]) < FLT_EPSILON:
        return 0.0
    return (-1.0 + sqrt(1.0 + (sqr(r*V[0]) + sqr(r*V[1]))/sqr(V[2]))) / 2.0

def G2_ggx(V, L, m, r):
    return 1.0 / (1.0 + ggx_lambda(V, r) + ggx_lambda(L, r))

def normalize(v):
    norm = np.linalg.norm(v)
    if norm == 0:
       return v
    return v / norm

def brdf_ggx(V, L, r):
    H = normalize(V + L)
    D = D_ggx(H, r)
    G2 = G2_ggx(V, L, H, r)
    J = 1.0 / max(4.0*abs(V[2])*abs(L[2]), DENOM_TOLERANCE)
    return G2 * D * J

def brdf_raab(V, L, r):
    V = np.array([-V[0], -V[1], V[2]])
    H = normalize(V + L)
    D = D_ggx(H, r)
    G2 = G2_ggx(V, L, H, r)
    J = 1.0 / max(4.0*abs(V[2])*abs(L[2]), DENOM_TOLERANCE)
    return G2 * D * J

def brdf_belcour(V, L, r):
    Vp = np.array([-V[0], -V[1], V[2]])
    B = normalize(Vp + L)
    D = D_ggx(B, r)
    G2 = G2_ggx(V, L, B, r)
    J = 1.0 / max(4.0*abs(V[2])*abs(L[2]), DENOM_TOLERANCE)
    return G2 * D * J


def brdf_lambert(V, L, r):
    return 1.0/math.pi

def directional_albedo_numerical(r, V, BRDF):

    #  Integrate over theta_i
    #  Do solid angle integral for:
    #    pho(theta_o, phi_o) = \int_0^{pi/2} dtheta_i sin(theta_i) cos(theta_i)
    #                          \int_0^{2pi}  dphi_i   f(theta_i, phi_i, theta_o, phi_o)
    Ntheta = 128
    Nphi   = 256
    theta_array = np.linspace(0.0, math.pi/2.0, Ntheta)
    phi_array   = np.linspace(0.0, math.pi*2.0, Nphi)

    theta_integrand = np.empty(Ntheta) # integrand for theta integral
    for n_theta_i in range(Ntheta):

        theta_i = theta_array[n_theta_i]
        cos_theta_i = np.cos(theta_i)
        sin_theta_i = np.sin(theta_i)

        phi_integrand = np.empty(Nphi) # integrand for phi integral

        # do integral over phi_i
        for n_phi_i in range(Nphi):
            phi_i = phi_array[n_phi_i]
            cos_phi_i = np.cos(phi_i)
            sin_phi_i = np.sin(phi_i)
            L = np.array([sin_theta_i*cos_phi_i,
                          sin_theta_i*sin_phi_i,
                          cos_theta_i])
            phi_integrand[n_phi_i] = BRDF(V, L, r)

        phi_integral = np.trapz(phi_integrand, phi_array)
        theta_integrand[n_theta_i] = phi_integral * np.sin(theta_i) * np.cos(theta_i)

    # do integral over theta_i
    theta_integral = np.trapz(theta_integrand, theta_array)
    return theta_integral



# Compute directional albedo, for various V theta angles

def plot_albedo(BRDF, plot_name):

    V_theta_N  = 32
    V_theta_array = np.linspace(0.0, 0.99 * math.pi / 2.0, V_theta_N)

    N_r = 5
    r_array = np.linspace(0.1, 1.0, N_r)

    for n_r in range(N_r):

        # Iterate over roughness
        r = r_array[n_r]
        print('Running for roughness %f' % r)

        # Compute directional albedo for a range of V
        rho_array = np.empty(V_theta_N)

        for n_theta_i in range(V_theta_N):

            V_theta_i = V_theta_array[n_theta_i]
            print('   running for theta_i %d/%d: %f' % (n_theta_i, V_theta_N, V_theta_i))
            V_mu_i = np.cos(V_theta_i)
            V = np.array([np.sin(V_theta_i), 0.0, V_mu_i])

            # Compute directional albedo at given V
            rho = directional_albedo_numerical(r, V, BRDF)
            print('   rho = %f' % rho)
            rho_array[n_theta_i] = rho

        label_str = r'$r$=%3.2f' % r
        grayscale = r / r_array[-1]
        Clin = colorsys.hsv_to_rgb(1, 0.9, grayscale)
        plt.plot(V_theta_array, rho_array, label=label_str, color=Clin, linestyle='solid')

    plt.xlabel (r'View angle $\theta_o$ / 90$^{\circ}$', fontsize=14)
    plt.ylabel (r'Albedo $E(\theta_o)$', fontsize=14)

    plt.legend(loc="lower right")
    plt.savefig('albedo_%s.png' % plot_name)

plot_albedo(brdf_ggx,     'ggx')
plot_albedo(brdf_belcour, 'belcour')
plot_albedo(brdf_ggx,     'raab')

#plt.show()



