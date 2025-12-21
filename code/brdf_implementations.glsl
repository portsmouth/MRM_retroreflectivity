// This file contains GLSL implementations of all of the relevant BRDFs.
// These functions can also be compiled as C++ using glm.

const float pi = 3.14159265f;
const float rcppi = 0.318309886f;

float sqr(float x)
{
    return x*x;
}

vec4 sample_cosine(float u1, float u2)
{
    float r = sqrt(u1);
    float theta = 2.0 * pi * u2;
    float x = r * cos(theta);
    float y = r * sin(theta);
    float z = sqrt(max(0.0, 1.0 - x*x - y*y));
    float pdf = abs(z) / pi;
    return vec4(x, y, z, pdf);
}

// Lambert BRDF
vec3 f_lambert(vec3 rho)
{
    return rho * rcppi;
}

vec4 sample_lambert(vec3 rho, vec3 woutputL, float u1, float u2,
                    vec3& f)
{
    vec4 s = sample_cosine(u1, u2);
    vec3 wi_local(s.x, s.y, s.z);
    f = f_lambert(rho);
    return s;
}

// full ON BRDF
vec3 f_ON(vec3 rho, float roughness, vec3 wi_local, vec3 wo_local)
{
    float sigma = roughness * (pi * 0.5f); // ON sigma

    // Apadters to match reference code:
    vec3 V = wi_local;
    vec3 L = wo_local;

    // Code below is from https://github.com/portsmouth/OpenPBR-viewer/blob/b7545a91dfcc854c3a711225ae0548621d5c3529/glsl/diffuse_brdf.glsl
    // TODO: Adjust to match the style of the other functions.

    // Full model
    float cosThetaI = V.z, sinThetaI = sqrt(1.0f - cosThetaI * cosThetaI);
    float cosThetaO = L.z, sinThetaO = sqrt(1.0f - cosThetaO * cosThetaO);
    float cosPhiDiff = 0.0f;
    if (sinThetaI > 0.0f && sinThetaO > 0.0f)
    {
        /* Compute cos(phiO-phiI) using the half-angle formulae */
        vec3 wi = V;
        vec3 wo = L;
        float sinPhiI = clamp(wi.y / sinThetaI, -1.0f, 1.0f),
              cosPhiI = clamp(wi.x / sinThetaI, -1.0f, 1.0f),
              sinPhiO = clamp(wo.y / sinThetaO, -1.0f, 1.0f),
              cosPhiO = clamp(wo.x / sinThetaO, -1.0f, 1.0f);
        cosPhiDiff = cosPhiI * cosPhiO + sinPhiI * sinPhiO;
    }

    float thetaI = acos(V.z),
          thetaO = acos(L.z),
          alpha = max(thetaI, thetaO),
          beta  = min(thetaI, thetaO);

    float sinAlpha, sinBeta, tanBeta;
    if (V.z > L.z)
    {
        sinAlpha = sinThetaO;
        sinBeta = sinThetaI;
        tanBeta = sinThetaI / max(1.0e-7f, cosThetaI);
    }
    else
    {
        sinAlpha = sinThetaI;
        sinBeta = sinThetaO;
        tanBeta = sinThetaO / max(1.0e-7f, cosThetaO);
    }

    float sigma2 = sigma * sigma;
    float tmp = sigma2 / (sigma2 + 0.09f),
          tmp2 = (4.0 / (pi * pi)) * alpha * beta,
          tmp3 = 2.0 * beta / pi;

    float C1 = 1.0f - 0.5f * sigma2 / (sigma2 + 0.33f),
          C2 = 0.45f * tmp,
          C3 = 0.125f * tmp * tmp2 * tmp2,
          C4 = 0.17f * sigma2 / (sigma2 + 0.13f);

    if (cosPhiDiff > 0.0f)
        C2 *= sinAlpha;
    else
        C2 *= sinAlpha - tmp3*tmp3*tmp3;

    /* Compute tan(0.5 * (alpha+beta)) using the half-angle formulae */
    float tanHalf = (sinAlpha + sinBeta) / (sqrt(1.0f - sinAlpha * sinAlpha) + sqrt(1.0f - sinBeta  * sinBeta));
    vec3 snglScat = rho * (C1 + cosPhiDiff * C2 * tanBeta + (1.0f - abs(cosPhiDiff)) * C3 * tanHalf);
    vec3 dblScat  = rho * rho * (C4 * (1.0f - cosPhiDiff*tmp3*tmp3));
    return (snglScat + dblScat) * rcppi; // Removed multiplication by cosThetaO
}

vec4 sample_ON(vec3 rho, float r, vec3 woutputL, float u1, float u2,
               vec3& f)
{
    vec4 s = sample_cosine(u1, u2);
    vec3 wi_local(s.x, s.y, s.z);
    f = f_ON(rho, r, vec3(s), woutputL);
    return s;
}

// QON BRDF
vec3 f_QON(vec3 rho, float roughness, vec3 wi_local, vec3 wo_local, bool interrefl)
{
    float sigma = roughness * (pi * 0.5f);                 // QON sigma
    float mu_i = wi_local.z;                               // input angle cos
    float mu_o = wo_local.z;                               // output angle cos
    float sigma2 = sigma * sigma;                          // QON sigma squared
    float Aq = !interrefl ? 1.0f - 0.5f * sigma2           // QON A coefficient
                                   / (sigma2 + 0.33f)
                          : 1.0f - 0.5f * sigma2           // discrepancy-compensated A
                                   / (sigma2 + 0.57f);
    float Bq = 0.45f * sigma2 / (sigma2 + 0.09f);          // QON B coefficient
    float s = dot(wi_local, wo_local) - mu_i * mu_o;       // QON s term
    float sovertq = s > 0.0f ? s / max(mu_i, mu_o) : 0.0f; // QON s/t
    vec3 f_q = (rho * rcppi) * (Aq + Bq * sovertq);        // single-scatt. BRDF
    return f_q;
}

vec4 sample_QON(vec3 rho, float r, vec3 woutputL, bool interrefl, float u1, float u2,
               vec3& f)
{
    vec4 s = sample_cosine(u1, u2);
    vec3 wi_local(s.x, s.y, s.z);
    f = f_QON(rho, r, vec3(s), woutputL, interrefl);
    return s;
}

// Fujii QON BRDF
vec3 f_fujiiQON(vec3 rho, float roughness, vec3 wi_local, vec3 wo_local)
{
    float sigma = roughness * (pi * 0.5f);                // QON sigma
    float mu_i = wi_local.z;                              // input angle cos
    float mu_o = wo_local.z;                              // output angle cos
    float sigma2 = sigma * sigma;                         // QON sigma squared
    vec3 Aq = vec3(1.0f - 0.5f * sigma2                   // Fujii QON A
                          / (sigma2 + 0.33f))
                        + 0.17f * rho * sigma2
                          / (sigma2 + 0.13f);
    float Bq = 0.45f * sigma2 / (sigma2 + 0.09f);         // QON B coefficient
    float s = dot(wi_local, wo_local) - mu_i * mu_o;      // QON s term
    float sovertF = s > 0.0f ? s / max(mu_i, mu_o) : s;   // Fujii QON s/t
    vec3 f_q = (rho * rcppi) * (Aq + vec3(Bq) * sovertF); // single-scatt. BRDF
    return f_q;
}

vec4 sample_fujiiQON(vec3 rho, float r, vec3 woutputL, float u1, float u2, vec3& f)
{
    float sigma = r * (pi * 0.5f);                // QON sigma
    vec4 s = sample_cosine(u1, u2);
    f = f_fujiiQON(rho, r, vec3(s), woutputL);
    return s;
}

const float constant1_FON = 0.5f - 2.0f / (3.0f * pi);
const float constant2_FON = 2.0f / 3.0f - 28.0f / (15.0f * pi);

// FON directional albedo (exact)
float E_FON_exact(float mu, float roughness)
{
    float sigma = roughness;                          // FON sigma prime
    float AF = 1.0f / (1.0f + constant1_FON * sigma); // FON A coeff.
    float BF = sigma * AF;                            // FON B coeff.
    float Si = sqrt(1.0f - (mu * mu));
    float G = Si * (acos(mu) - Si * mu)
            + (2.0f / 3.0f) * ((Si / mu) * (1.0f - (Si * Si * Si)) - Si);
    return AF + (BF * rcppi) * G;
}

// FON directional albedo (approx.)
float E_FON_approx(float mu, float roughness)
{
    float sigma = roughness; // FON sigma prime
    float mucomp = 1.0f - mu;
    const float g1 = 0.0571085289f;
    const float g2 = 0.491881867f;
    const float g3 = -0.332181442f;
    const float g4 = 0.0714429953f;
    float GoverPi = mucomp * (g1 + mucomp * (g2 + mucomp * (g3 + mucomp * g4)));
    return (1.0f + sigma * GoverPi) / (1.0f + constant1_FON * sigma);
}

// FON BRDF
vec3 f_FON(vec3 rho, float roughness, vec3 wi_local, vec3 wo_local)
{
    float sigma = roughness;                              // FON sigma prime
    float mu_i = wi_local.z;                              // input angle cos
    float mu_o = wo_local.z;                              // output angle cos
    float s = dot(wi_local, wo_local) - mu_i * mu_o;      // QON s term
    float sovertF = s > 0.0f ? s / max(mu_i, mu_o) : s;   // FON s/t
    float AF = 1.0f / (1.0f + constant1_FON * sigma);     // FON A coeff.
    return (rho * rcppi) * AF * (1.0f + sigma * sovertF); // single-scatter
}

vec4 sample_FON(vec3 rho, float r, vec3 woutputL, float u1, float u2,
               vec3& f)
{
    vec4 s = sample_cosine(u1, u2);
    vec3 wi_local(s.x, s.y, s.z);
    f = f_FON(rho, r, vec3(s), woutputL);
    return s;
}

// EON directional albedo (approx)
vec3 E_EON(vec3 rho, float roughness, vec3 wi_local, bool exact)
{
    float sigma = roughness;                           // FON sigma prime
    float mu_i = wi_local.z;                           // input angle cos
    float AF = 1.0f / (1.0f + constant1_FON * sigma);  // FON A coeff.
    float EF = exact ? E_FON_exact(mu_i, sigma):       // FON wi albedo (exact)
                       E_FON_approx(mu_i, sigma);      // FON wi albedo (approx)
    float avgEF = AF * (1.0f + constant2_FON * sigma); // average albedo
    vec3 rho_ms = (rho * rho) * avgEF / (vec3(1.0f) - rho * (1.0f - avgEF));
    return rho * EF + rho_ms * (1.0f - EF);
}

// EON BRDF
vec3 f_EON(vec3 rho, float roughness, vec3 wi_local, vec3 wo_local, bool exact)
{
    float sigma = roughness;                                   // FON sigma prime
    float mu_i = wi_local.z;                                   // input angle cos
    float mu_o = wo_local.z;                                   // output angle cos
    float s = dot(wi_local, wo_local) - mu_i * mu_o;           // QON s term
    float sovertF = s > 0.0f ? s / max(mu_i, mu_o) : s;        // FON s/t
    float AF = 1.0f / (1.0f + constant1_FON * sigma);          // FON A coeff.
    vec3 f_ss = (rho * rcppi) * AF * (1.0f + sigma * sovertF); // single-scatter lobe
    float EFo = exact ? E_FON_exact(mu_o, sigma):              // FON wo albedo (exact)
                        E_FON_approx(mu_o, sigma);             // FON wo albedo (approx)
    float EFi = exact ? E_FON_exact(mu_i, sigma):              // FON wi albedo (exact)
                        E_FON_approx(mu_i, sigma);             // FON wi albedo (approx)
    float avgEF = AF * (1.0f + constant2_FON * sigma);         // avg. albedo
    vec3 rho_ms = (rho * rho) * avgEF / (vec3(1.0f) - rho * (1.0f - avgEF));
    const float eps = 1.0e-7f;
    vec3 f_ms = (rho_ms * rcppi) * max(eps, 1.0f - EFo)        // multi-scatter lobe
                                 * max(eps, 1.0f - EFi)
                                 / max(eps, 1.0f - avgEF);
    return f_ss + f_ms;
}

///////////////////////////////////////////////////////////////////
// GGX BRDF
///////////////////////////////////////////////////////////////

float D_GGX(vec3 M, float r)
{
    const float tol = 1.0e-7f;
    float ax = r * r; // GGX alpha
    float ay = r * r; // GGX alpha
    float Ddenom = pi * ax * ay * sqr(sqr(M.x/ax) + sqr(M.y/ay) + sqr(M.z));
    return 1.0f / max(Ddenom, tol);
}

float lambda_ggx(vec3 V, float r)
{
    const float tol = 1.0e-7f;
    if (abs(V.z) < tol)
        return 0.0f;
    float a = sqr(r);
    return (-1.0 + sqrt(1.0 + (sqr(a*V.x) + sqr(a*V.y))/sqr(V.z))) / 2.0;
}

float G2_ggx(vec3 V, vec3 L, float r)
{
    // separable masking-shadowing function
    // return 1.0f/(1.0f + lambda_ggx(V, r)) * 1.0f/(1.0f + lambda_ggx(L, r));

    // height-correlated masking-shadowing function
    return 1.0f / (1.0f + lambda_ggx(V, r) + lambda_ggx(L, r));
}

vec3 Schlick(vec3 F0, float mu)
{
    return F0 + pow(1.f - mu, 5.f)*(vec3(1.f) - F0);
}

vec3 FresnelConductorF82(float mu, vec3 F0, vec3 F82)
{
    const float mu_bar = 1.f/7.f;
    const float denom = mu_bar * pow(1.f - mu_bar, 6.f);
    vec3 Fschlick_bar = Schlick(F0, mu_bar);
    vec3 Fschlick     = Schlick(F0, mu);
    return Fschlick - mu * pow(1.f - mu, 6.f) * (vec3(1.f) - F82) * Fschlick_bar / denom;
}

vec3 f_GGX(vec3 rho, float r, vec3 wi_local, vec3 wo_local)
{
    vec3 V = wi_local;
    vec3 L = wo_local;
    vec3 H = normalize(V + L);
    float D = D_GGX(H, r);
    float G2 = G2_ggx(V, L, r);
    const float tol = 1.0e-7f;
    float J = 1.f / max(4.f*abs(V.z)*abs(L.z), tol);
    vec3 F0  = vec3(0.19f);
    vec3 F82 = vec3(0.f);
    vec3 F = FresnelConductorF82(abs(dot(V, H)), F0, F82);
    return rho * F * G2 * D * J;
}

///////////////////////////////////////////////////////////////////
// Minimal Retro-reflective Microfacet (MRM) GGX BRDF
///////////////////////////////////////////////////////////////

vec3 f_MRM(vec3 rho, float r, vec3 wi_local, vec3 wo_local)
{
    vec3 V = wi_local;
    vec3 L = wo_local;
    V = vec3(-V.x, -V.y, V.z); // retro-reflect V
    vec3 B = normalize(V + L); // "back-vector" B
    float D = D_GGX(B, r);
    float G2 = G2_ggx(V, L, r);
    const float tol = 1.0e-7f;
    float J = 1.f / max(4.f*abs(V.z)*abs(L.z), tol);
    vec3 F0  = vec3(0.19f);
    vec3 F82 = vec3(0.f);
    vec3 F = FresnelConductorF82(abs(dot(V, B)), F0, F82);
    return rho * F * G2 * D * J;
}

///////////////////////////////////////////////////////////////////
// EON CLTC sampling
///////////////////////////////////////////////////////////////////

mat3 orthonormal_basis_ltc(vec3 w)
{
    float lenSqr = 1.f - w.z*w.z;
    vec3 X = lenSqr > 0.f ? vec3(w.x, w.y, 0.f) * inversesqrt(lenSqr) : vec3(1.f, 0.f, 0.f);
    vec3 Y = vec3(-X.y, X.x, 0.f); // cross(Z, X)
    return mat3(X, Y, vec3(0.f, 0.f, 1.f));
}

vec4 ltc_terms(float mu, float r)
{
    float a = 1.0 + r*(0.303392f + (-0.518982f + 0.111709f*mu)*mu + (-0.276266f + 0.335918f*mu)*r);
    float b = r*(-1.16407f + 1.15859f*mu + (0.150815f - 0.150105f*mu)*r)/(mu*mu*mu - 1.43545f);
    float c = 1.0f + (0.20013f + (-0.506373 + 0.261777f*mu)*mu)*r;
    float d = r*(0.540852f + (-1.01625f + 0.475392f*mu)*mu)/(-1.0743f + mu*(0.0725628f + mu));
    return vec4(a, b, c, d);
}

vec4 cltc_sample(vec3 woutputL, float r, float u1, float u2)
{
    vec4 abcd = ltc_terms(woutputL.z, r);  // coeffs of LTC $M$
    float a = abcd.x;
    float b = abcd.y;
    float c = abcd.z;
    float d = abcd.w;
    float R = sqrt(u1); float phi = 2.0f * pi * u2;          // CLTC sampling
    float x = R * cos(phi); float y = R * sin(phi);          // CLTC sampling
    float vz = 1.0 / sqrt(d*d + 1.f);                        // CLTC sampling factors
    float s = 0.5 * (1.f + vz);                              // CLTC sampling factors
    x = -mix(sqrt(1.f - y*y), x, s);                         // CLTC sampling
    vec3 wh = vec3(x, y, sqrt(max(1.f - (x*x + y*y), 0.f))); // $w_h$ sample via CLTC
    float pdf_wh = wh.z / (pi * s);                          // PDF of $w_h$ sample
    vec3 wi = vec3(a*wh.x + b*wh.z, c*wh.y, d*wh.x + wh.z);  // $M w_h$ (unnormalized)
    float len = length(wi);                                  // $|M w_h| = 1/|M^{-1} w_h|$
    float detM = c*(a - b*d);                                // $|M|$
    float pdf_wi = pdf_wh * len*len*len / detM;              // $w_i$ sample PDF
    mat3 fromLTC = orthonormal_basis_ltc(woutputL);          // transform $w_i$ to world space
    wi = normalize(fromLTC * wi);                            // transform $w_i$ to world space
    return vec4(wi, pdf_wi);
}

float cltc_pdf(vec3 woutputL, vec3 winputL, float r)
{
    mat3 toLTC = transpose(orthonormal_basis_ltc(woutputL));                 // transform $w_i$ to LTC space
    vec3 wi = toLTC * winputL;                                               // transform $w_i$ to LTC space
    vec4 abcd = ltc_terms(woutputL.z, r);  // coeffs of LTC $M$
    float a = abcd.x;
    float b = abcd.y;
    float c = abcd.z;
    float d = abcd.w;
    float detM = c*(a - b*d);                                                // $|M|$
    vec3 wh = vec3(c*(wi.x - b*wi.z), (a - b*d)*wi.y, -c*(d*wi.x - a*wi.z)); // $\mathrm{adj}(M) wi$
    float lensq = dot(wh, wh);                                               // $|M| |M^{-1} wi|$
    float vz = 1.f / sqrt(d*d + 1.f);                                        // CLTC sampling factors
    float s = 0.5f * (1.f + vz);                                              // CLTC sampling factors
    float pdf = (detM*detM)/(lensq*lensq) * max(wh.z, 0.f) / (pi * s);            // $w_i$ sample PDF
    return pdf;
}

vec3 uniform_lobe_sample(float u1, float u2)
{
    float z = u1;
    float R = sqrt(1.0 - z*z); float phi = 2.0 * pi * u2;
    float x = R * cos(phi); float y = R * sin(phi);
    return vec3(x, y, z);
}

vec4 sample_EON_COS(vec3 rho, float r, vec3 woutputL, bool exact,
                     float u1, float u2, vec3& f)
{
    vec4 s = sample_cosine(u1, u2);
    f = f_EON(rho, r, vec3(s), woutputL, exact);
    return s;
}

vec4 sample_EON_CLTC(vec3 rho, float r, vec3 woutputL, bool exact,
                     float u1, float u2, vec3& f)
{
    float mu = woutputL.z;
    float P_u = pow(r, 0.1) * (0.162925 + mu*(-0.372058 + (0.538233 - 0.290822*mu)*mu));
    float P_c = 1.0 - P_u;
    vec4 wip; float pdf_C;
    if (u1 <= P_u) {
        u1 = u1 / P_u;
        vec3 wi = uniform_lobe_sample(u1, u2);
        pdf_C = cltc_pdf(woutputL, wi, r);
        wip = vec4(wi, 1.f); }
    else {
        u1 = (u1 - P_u) / P_c;
        wip = cltc_sample(woutputL, r, u1, u2);
        pdf_C = wip.w; }

    const float pdf_U = 1.0 / (2.0 * pi);
    wip.w = P_u*pdf_U + P_c*pdf_C;
    f = f_EON(rho, r, vec3(wip), woutputL, exact);
    return wip;
}

///////////////////////////////////////////////////////////////////
