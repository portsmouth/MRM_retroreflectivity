# retroreflectivity

## Summary
This document proposes modification to the conductor, dielectric and generalized Schlick BSDFs to allow reflection in the retro-reflective direction. The modification is based on the article “BRDF Measurements and Analysis of Retroreflective Materials” (Belcour, et al 2014)[^Belcour2014]. Retro-reflection is phenomenon where by light is scattered back in the direction of the light source. It is commonly observed in materials used for safety such as traffic signs, road markings, safety clothing, vehicle license plates, and bicycle reflectors. Automotive, Industrial and apparel companies developing digital twins for simulation and visualization have expressed a need for materials that exhibit such behavior. The following describes BSDF model along with changes to implement a GGX microfacet distribution. Changes to the existing BSDF parameterization is also described along with the intested behavior.

### Parameters
The BSDF proposes a simple set of anisotropic roughness values, single scattering tint of the BSDF, and a reference tangent coordinate frame.

| Parameters    |   Type    |   Range                   |   Default             |   Description                             |
| ------------- | --------- | --------------------------| ----------------------| ----------------------------------------- |
| roughness     | float2    | [0, 1]                    | (0.3,0.3)             | Roughness coefficents                     |
| tint          | color3    | [(0,0,0), (1,1,1)]        | (1,1,1)               | Color weight to tint the reflected light  |
| tangent       | float3    | [(-1,-1,-1) to (1,1,1)]   | World Space Tangent   | Tangent vector of the surface             |

## BSDF
**retro_reflection_bsdf:** 

Constructs a empirical glossy scattering lobe based on Belcour et. al.[^Belcour2014] The paper describes a Back vector as the half vector of the reflected view direction and the light direction. It applies the Back vector to several different empirical models, but for this implementation we focus on the Beckmann microfacet model. Belcour replaces the half vector H with a new Back vector B, to create a new retro-reflective BSDF. 

For our implmentation we replace Beckmann with GGX [^Raab2025] so that it behaves similarly to the existing dielectric and generalized Schlick BSDFs.

The Back vector is defined by equation 3:

$$v' = 2(n \cdotp v)n - v$$

$$b = {v' + l\over ||v' + l||}$$  

The standard form of the Beckmann distribution BRDF[^Walter2007] is given by:

$$f = {D(h)F(v \cdotp h)G(l,v,h)\over 4(l \cdotp n)(v \cdotp n)}$$

Replacing H with B produces the retro reflective BRDF:

$$f_b = {D_b(b)F_b(v \cdotp b)G(l,v,h)\over 4(l \cdotp n)(v \cdotp n)}$$

where the distribution term D_b is:

$$D_b = {{1 \over πα^2cos^4_{RB}}} {exp (} {cos^2_{RB} - 1 \over α^2cos^2_{RB}}{)}$$

## References
[^Belcour2014]: Laurent Belcour, Romain Pacanowski, Marion Delahaie, Aude Laville-Geay, Laure Eupherte *BRDF Measurements and Analysis of Retroreflective Materials* Journal of the Optical Society of America, 2014, JOSA A, 31 (12), pp.2561-2572. https://hal.science/hal-01083366/

[^Raab2025]: Matthias Raab NVIDIA *Properties of the Back-vector Modification for Microfacet BRDF Models* 

[^Walter2007]: B. WalterS. Marschner, H. Li, and K. Torrance, *Microfacet Models for Refraction through Rough Surfaces* Proceedings of Eurographics Symposium on Rendering, 2007 pp. 195-206. https://www.graphics.cornell.edu/~bjw/microfacetbsdf.pdf