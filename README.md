# retroreflectivity

## Summary
This document proposes modification to the conductor, dielectric and generalized Schlick BRDFs to allow reflection in the retro-reflective direction. The modification is based on the articles "BRDF Measurements and Analysis of Retroreflective Materials" (Belcour, et al 2014)[^Belcour2014] and "Properties of the Back-vector Modification for Microfacet BRDF Models" [^Raab2025]. Retro-reflection is phenomenon where by light is scattered back in the direction of the light source. It is commonly observed in materials used for safety such as traffic signs, road markings, safety clothing, vehicle license plates, and bicycle reflectors. Automotive, Industrial and apparel companies developing digital twins for simulation and visualization have expressed a need for materials that exhibit such behavior. The following proposes a change to existing BRDF models to allow reflection in the retro-reflective direction. 

### Additional Parameters to conductor_bsdf, dielectric_bsdf, and generalized_schlick_bsdf
| Parameters        |   Type    |   Default     |   Description             |
| ----------------- | --------- | --------------| ------------------------- |
| retroreflective   | bool      | false         | Enables Retro-Reflection  |

## BRDF
When enabled, the proposed BRDFs, construct a empirical glossy scattering lobe based on Belcour et. al.[^Belcour2014]. The paper describes a Back vector as the half vector of the reflected view direction and the light direction. It applies the Back vector, replacing H to several different empirical models, but for this implementation we focus on the microfacet model.

The Back vector is defined by Equation 3:

$$v' = 2(n \cdotp v)n - v$$

$$b = {v' + l\over ||v' + l||}$$  

The standard form of the microfacet distribution BRDF[^Walter2007] is given by Eq 20:

$$f = {D(h)F(v \cdotp h)G(l,v,h)\over 4(l \cdotp n)(v \cdotp n)}$$

Replacing H with B produces the retro reflective BRDF[^Belcour2014] Section 4.B.2:

$$f_b = {D_b(b)F_b(v \cdotp b)G(l,v,h)\over 4(l \cdotp n)(v \cdotp n)}$$

Belcour demonstrates this using a Beckmann distribution however for our implmentation we consider the work of Raab[^Raab2025], which discusses the plausibility and generalizes the Back vector substitution for standard microfacet BRDFs. Raab provides proofs for symmetry and energy conservation, as well as albedo equivalence in the retro and forward directions and sites plausible BRDFS that use Smith and V-Cavity masking functions, with Phong, Beckmann and GGX Distributions. Eq 7 - Eq 16.

Raab notes:

> Given an implementation of a regular microfacet BRDF, extending it to retro-reflection is straightforward:
> * Evaluation merely needs to replace v with v′ upfront.
> * Similarly, importance sampling of l given v can be realized by replacing v with v′ upfront and
> then importance sampling the regular microfacet BRDF. This may include low variance sampling
> using the domain of visible microfacets [Hd14].
> * As the albedos of standard BRDF and retro-reflective BRDF are essentially identical, compensating
> for energy loss in the sense of [KSK01] can be realized using the same data tables.
>
> — [Raab2025]

The proposed changes to **conductor_bsdf, generalized_schlick_bsdf, and dielectric_bsdf,**  use the GGX microfacet distribution [^Heitz2014] and adds the boolean **retrorefelctive**. When enabled, **v'** is substituted for **v**.


## References
[^Belcour2014]: Laurent Belcour, Romain Pacanowski, Marion Delahaie, Aude Laville-Geay, Laure Eupherte *BRDF Measurements and Analysis of Retroreflective Materials*, Journal of the Optical Society of America, 2014, JOSA A, 31 (12), pp.2561-2572. https://hal.science/hal-01083366/

[^Raab2025]: Matthias Raab NVIDIA *Properties of the Back-vector Modification for Microfacet BRDF Models* 

[^Walter2007]: B. WalterS. Marschner, H. Li, and K. Torrance, *Microfacet Models for Refraction through Rough Surfaces*, Proceedings of Eurographics Symposium on Rendering, 2007 pp. 195-206. https://www.graphics.cornell.edu/~bjw/microfacetbsdf.pdf

[^Heitz2014]: Eric Heitz, *Understanding the Masking-Shadowing Function in Microfacet-Based BRDFs*, Journal of Computer Graphics Techniques (JCGT), vol. 3, no. 2, 48-107, 2014