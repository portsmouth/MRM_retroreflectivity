# retroreflectivity

## Summary
This document proposes a new MaterialX retro reflective BRDF node based on the article “BRDF Measurements and Analysis of Retroreflective Materials” (Belcour, et al 2014) [1]. Retro-reflection is phenomenon where by light is scattered back in the direction of the light source. It is commonly observed in materials used for safety such as traffic signs, road markings, safety clothing, vehicle license plates, and bicycle reflectors. Automotive, Industrial and apparel companies developing digital twins for simulation and visualization have expressed a need for materials that exhibit such behavior. The following describes a new node retro_reflection_bsdf and describes a modification to use a GGX microfacet distribution.


## BSDF
retro_reflection_bsdf: Constructs a glossy scattering lobe based on Belcour et. al.[^Belcour2014]. The paper describes a Back vector which is the half vector of the reflected view direction and the light direction. It applies the Back vector to a Beckmann distribution, replacing H with B, to create a new empirical retro-reflective BSDF. For out implmentation Beckmann has been replaced with GGX with supporting proof from Raab. [^Raab2025]

### Parameters
| Parameters    |   Type    |   Range                   |   Default             |   Description                             |
| ------------- | --------- | --------------------------| ----------------------| ----------------------------------------- |
| roughness     | float2    | [0, 1]                    | (0.3,0.3)             | Roughness coefficents                     |
| tint          | color3    | [(0,0,0), (1,1,1)]        | (1,1,1)               | Color weight to tint the reflected light  |
| tangent       | float3    | [(-1,-1,-1) to (1,1,1)]   | World Space Tangent   | Tangent vector of the surface             |


## References
[^Belcour2014]: Laurent Belcour, Romain Pacanowski, Marion Delahaie, Aude Laville-Geay, Laure Eupherte *BRDF Measurements and Analysis of Retroreflective Materials* Journal of the Optical Society of America, 2014, JOSA A, 31 (12), pp.2561-2572. https://hal.science/hal-01083366/

[^Raab2025]: Matthias Raab NVIDIA *Properties of the Back-vector Modification for Microfacet BRDF Models* 