# Injector Design

The following scripts can be used to design configurations of coaxial swirl and jet injectors. This also includes code to analyze the dynamic stability of swirl injector elements and to calculate the resonant frequencies of the combustion chamber.

## Injector_Design.m
  - Choose the injector configuration [op_cl1, op_cl2, prop, eta1, eta2, mix, tao].
    
    This allows for the design of coaxial swirl elements (open or closed), coaxial shear elements, and any other combination of coaxial jet and swirl orifices.

  - Input the fluid parameters [mdot1, mdot2, deltap1, deltap2, rho1, rho2, n1, n2, xi1, xi2, wall].

    Set the hydraulic loss coefficients to zero for potential flow calculations (default value is 0.5). Pick the inner wall thickness based on structural requirements, usually 0.2mm - 1.0mm.

  - Input [alpha1, d, cd1, cd2].
    
    If a closed inner swirl configuration is chosen, input the inner element total spray angle. For internal mixing configuration 60-80 degrees is standard, for external mixing configurations larger angles may be chosen. If external mixing is chosen, input [d], the difference between the inner and outer spray angles. If jet elements are chosen, input the estimated discharge coefficients.

  - Select [alpha2].

    If an internal mixing outer closed swirl configuration is chosen, input the base outer element spray angle. The total spray angle of the internally mixed coaxial injector is usually 30-40 degrees less than the base angle of the outer element. 

  - Run script and open command window.

    For outputs to be valid the elements must be hydraulically independant, meaning that the inner element must be accomodated within the gas vortex of the outer element and the fluid from the inner element must arrive further than 2-3mm downstream of the inlets to the outer element for internal mixing configurations. The gas vortex calculations assume that all sections of the outer injector are adjacent to the nozzle of the inner injector. For the shear coaxial configuration, the inner-to-outer momentum flux ratio is also an output, ensure this is a suitable value.

## Dynamic_Response.m
  - Input swirl injector geometry inlcuding all length parameters [rn, rv, rt, ln, lv, lt, theta, n].

    The lengths can be initially chosen based on rough guidelines, lt=(3-4)rt, ln=(2)rn, lv=(2-3)rin, and then updated to generate a more stable response.

  - Select artificial viscosity coefficient and range for analysis.

    Input the throttle values, frequency range, and also the nominal pressure drop.

  - Adjust resonant mode calculation section if needed, run script, and open command window.

    An important parameter for ensuring stability is the ratio of injector pressure drop to combustion chamber pressure. As a rule of thumb, this should always remain above 8-10% to reduce the risk of low-frequency instability. Nominal values are 20-30%, but can be higher depending on feed system capabilities. 

## Combustion_Chamber.m
  - Input combustion chamber parameters and frequency threshold for outputs [a, R, L, threshold].

    The local speed of sound will change at different throttles. Calculations can be repeated to analyze more throttle values.

  - Run script and open command window.

    Return to Dynamic_Response.m and ensure that the injector resonant modes do not overlap with the combustion chamber modes. Stability can also be improved by adding baffles to the injector plate or even acoustic cavities.

## Swirl_Plots.m
  - Visualization script. Input range of parameters and check output plots to see how varying parameters affects injector characteristics.

## Inverse_Design.m
  - Input geometry and fluid characteristics. Can be used to reverse-engineer swirl elements.

## Injector_Tolerences.m
  - Input the initial geometry and target mass flow rate.

    Adjust tolerences based on design and feed system capabilities.

  - Input geometry range based on manufacturing methods and tolerences.

    Generally 0.1mm precision is acceptable. This corresponds to intervals of 0.05 (*1e-3) in the radial values.

  - Run script and open command window.

    If this script is used, check that the new geometry still forms a valid coaxial configuration.

## Unlike_Triplet.m and Unlike_Doublet.m
  - Input fluid parameters, collision angle, and/or tilt angles.

    The default collision angle is 60 degrees but other values can be chosen. The default [D] value is 5 but values between 3 and 7 are common, with lower values generally providing better mixing at the risk of melting the injector. Jet injector diameters usually fall within the range of 0.8-2.5mm.

  - Run script and open command window.

    Check that the area and momentum ratios are suitable values. If this is not the case, using a different type of impinging injector may be better. In many cases, an injector plate may consist of more than one type of impinging injector configuration. Non-zero resultant angles and non-unity momentum ratios can be used to shape the mass flux distribution of the spray if desired.













    

