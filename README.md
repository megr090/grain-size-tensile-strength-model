# How to Use this Code
Meghana Ranganathan [meghanar@mit.edu] 

This repository contains code that pertains to the recently submitted paper "Recrystallization of ice enhances the creep and vulnerability to fracture of ice shelves", submitted to PNAS 3.6.21. It contains a steady-state grain size model applicable to glacier ice.

How to use this code: There are five things you can do with these codes:
(1) Compare the grain size model with data from three ice cores: WAIS Divide, GRIP, GISP2; you can just run the files "GrainSizeModel_GISP2.m", "GrainSizeModel_GRIP.m", "GrainSizeModel_WAISDivide.m" and it should produce all the plots from the paper and supplement

(2) Test the sensitivity of the grain size model to various parameters; you can run "ParameterSensitivityStudy_ActivationEnergies.m" to test the sensitivity to activation energy, "ParameterSensitivityStudy_D.m" to test the sensitivity to a characteristic length-scale, "ParameterSensitivityStudy_n.m" to test the sensitivity to the flow-law exponent, "ParameterSensitivityStudy_p.m" to test the sensitivity to the grain growth exponent, "ParameterSensitivityStudy_tc.m" to test the sensitivity to the critical temperature, and "ParameterSensitivityStudy_Theta.m" to test the sensitivity to the energy partitioning parameter. These should produce all the plots from the supplement.

(3) Run idealized versions of the steady-state grain size model: "GrainSizeEstimate_main.m" is the main code in which you can alter various parameters and compute ice temperature and grain size. Running this code will also plot the time to steady state in the idealized case and a figure of ice temperature, grain size, and stress. This code uses the following:
	-defineActivationEnergies.m
	-findGrainSize.m
	-findIceTemperature.m
	-plotGrainSizeIceTemperature.m

(4) Run an idealized version of the steady-state grain size model and compute tensile strength (which replicates Figure 3b of the paper) - this requires just running "computeplotShearMarginCase_Fracture.m"

(5) Compute grain size and ice temperature for Antarctic Ice Streams: "RunAntarcticIceStreams.m" is the main code in which you can alter parameters and choose the ice stream you want to look at. The code will then import various data through "readData.m" (surface velocity, surface strain-rates, surface mass balance, ice thickness) and compute ice temperature and grain size for the full domain defined in "readData.m". Then this plots some of these fields in "plotShearMarginProperties.m". This code uses the following:
	-findPropertiesofShearMargin.m
	-findShearMarginProperties_Antarctica.m
	-readData.m
	-findGrainSize.m
	-findIceTemperature.m
	-defineActivationEnergies.m
