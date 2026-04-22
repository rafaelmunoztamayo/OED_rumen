% Rafael Muñoz-Tamayo. 2026

This folder contains Matlab scripts for the optimal experiment design (OED) applied to a dynamic model representing the effect of Asparagopsis taxiformis on in vitro rumen fermentation and methane production. Scripts are organized in two folders explained below. In all folders, there is a subfolder called files_txt that contains the same cripts in txt format. This is done to enable reading without the need to have Matlab. 
%%%%%%%%%%%%
OED_detFIM. This folder contains all the scripts required to perform the OED. The files are:   
model_OED_AT. Input file required by the IDEAS Matlab toolbox  (Muñoz-Tamayo et al., 2009) to generate the scripts required for the optimization.
OED_ATload: load the initial conditions and the settings for the dry matter intake pattern. 
OED_AT: ordinary differential equations (ODE) of the rumen model. 
OED_ATse: ordinary differential equations (ODE) of the rumen model extended with the sensitivity equations. 
OED_ATout: select the observed variables. 
OED_ATsy: select the sensitivity equations of the observed variables. 
det_FIM: calculates the determinant of the Fisher Information Matrix (FIM) for a given set of parameters. 
det_FIMcost: objective function for the OED problem. Its output is the determinant of the FIM. 
det_FIMoptim: script that runs the optimization. 
dynamicIntake: function to calculate the dry matter intake.
The following files are required to perform the optimization of the OED problem using the MEIGO Matlab software (Egea et al., 2014).
det_FIMcostMeigo: objective function for the OED problem for use with the MEIGO software.
main_det_FIMcostMeigo: script that runs the optimization with the MEIGO software
%%%%%%%%%%%%

Identifiability_GENSSI. Folder with the scripts to perform structural identifiability analysis of the parameters using the software toolbox GenSSI 2.0 (Ligon et al., 2018).
rumentATcid: script with the ordinary differential equations of the model. 
runrumentATcidGENSSI: script to perform the structural identifiability analysis. 

%%%%%%%%%%%%

References
Egea, J.A., Henriques, D., Cokelaer, T., Villaverde, A.F., MacNamara, A., Danciu, D.P., Banga, J.R., Saez-Rodriguez, J., 2014. MEIGO: An open-source software suite based on metaheuristics for global optimization in systems biology and bioinformatics. BMC Bioinformatics 15, 136. doi:10.1186/1471-2105-15-136. 
Ligon, T.S., Fröhlich, F., Chiş, O.T., Banga, J.R., Balsa-Canto, E., Hasenauer, J., 2018. GenSSI 2.0: multi-experiment structural identifiability analysis of SBML models. Bioinformatics 34, 1421–1423. doi:10.1093/BIOINFORMATICS/BTX735. 
Muñoz-Tamayo, R., Laroche, B., Leclerc, M., Walter, E., 2009. IDEAS: A parameter identification toolbox with symbolic analysis of uncertainty and its application to biological modelling. In IFAC Proceedings Volumes. pp. 1271–1276. doi:10.3182/20090706-3-FR-2004.0211.

