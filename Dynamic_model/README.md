# Cassava_model

This cassava leaf model is developed by Yu Wang(yuwangcn@illinois.edu), IGB, University of Illinois. The metabolic model part is developed based on the general C3 photosynthesis model (Zhu et al., 2007). If you want to know more details, please check the file’ Appendix cassava model equations and parameters.docx’

Quick start

Measured light intensity simulation:

RAC3leafMetaDriveLight(Lightinputfile,Pst,PRca)

An example command: RAC3leafMetaDriveLight('bon18182.dat',0,0)

High low light simulation:

RAC3leafMetaDriveLight2(Lightinputfile,Pst,PRca)

An example command: AC3leafMetaDriveLight2('Light150015.txt',0,0)


Input parameters:

Pst:
Pst=0 stomata conductance is calculated by steady state Ball-Berry model
Pst=1 Time dependent gs response, using ki and kd
PRca: 
PRca=0 Rubisco always actived
PRca=1 Consider Rubisco activation process

Cassava parameter file: 

cassavaP2.txt.

Cultivar	Vcmax25	Jmax25	kd	Ki	Ball-Berry Slop	Ball-Berry Intercept


Output: 

Transpiration(mol m-2 day-1)	 Photosynthesis(mol m-2 day-1)

