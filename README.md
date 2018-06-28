# LightSheetImageAnalysis
This code assumes that raw fluorescence intensity traces have been generated using software (e.g. Imaris, Matlab) and that this has been exported to a Matlab variable "trace". Using this variable, Matlab can compute the normalized fluorescence:
Delta F = (F-F0)/F0
in the code "computeDeltaF.m". This function can be used to calculate a global value of F0 or a local value of F0.
The function "cellSynchrony.m" is used to identify calcium peaks and compute similarity indices to assess which cells are synchrnous with a template cell.
