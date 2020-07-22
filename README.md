# DOvEEgene
Endometrial and Ovarian Cancer prediction with machine learning methods


CreatInput.r:
Load data and create summarized input X matrix.

Prediction.r:
Make prediction of the test X matrix (from creatInput file) based on the Phase II trained model.

func_collapse.r:
Internal functions to collapse variants level data to individual level data. It will be called by CreatInput.r.

ROC_curve:
The ROC curve of the prediction performance on the new variants.


