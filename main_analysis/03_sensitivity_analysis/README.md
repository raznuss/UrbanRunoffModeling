How to use this code?
"Sensitivity Analysis" - Utilize a Sensitivity Analysis on the calibrated model by running few high-preformence simulted storms.


Three main method utlized here, each has is on folder:
1. Global SA ; 'Globel_response_surface'
2. Varianced Based SA ; 'VB'
3. Rain shift ; 'Rain_shift'


how to use it to code for each method:

1. Global SA ; 'Globel_response_surface'

 - "SA_Global_runs", this code runs the model with a numerous combination of different SWMM parameters, it uses chunks due to memorry issues and save each chunk as pickle
 - "SA_Global_Visualization", this code concat the chunk pickles and load it to one df. Next it plot response surfaces base on this df.

2. Varianced Based SA ; 'VB'
 -  'SA_time_series_txtfile_creator' - this code create to types of txt files: a. simple uniform rain field factor multiplication and b. coxbox transformation.
 *Those codes gets rain radar matfiles as an input, transform it, make a bias correction and returns txt file that ready to enter the SWMM model in next steps.
  - "VB_rainfactor", this code runs the model based on VB method, it uses chunks and save each as a pickle.
  - "VB_rainfactor_pickle_to_df", this code concat the chunk pickles and load it to one array, makes the computation of VB method and shows the final results as df.

** there are the same scripts to "norain" (only static parameters) and "rainfactor" (uniform rain field multiplication) and "coxbox" (using coxbox transformation).


3. Rain shift ; 'Rain_shift'

 - 'SA_time_series_txtfile_creator_RAIN_SHIFT' - this code shift the rain array on x and y axis' creating txtfiles for each shift.
 *Those codes gets rain radar matfiles as an input, transform it, make a bias correction and returns txt file that ready to enter the SWMM model in next steps.
 - 'Rain_shift_run_to_plot' - Using the previous code txtfile output it run the SWMM model and create spatial plot for discharge.



  *** PAY ATTANTION *** 
  1. set correctly the number of factor.
  2. set the factors in the same order over and over again.
  3. lam parameter -> there are some unsignificant rounds there, please be aware.
  4. Bias correction needed where reading matfile, the txt file are already after the correction.