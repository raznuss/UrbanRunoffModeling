How to use the code in this file!

The goal is to take a list of flood events and extract radar files that relate to the events.

The folder should contain:
functions: 'statCalc','CreateASCII'
gauge data: 'gaugeData_updated_042023'


the steps:
do it one time for discharge bigger than 13 and one time for discharge value 6.5-13:

1. Create xlxs file with list of flood events. for ex 'high_dischrage_and_dates_raanana'.
2. Use'RainGaugesCorrelationRaanana' to find how good is the correlation between radar data and gauge data in this events and area.
The corrolection and bias correction based on this domain:
xcordmin = 185000;
xcordmax = 197000;
ycordmin = 668000;
ycordmax = 685000;

3. Use 'BoundStormsInRadarData' to create the 'event_sturcture' file that select the relevent times for radar data extraction
4. Use 'CreateRadar_MATandASCII' to create the radar files by 'event_sturcture' file.
5. Use 'MissingRadarTimeCorrection' to correct mat files where time is missing, by doing linear interpolation.
PAY ATTANTION!! next steps includs simple manual work:
** The storm '20191227' need to be copy manuualy from the row matfile folder to the correct matfile folder because she does not need coreection.
** DELETE the storms: 20151025, 2016/12/01, 2016/12/02, 2018/01/26 matfiles from the row and correct folder

6. run Basin_Radar_Overleap, Bias correction needed where reading matfile, the txt file are after the correction
this is the logic that lead me to the Bias values, most of them from the matlab fig ploted in Use'RainGaugesCorrelationRaanana' code.
EVENTS_NAMES_L 

[('20120113.mat', 1.56),
 ('20130106.mat', 1.13),
 ('20141214.mat', 0.92),
 ('20151007.mat', 0.75),
 ('20151027.mat', 0.6),       אין ביאס, מציבים ביאס כמו ביום שאחריו, כי זה ביאס פר אירוע
 ('20151028.mat', 0.6),
 ('20151107.mat', 0.61),
 ('20151214.mat', 1.01),
 ('20151218.mat', 0.87),
 ('20160108.mat', 1.24),
 ('20160126.mat', 1.46),
 ('20161213.mat', 1.26),
 ('20161219.mat', 1.17),
 ('20161227.mat', 1.2),
 ('20171121.mat', 0.9),
 ('20171224.mat', 1.33),
 ('20180101.mat', 1.29),
 ('20180105.mat', 1.3),      נבדק נקודתית, הביאס הוא 1.3
 ('20180114.mat', 1.27),
 ('20180123.mat', 1.43),
 ('20180125.mat', 1.48),
 ('20180213.mat', 0.7),
 ('20181207.mat', 0.81),
 ('20191213.mat', 0.87),
 ('20191227.mat', 0.79),
 ('20200119.mat', 1.02)]












