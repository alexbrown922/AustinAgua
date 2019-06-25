FRAISE 02/2011

Directory should include the following subdirectories and files:

main.R
data/Qdown_LODZ_172.txt
data/QF_LUCY_LODZ_172
util/FRAISE_Calculations
util/Functions.R
outputs/SAMPLE_FRAISE_OUTPUT.txt

1) Edit the R script "main.R" (User input section), 
   providing site characteristics and estimates for mean day time values of
   incoming radiant energy (QDOWN_MEAN) and anthropogenic heat flux (QF_MEAN). 
   In the current default settings both estimates are taken from data files in the data/ directory
   If only interested in flux ratio values (and UZE) simply assign QDOWN_MEAN <-0 and QF_MEAN <- 0

2) Source the script "main.R" from R

3) An output table similar to the file "SAMPLE_FRAISE_OUTPUT.txt" should be created in the directory outputs/

contact: thomas.loridan@gmail.com