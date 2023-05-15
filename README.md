# This research archive contains all relevant information to reproduce the results in the study.
This study with the title *The performance of interrupted time series designs with a limited number of time points: Learning losses due to school closures during the COVID pandemic* consists of two ports. The first part is an motivating example using data from CITO. The second part is a simulation study based on the results of the motivating example. The full manucript can be found in the Manuscript folder.

## Motivating example
The data used in the motivating example is stems from CITO. The data is not available to the public. Refer to the data text file for information on how to obtain the data. To reproduce the results of the motivating example, the code in the scripts folder needs to be run. First the rmd file in the *1.1 pre-processing* folder can be run to transform the raw to data to the clean dataset with the necessary variables for a ITS design. Then the *1.2 Analysis* folder should be used, the rmd file named *1.2.1 Main_analyses* can be used to obtain the models presented in the study. *File 1.2.2* can be used to see extra analysis which considered, but not included in the study.

## Simulation study
To reproduce the results of the simulation study the code *2.1_Simulation* file should be run. The parameters used in the studies are given in the file and can be adjusted if need be. The functions behind the simulation study can be found in the Functions folder. The output of this code will be saved in in the output folder (results). The *2.2_Visualisation* can then be used to use the visualize the results. The figures can be created per simulation scenario.  

## Ethical approval
Furthermore, The study is approved by the Ethical Review Board of the Faculty of Social and Behavioural Sciences of Utrecht University. The approval is based on the documents sent by the researchers as requested in the form of the Ethics committee and filed under number 22-1805. The approval is valid through 01 June 2023. The approval of the Ethical Review Board concerns ethical aspects, as well as data management and privacy issues (including the GDPR).

## Technical Requirements
The technical requirements can be found in the Requirements.txt file. 

## License
This repository is accessible to the public on GitHub under the license type of GNU General Public License v3.0.

## Contact
The maintenance and public accessibility of this repository are managed by Florian van leeuwen. If you have any questions please feel free to contact me via email: f.dammesvanleeuwen@gmail.com.
