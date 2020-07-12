# occulo

This repository contains data analysis scripts and data generated for the Closed-Eyes study.
This project is collaboration betwen the UBC Vision Lab and the Sleep Disorders Monitoring group in Mechanical Engineering.

## Instructions

Make sure your participant data is formatted and named correctly. XDF files should be located in a separate directory, structured as follows:

```
├───Data
    ├───P19
    │       P019_C1.xdf
    │       P019_C2.xdf
    │       P019_C3.xdf
    │       P019_C4.xdf
    │
    ├───P20
    │       P020_C1.xdf
    │       P020_C2.xdf
    │       P020_C3.xdf
    │       P020_C4.xdf
...
```

Ensure that no other files are stored in these directories.\

### Data Analysis
For data analysis, run **run_all.m**.

Relevant parameters are located under the Parameters section. These are the only values you will need to adjust.

```
%-------------------------------
% Parameters
opbci_threshold = [-inf inf]; % Threshold for filtering saccade
elink_threshold = [-inf inf];
outlier_iqr_scale = 1.5;

% Blinks
remove_blinks = false;

% For WH filter Q = 0.0005
% LR filter Q = 0.000005
% Else Q = 0.012
filter_types = {@Bandpass, @KFWesthInputAll, @LinearRecipoAll};
filter = filter_types{1};
Q = 0.000005;

excel_fname = "test.xlsx";
plot_folder_name = "test";
write_excel = true;
save_plots = true;
save_csv = false;

run_all_participants = false;
participants_to_run = ["P19", "P23"];
```
#### Parameters
* **opbci_threshold/elink_threshold**: saccade values beyond this range are filtered out. No threshold values for now.
* **outlier_iqr_scale**: For identifying outliers outside [value*interquartile range]. Default is 1.5.
* **remove_blinks**: Flatten blinks in signal. Leave at false.
* **filter**: Bandpass, LinearReciprocal and Westheimer.
* **Q**: Model noise for filter. Suggested values in comments.
* **excel_fname, plot_folder_name**: Names of Excel and Plot folder names if data is saved. Plot folder will automatically be created if it does not exist.
* **write_excel/save_plots**: Save Excel file and/or plots. Plots are saved in the specified folder.
* **save_csv**: Used when plotting single saccades.
* **run_all_participants/participants_to_run**: You can choose to run on all participants on on select participants. If **run_all_participants** is false, only the particpants in the **participants_to_run** array will be analysed.

#### Running the Script
Run **run_all.m**. When prompted, select the directory containing participant data (e.g. Data folder). The script can take up to 10 min to execute, if all participants are analysed.

### Plotting Saccades

To plot saccades, first run **run_all.m** with **save_csv** parameter set to **True** for the relevant participants. This will generate several csv files in the current directory, labelled with participant and condition. Then, you can run:
* **plot_indiv_saccades.m**: plots individual saccades of a chosen participant and condition. Eyelink, filtered and raw signals are overlayed. When prompted, navigate to and select the relevant xdf file, and then select the corresponding csv file. 

OR
* **plot_indiv_saccades_separate.m**: plots individual saccades of a chosen participant and condition. Eyelink, filtered and raw signals are plotted separately. When prompted, navigate to and select the relevant xdf file. Ensure csv files are in the same directory as this script.
 






