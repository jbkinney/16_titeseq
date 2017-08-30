# 16_titeseq
Processed data and analysis code for Adams et al., 2016 bioRxiv (http://www.biorxiv.org/content/early/2016/11/15/036335)  
## Data generation
To generate data based figures, run:  
```
./make_figs.sh
```
Figures are saved in the pdf folder.  
Figure 6 relies on pymol to generate the colored antibody structure. To generate:  
```
run ./make_figs.sh  
```
or
```
python figure_6_structure.py
```
Then,  
Open pymol  
From pymol open ./structure/figure_6.pse  
Run ./structure/pymol_color.py  
I then used the command
```
ray 800,800
```
to generate the png figure.  

The pymol scipt may not run if there are any " " characters in the filename.  

## Processed data  
Processed experimental data is stored in ./data  
Annotation for data follows:  
all_aas.txt - txt file of atoms in the 4-4-20 crystal structure  

bar_codes.xlsx - barcodes used in the tite-seq replicate experiments + which cell populations they were used to identify.  

CDR_library_July_5_2013_sequences.txt - list of CDR1H and CDR3H mutations used to create library.  

CDR1and3_Feb23_2015_pJK36_001.fcs - flow cytometry measurements of wt yeast displayed antibodies measured under saturating conditions.  

CDR1and3_Feb23_2015_pJK37_001_003.fcs - flow cytometry measurements of opt antibodies measured under saturating conditions.  

CDR1and3_Feb27_2015_pJK36_001.fcs - flow cytometry measurements of wt yeast displayed antibodies measured under saturating conditions. shown in figure S5.  

CDR1and3_Feb27_2015_pJK37_001_003.fcs - flow cytometry measurements of opt antibodies measured under saturating conditions. shown in figure S5.  

CDR1and3_Mar3_2015_pJK36_001.fcs - flow cytometry measurements of wt yeast displayed antibodies measured under saturating conditions.   

CDR1and3_Mar3_2015_pJK37_001_003.fcs  - flow cytometry measurements of super-optimized (OPT) yeast displayed antibodies measured under saturating conditions.   

CDR3_Feb12_2015_pJK36_001.fcs - flow cytometry measurements of opt antibodies measured under saturating conditions.  

CDR3_Feb12_2015_pJK37_002.fcs - flow cytometry measurements of wt yeast displayed antibodies measured under saturating conditions.   

CDR13_positions.txt - positions of just CDR1 and CDR3H atoms in the crystal structure. Shown in figure 6.  

fcs1 - flow cytometry measurements of the first tite-seq FACS sort. Shown in figure 3  
fcs2 - flow cytometry measurements of the second tite-seq FACS sort. Shown in figure S2  
fcs3 - flow cytometry measurements of the third tite-seq FACS sort. Shown in figure S3  

ideal_0_0.001.csv - simulated fits of poorly sampled data. Shown in figure S8.  
ideal_0_1.0.csv - simulated fits of sampled data. Shown in figure S8.  
ideal_0_10000.0.csv - simulated fits of massively sampled data. Shown in figure S8.  

replicate_1.csv - tite-seq results for first replicate. Shown in figure 4-6, S4, S6, S9-S12. Used to generate figure S7, S8.   
replicate_2.csv - tite-seq results for second replicate. Shown in figure 4-6, S4, S6, S9-S12.  
replicate_3.csv - tite-seq results for third replicate. Shown in figure 4-6, S4, S6, S9-S12.  

sort_counts_16.4.15.txt - FACS cell sort counts for first replicate. Shown in figure 3. Used to generate figure S7, S8.  
sort_counts_16.4.19.txt - FACS cell sort counts for second replicate. Shown in figure S2  
sort_counts_16.4.21.txt - FACS cell sort counts for third replicate. Shown in figure S3  

titration_curves.csv - flow cytometry titration curves and KD fits for isolated clones. Shown in figure 4, S6  

## Testing the fitting algorithm on simulated tite-seq data

To generate the fitting simulations for figure S8, run:  
```
python validation_simulations.py 0 0.001
python validation_simulations.py 0 1
python validation_simulations.py 0 10000
```
The output of these files are saved in ./data/ideal_0_0.001.csv  
./data/ideal_0_1.0.csv  
./data/ideal_0_10000.0.csv  

respectively. The first input argument (0) specifices the tite-seq replicate to simulate. This can be 0, 1, or 2. The second input argument specifices the sampling rate of the experiment. Higher numbers decrease the variability of fits. Each round of simulation fitting took 5 days to complate on the blacknblue cluster at Cold Spring Harbor Laboratory, August 2016.

## Running the binding affinity algorithm

To run the KD inference algorithm, you can include the line   
```
import x_star from KD_fit_log_poiss
```

in a python script, or you can run this script as  
```
python KD_fit_log_poiss.py input.dat output.dat
```
The input.dat file is fairly structured. An example of how to create one is created as example.dat when you run  
```
python KD_fit_log_poiss.py
```
without any inputs. Default test simulation figures are stored in ./fits/.
