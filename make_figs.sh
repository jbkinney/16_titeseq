#!/bin/bash

echo 'Removing old figures...'
rm pdfs/*

echo 'Running figure_1.py...'
python figure_1_illustration.py	

echo 'Running figure_3.py...'
python figure_3_procedure.py

echo 'Running figure_4.py...'
python figure_4_accuracy.py

echo 'Running figure_5.py...'
python figure_5_landscape.py

echo 'Running figure_6.py...'
python figure_6_structure.py

echo 'Running figure_S4_sensitivity.py...'
python figure_S4_sensitivity.py

echo 'Running figure_S5_facs.py...'
python figure_S5_facs.py

echo 'Running figure_S6_curves.py...'
python figure_S6_curves.py

echo 'Running figure_S7_facssim.py...'
python figure_S7_facssim.py

echo 'Running figure_S8_pipelinevalidation.py...'
python figure_S8_pipelinevalidation.py

echo 'Running figure_S9_histograms.py...'
python figure_S9_histograms.py

echo 'Running figure_S10_reproducibility.py...'
python figure_S10_reproducibility.py

echo 'Running figure_S11_composition.py...'
python figure_S11_composition.py

echo 'Running figure_S12_multipoint.py...'
python figure_S12_multipoint.py

echo 'Running figure_S13_synonymous_mean_variance.py...'
python figure_S13_synonymous_mean_variance.py

