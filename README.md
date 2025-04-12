# Supplementary Materials for "A hybrid mixture approach for clustering and characterizing cancer data"

This repository contains supplementary figures and plots referenced in the manuscript titled\
**"A hybrid mixture approach for clustering and characterizing cancer data"** by *Kazeem Kareem* and *Fan Dai*, submitted to Bioinformatics.

Due to space constraints in the main paper, the figures here provide additional insight into model performance and analyses.

---

## 📄 Contents

### 1. `Supplementary_materials_for_gmmfad_paper.pdf`

- Additional plots providing further insights into the simulation and real data results.  It also include codes for fitting the real datasets. 

### 2. `main`

- Contains all the source codes for the algorithm and simulations in the paper.
    - `breast_cancer_data_constq_variedq.R`-- contans code for fitting GMMFAD and GMMFADq on the breast cancer dataset.
    - `em4gmm.R` -- This is the primary source code, containing the functions GMMFAD  and other relevant functions.
    - `gmfad_emmix_emEM_p350_sim_time.R`-- Code for the simulation setup for comparing the time speedup, frobsenius errors, and ARI for GMMFAD and EMMIX algorithms for n=p=150. 
    - `gmfad_emmix_emEM_sim_time.R`-- Code for the simulation setup for comparing the time speedup, frobsenius errors, and ARI for GMMFAD and EMMIX algorithms for n=300, p=10. 
    - `gmfad_qq` -- This is the primary source code, containing the functions GMMFADq, and other relevant functions.
    - `lymphoma gene data (1).R`--  Code for fitting the lymphoma data with GMMFADq
    - `sim_BIC_varied_q.R` -- Code for the simulation setup for comparing the correctness of model selection by GMMFADq for data parameters n=300, p=10. 
    - `sim_gmfad_emmix_BIC_2.R`-- Code for the simulation setup for comparing the correctness of model selection by GMMFAD for data parameters n=300, p=10. 
    - `wdbc.data` -- Wisconsin breast cancer (Diagnostic) data with label.
    - `wdbc.names` -- Information related to the features of the Wisconsin breast cancer data.



## 🔗 How to Cite

If you reference any materials from this repository, please cite the main manuscript and optionally include the GitHub link in your supplementary section or appendix.

---

## 📬 Contact

For questions, suggestions, or collaboration inquiries, feel free to reach out:

**Kazeem Kareem**\
PhD Candidate, Statistics\
Michigan Technological University\
Email: kareemkazeem718\@gmail.com

---

