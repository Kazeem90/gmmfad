# Supplementary Materials for "A hybrid mixture approach for clustering and characterizing biomedical data"

This repository contains source code for the paper\
Kareem KA, Dai F (2026) A hybrid mixture approach for clustering and characterizing biomedical data. PLoS ComputBiol 22(11): e1014974. [doi:10.1371/journal.pcbi.1014974](https://doi.org/10.1371/journal.pcbi.1014974)

---

## 📄 Contents

### `main`

- Contains the source code for the algorithms, simulation studies and applications in the paper.
    - `breast_cancer_data_constq_variedq.R`-- Code for fitting GMMFAD and GMMFAD-q on the breast cancer dataset.
    - `data_pca_kmeans_gmm.R`-- Code for PCA-based clustering baselines with data-driven PC selection.
    - `em4gmm.R` -- This is the primary source code, containing the functions GMMFAD  and other relevant functions.
    - `gmfad_emmix_emEM_p350_sim_time.R`-- Code for the simulation setup for comparing the time speedup, frobsenius errors, and ARI for GMMFAD and EMMIX algorithms for n=p=150. 
    - `gmfad_emmix_emEM_sim_time.R`-- Code for the simulation setup for comparing the time speedup, frobsenius errors, and ARI for GMMFAD and EMMIX algorithms for n=300, p=10. 
    - `gmfad_qq` -- This is the primary source code, containing the functions GMMFAD-q, and other relevant functions.
    - `human-brain-scRNA-seq-data-preprocessing.R`-- Code for constructing the compact sparse-omics benchmark dataset using two adult human brain cell classes.
    - `human-brain-scRNA-seq-data-variedq.R`-- Code for fitting the human brain scRNAseq data with GMMFAD-q.
    - `lymphoma_gene_data.R`--  Code for fitting the lymphoma data with GMMFAD-q.
    - `sim_BIC_varied_q.R` -- Code for the simulation setup for comparing the correctness of model selection by GMMFAD-q for data parameters n=300, p=10. 
    - `sim_gmfad_emmix_BIC.R`-- Code for the simulation setup for comparing the correctness of model selection by GMMFAD for data parameters n=300, p=10. 
    - `wdbc.data` -- Wisconsin breast cancer (Diagnostic) data with label.
    - `wdbc.names` -- Information related to the features of the Wisconsin breast cancer data.



## 🔗 How to Cite

If you reference any materials from this repository, please cite the main paper and optionally include the GitHub link in your supplementary section or appendix.

---

## 📬 Contact

For questions, suggestions, or collaboration inquiries, feel free to reach out:

**Kazeem Kareem**\
PhD, Statistics\
Michigan Technological University\
Email: kareemkazeem718\@gmail.com

---

