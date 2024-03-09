# PIPI2
PIPI2: Sensitive Tag-based Database Search to Identify Peptides with Multiple PTMs

Contact: slaiad@connect.ust.hk
https://bioinformatics.hkust.edu.hk/
 
## Requirement
- Java 1.8.
Download and install JDK 8 from https://www.oracle.com/java/technologies/downloads/#java8

Usage:

with Java 8 (recommended):
java -Xmx8g -jar PIPI.jar parameter.def spectra_file output_directory

with Java of higher versions:
java -Xmx8g -jar --add-opens java.base/java.util=ALL-UNNAMED --add-opens java.base/java.lang=ALL-UNNAMED PIPI.jar parameter.def spectra_file output_directory

## A test data set
On Zenodo https://zenodo.org/records/10115308

## About
Peptide identification is important to protein inference in bottom-up proteomics. Post-translational modifications (PTMs) are crucial in regulating cellular activities. Many database search methods have been developed to identify peptides with PTMs and characterize the numbers, types, and sites of the PTMs. However, the existence of PTMs on peptides hinders the peptide identification rate and the PTM characterization precision, especially for peptides with multiple PTMs. To address this issue, we present a sensitive open search engine, PIPI2, to identify peptides with multiple PTMs and characterize the PTM patterns.



