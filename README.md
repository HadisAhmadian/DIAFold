# Enhancing Protein Structure Prediction via Physicochemical Feature-Based Homology Search

This project evaluates the impact of amino acid physicochemical features in optimizing homology search for protein tertiary structure prediction.

## Overview

We propose a modified homology search pipeline (`O-Physicho`) that leverages DIAMOND with reduced amino acid alphabets based on physicochemical properties. Integrated into a ColabFold-based prediction framework, the method achieves:

- **~5.9× faster** homology search
- **~37× fewer false positives**
- **Equal or improved** structural prediction accuracy (TM-Score, GDT, RMSD)

## Pipeline

1. **Homolog Search:** DIAMOND with feature-based prefiltering  
2. **MSA Construction:** Star alignment from DIAMOND pairwise output  
3. **Structure Prediction:** ColabFold using the generated MSA  

## Datasets

- **Homology Search Benchmark:** 1.7M SCOPe-annotated sequences vs. UniRef50  
- **Structure Evaluation:** 47 single-chain proteins from CASP14/15  

## Key Findings

- Reduced false positives without loss of homolog coverage  
- Improved structure prediction precision under lower sequence noise  
- Maintains or exceeds performance of MMseqs2-based pipelines

## Citation

If using this method, please cite the corresponding paper:

> *Enhancing Protein Structure Prediction: Evaluating the Role of Amino Acid Physicochemical Features in Homology Search* (2025)

