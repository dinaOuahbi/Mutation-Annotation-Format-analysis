# **Pancreatic Cancer from TCGA**  
**Analysis of Mutational Data Based on a Deep Learning Survival Model**  

## **Overview**  

This project explores the mutational landscape of pancreatic cancer patients from TCGA, following a survival analysis using a **deep learning model based on tile clusters** (Multiple Instance Learning).  

- A **linear predictor** was derived to summarize the contribution of tile clusters in predicting patient survival.  
- The predictor was **dichotomized using Maxstat** ([reference](http://www.sthda.com/french/wiki/maxstat-et-courbe-de-survie-pour-une-variable-continue-avec-rquery)).  
- Based on this threshold, patients were categorized into two groups:  
  - **HIGH** → Good prognosis  
  - **LOW** → Poor prognosis  

## **Objective**  

We aim to analyze the **mutational data** of these two patient groups by:  
- Identifying mutation types.  
- Assessing mutational load per patient and per group.

![Image of aciduino on protoboard](https://github.com/dinaOuahbi/Mutation-Annotation-Format-analysis/blob/main/oncoplot.png)
