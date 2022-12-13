## SRX-DTI

Predicting drug-target interactions based on fusing multiple features with data balancing and feature selection techniques
Predicting drug-target interaction (DTI) is an important research area in the field of drug discovery.  This framework proposes a novel drug–target interaction prediction method called SRX-DTI. First, we extract various descriptors from the protein sequences; and the drug is encoded as an FP2 molecular fingerprint. For handling the class imbalance problem, we propose the One-SVM-US technique to deal with imbalanced data. We also develop the FFS-RF algorithm, a forward feature selection algorithm, and coupled it with the random forest (RF) classifier to maximize the predictive performance. The forward feature selection algorithm adds new features to a set of selected features as far as the predictive power is improved. This feature selection algorithm removes the irrelevant features to obtain the best optimal features. Finally, the balanced dataset with optimal features is given to the XGBoost classifier to identify DTIs. 

![image](https://user-images.githubusercontent.com/72028345/204578716-30f41a3e-0f22-4881-82dc-f0af97e1eb52.png)

## About data
In this research, four golden standard datasets, including enzymes (EN), G-protein-coupled receptors (GPCR), ion channel (IC), and nuclear receptors (NR) released by [Yamanishi et al](http://web.kuicr.kyoto-u.ac.jp/supp/yoshi/drugtarget/).  All these datasets are freely available from [http://web.kuicr.kyoto-u.ac.jp/supp/yoshi/drugtarget/](http://web.kuicr.kyoto-u.ac.jp/supp/yoshi/drugtarget/). 

## Environment Settings
- Python version:  '3.9' or higher

- You have to install the required libraries

## To run the code
- Run ./feature extraction/00-AAC.py: extract AAC descriptor (for other descriptors, just change to related python code).  
- Run ./NR-run/run.py: make balanced dataset, and feature selection.
- Run ./DTI prediction/DTI_predictor.py: predict drug-target interactions, and evaluate the results with five cross-validation.

# Citation
If you use our source code, dataset, and experiments for your research or development, please cite our paper:
Khojasteh, H., Pirgazi, J. (2022). Improving prediction of drug-target interactions based on fusing multiple features with data balancing and feature selection techniques. bioRxiv 2022.12.07.519302.
https://www.biorxiv.org/content/10.1101/2022.12.07.519302v2

# Contact
If you have any questions, do not hesitate to contact me by `khojasteh@znu.ac.ir` or `khojasteh.hb@gmail.com`, I will be happy to assist.
