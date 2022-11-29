## SRX-DTI

Predicting drug-target interactions based on fusing multiple features with data balancing and feature selection techniques
Predicting drug-target interaction (DTI) is an important research area in the field of drug discovery.  This framework proposes a novel drug–target interaction prediction method called SRX-DTI. First, we extract various descriptors from the protein sequences; and the drug is encoded as an FP2 molecular fingerprint. For handling the class imbalance problem, we propose the One-SVM-US technique to deal with imbalanced data. We also develop the FFS-RF algorithm, a forward feature selection algorithm, and coupled it with the random forest (RF) classifier to maximize the predictive performance. The forward feature selection algorithm adds new features to a set of selected features as far as the predictive power is improved. This feature selection algorithm removes the irrelevant features to obtain the best optimal features. Finally, the balanced dataset with optimal features is given to the XGBoost classifier to identify DTIs. 

![image](https://user-images.githubusercontent.com/72028345/204578716-30f41a3e-0f22-4881-82dc-f0af97e1eb52.png)


## Environment Settings
- Python version:  '3.9' or higher

- You have to install the required libraries


## To run the code 

# Citation
If you use our source code, dataset, and experiments for your research or development, please cite the following paper:

# Contact
If you have any questions, do not hesitate to contact me by `khojasteh@znu.ac.ir` or `khojasteh.hb@gmail.com`, I will be happy to assist.
