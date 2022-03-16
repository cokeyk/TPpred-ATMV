# TPpred-ATMV: Therapeutic peptide prediction by adaptive multi-view tensor learning model

## Requirements
This project is built on python, matlab and java. Therefore, you need to configure the above three environments before using it on win, linux or mac.
Note that the matlab must be compatible with the corresponding python verion, visit https://www.mathworks.com/content/dam/mathworks/mathworks-dot-com/support/sysreq/files/python-compatibility.pdf for more details.\
The main python libraries used in this project:

```
numpy
pandas
scipy
sklearn
matlab engine for python
```
## Downloads

You can download the models we trained and the tools used in this project by the following URL.\
https://drive.google.com/drive/folders/1pJAChvcrivxWZWnPaY98b0T9j1EqXqP6?usp=sharing. Your need to download both `tool` and `matlab_models`, and then put them into the same directory as `predict.py`.

## Predict
`predict.py` is used to performe predict task. For example:
```python
python predict.py -src test -fasta AAPT.txt -type AAP -th 0.5 -o result
```
## References
In this project we used PEPred-Suite[1] to generate 5 features: `Bit20NTCT(2)40 CTD3 GAP5 Ngram(N=1)1 OVNT(5)84`. We used BioSeq-Analysis2.0[2] to generate 6 features `DP DR DT Kmer PC-PseAAC-General Top-n-gram`. And we use our code to generate the rest features.

[1]L. Wei, C. Zhou, R. Su, and Q. Zou, PEPred-Suite: improved and robust prediction of therapeutic peptides using adaptive feature representation learning, Bioinformatics, vol. 35, no. 21, pp. 4272-4280, Nov 1 2019.

[2]Liu, B., et al. BioSeq-Analysis2.0: an updated platform for analyzing DNA, RNA, and protein sequences at sequence level and residue level based on machine learning approaches. Nucleic Acids Research 2019;47(20):e127.
