# FFCM-MRF: An accurate and generalizable cerebrovascular segmentation pipeline for humans and rhesus monkeys based on TOF-MRA
FFCM-MRF is an open-source MATLAB GUI for ToF-MRA data preprocessing, cerebrovascular segmentation and quantification. It is capable of accurately segmenting intracranial cerebrovascular structures, despite variations in scanner platforms, sequences, and species.

For detailed information using FFCM-MRF, please refer to the manual.pdf file.


Due to git LFS limitations on large file downloads, we deleted the model file that was originally located in “./HD-BET/hd-bet_params”, and modified it to automatically detect whether it exists when executing the de-cranialization command, and download it automatically if the file does not exist.

## :books: Citation
Please cite the [following paper](https://doi.org/10.1016/j.compbiomed.2024.107996) when using FFCM-MRF:

Yue Cui, Haibin Huang, Jialu Liu, Mingyang Zhao, Chengyi Li, et al. FFCM-MRF: An accurate and generalizable cerebrovascular segmentation pipeline for humans and rhesus monkeys based on TOF-MRA, Computers in Biology and Medicine, 170: 107996, 2024
