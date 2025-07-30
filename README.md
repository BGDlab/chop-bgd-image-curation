# chop-bgd-image-curation
Customized HeuDiConv and CuBIDS configurations for the curation of brain imaging data from CHOP for the BGD lab

## HeuDiConv
HeuDiConv provides a framework to convert raw DICOM images to NIFTI using a custom heuristic. `heuristic.py` was created to specifically convert CHOP clinical brain imaging data to the BIDS format. Dictionaries were created specific to data from CHOP to standardize acquisition (`acq_dict.json`, `label_dict.json`) and modality (`mod_dict.json`) names such as the protocol name 'routine_brain/cor_mts_post' corresponding to a T1w scan. Further, the heuristic aims to identify and label a standardize MPRAGE sequence that was harmonized across scanner sites in the CHOP network.



## CuBIDS
CuBIDS provides a framework to address the curation of highly hetereogenous data using a custom heuristic. `config.yml` was created to identify scans with similar sequence parameters allowing for specified deviation from common parameter values. 

`parameter_ranges.tsv` CHOP common sequences and parameters