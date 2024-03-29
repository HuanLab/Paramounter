Note

1. mzdiff is used as the mass tolerance to dereplicate the features (similar m/z values and retention times) extracted by XCMS CentWave. Suggested default value: 0.001 or 0.01. However, some true positive metabolic features with mass differences smaller than that value may be removed by mistake. Therefore, if a user wants to disable the dereplication function, set the mzdiff to be any negative value.
2. Paramounter tunes an optimized peak height to maximize the number of true positive features. A drawback of that optimized value is the higher rate of false positive features and the likelihood of software crash. To address this issue, users can try a higher peak height threshold to reduce the number of false positive features (e.g., 2X the optimized peak height threshold).
