# Analyses for Exp 2 of Jones et al., Psych Science, 2024

## 1. Behavioral analyses (Table 1)
`behavioral analyses.ipynb`

## 2. IEM results (Figs 7,8)
`IEM_alpha_analyses.ipynb`

## 3. mvLoad Decoding (Fig 9)
`decode_load.ipynb` 

## 4. Hyperplane Distances (Fig 10)
`hyperplane.ipynb`

## 5. RSA (Figs 11, 12)
`RSA.ipynb`

## Controls for Eyemovements (Fig 13)
`IEM_alpha_analyses-eyeGazeBins`

the following notebooks include sections at the end excluding subjects with seemingly informative eye movements:  
- `IEM_alpha_analyses.ipynb`  
- `decode_load.ipynb`  
- `hyperplane.ipynb`
- `RSA.ipynb`


## Helper functions
`utils.py` includes functions used across the notebooks to clean up the behavioral data for easy analysis and condition labeling.  
`EEG_Decoder/` contains a copy of the package [https://github.com/AwhVogelLab/EEG_Decoder] that has been slightly modified to work with the condition structure of this experiment.