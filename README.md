# Proteolytic Module for _Bacteroides caccae_ Metabolic Model

This repository contains scripts and resources for developing and implementing a proteolytic module in a _Bacteroides caccae_ ATCC 43185 metabolic model.

## Project Description

The aim of this project is to enhance the _B. caccae_ metabolic model growth and metabolic production prediction accuracy by incorporating a proteolytic module. The proteolytic module allows for protein breakdown and provides a more comprehensive representation of the organism's metabolic capabilities.

## Dependency

- MATLAB: Make sure you have MATLAB installed on your system to run the scripts and perform the simulations.

- [CobraToolBox](https://github.com/opencobra/cobratoolbox.git): The CobraToolBox, a MATLAB-based toolbox, is required for working with the B. caccae metabolic model and performing FBA simulations.

- Homemade Scripts: In addition to the CobraToolBox, several homemade scripts were developed specifically for this project. These scripts modify the functionality of the CobraToolBox and provide additional custom features.

## Usage 

### Optimization: 

- Optimization: optimize the exchange reactions' lower and upper bounds using an in vitro training data set using a GLM and a P20 in silico media. The optimization allows to obtain a parametrized _B. caccae_ metabolic model.

```matlab
optimisation
```

Output: optimized _B. caccae_ metabolic model ready to be used for dynamicFBA 

Requirements: parameters files (list) 

### Perform dynamic Flux Balance Analysis: 

- Dynamic FBA Simulation: the modified dynamic FBA function incorporate additional regulations for improved growth and metabolic production prediction accuracy. dFBA analysis is performed on 3 different media (GLM, P2 and P20)

 - Comparison with laboratory data:

```matlab
Test_model_after_optimisation_3_media
```

Output: 

Requirements: parameters files (list)

## Licence 
The MIT License (MIT)

Copyright (c) 2023 Biomathematica

Permission is hereby granted, free of charge, to any person obtaining a copy of this software and associated documentation files (the "Software"), to deal in the Software without restriction, including without limitation the rights to use, copy, modify, merge, publish, distribute, sublicense, and/or sell copies of the Software, and to permit persons to whom the Software is furnished to do so, subject to the following conditions:

The above copyright notice and this permission notice shall be included in all copies or substantial portions of the Software.

## Citation 

A completer et lien vers le pre print éventuellement 

## Contact Information

For any questions or inquiries regarding this project, please contact [provide your contact information].
