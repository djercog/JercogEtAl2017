# Model for cortical UP-DOWN states in a bistable inhibitory-stabilized network (Jercog et al 2017)

This program is free software: you can redistribute it and/or modify it under the terms of the GNU General Public License as published by the Free Software Foundation, either version 3 of the License, or (at your option) any later version. This program is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License for more details (http://www.gnu.org/licenses/).

Contains both rate and spiking network model implementations from the paper Jercog et al 2017. Both models generate the time courses for E & I population rates and adaptation presented. In addition, the Spiking network generates individual spiking activity in rasters. A complete description of the model is provided in the methods of the paper.

Rate model (Matlab 2009 - Linux).<br>
  Include following files in a separated folder:<br>
  singleSimEIModel.m: main function (also plots the output traces).<br>
  eiModelParams.m: model parameters.<br>
  eiModel.mexa64: Compiled MEX-function from source eiModel.cpp (Unix).<br>
  eiModel.cpp: c++ mex-function source.<br>
  rate_readme.txt: description.<br>
<br>
Spiking Model (g++ 6.4.3 - Linux):<br>
  Include following files in a separated folder:<br>
  main.cpp: main function, generate output bin files.<br>
  functions.h:   auxiliar functions.<br>
  functions.cpp: auxiliar functions.<br>
  spiking_readme.txt: description.<br>
<br>
Further instructions for usage of each model in their corresponding readme-files.

<p>

<a href="https://elifesciences.org/articles/22425">
"UP-DOWN cortical dynamics reflect state transitions in a bistable network"<br>
Jercog D, Roxin A, Barthó P, Luczak A, Compte A & de la Rocha J<br>
Elife 2017; DOI: 10.7554/eLife.22425
</a>
