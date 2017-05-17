Contributed by Daniel Jercog.

This model code generates the data to obtain r_E and r_I population rates (countE.bin and countI.bin files, respectively), adaptation current I_a (Ada.bin file), timestamp of kicks impinged to the network (kicks.bin), spike count used for U and D period detection (countdet.bin), E and I cells spike timestamps for rastergram generation (rasterE.bin and rasterI.bin) from spiking network simulations.

The model is a network composed of leaky integrate-and-fire all-to-all connected neurons. Adaptation affects the E population. The network is in a bistable regime where external/internal fluctuations can cause transitions between two attractor states, and adaptation modulates the probability of those fluctuations to cause transitions.


___________________________

To compile:
  $g++ *.cpp -O3 -lm -o execSim

To execute the code model, 2 additional input parameters are required:
  - To use Fig 7 parameters set value to 1 (set to 0 for Fig 7 Supp 2 parameters)
  - To save rasters, second parameter in the call must be set to 1 (0 will not generate raster files)

Example call:
  $./execSim 1 1

___________________________

Output files:

Some model parameters used in are part of the [prefix] of output filenames. 
Thus, generated output files are:

[prefix]countE.bin: Vector of int32 of E cells spk count (bin size defined in countBin, in ms).
[prefix]countI.bin: Vector of int32 of I cells spk count (bin size defined in countBin, in ms).
[prefix]Ada.bin: Vector of int32 of 1000*adaptation current. 
[prefix]kicks.bin: Vector of double composed by pairs of values (timestamp of kick in units of dt, kick amplitude).
[prefix]countdet.bin: Vector of int32 of spk count for U and D detection (bin size defined in countBin, in ms).

In addition, if rasters are saved:
[prefix]rasterE.bin: Vector of int composed by pairs of values (timestamp of spk, id of cell).
[prefix]rasterI.bin: Vector of int composed by pairs of values (timestamp of spk, id of cell).


Depending on your CPU, the execution can take several minutes.
