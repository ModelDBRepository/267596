Here are instructions to generate figures of the Muller, Abbott and Sawtell (2023), Current Biology paper titled: "A Mechanism for Differential Control of Axonal and Dendritic Spiking Underlying Learning in a Cerebellum-like Circuit" 

1) Install the NEURON module (https://www.neuron.yale.edu/neuron/download).

2) Download all the associated .mod files into a folder and compile the folder with mknrndll (drop the folder onto mknrndll in the NEURON folder).

3) Write in the directory for the folder just generated in the second line of the main code files. (neuron.load_mechanisms('/Users/..../[name of folder containing the compiled files]')).

4) Download '.ipynb' files and 'MG_constructed_by_Nate_Sawtell' files and make sure they are in the same folder.

5)Run the '.ipynb' files (Python notebook file). Follow the instructions for specific plots. Simulation may take a few minutes.

6) To generate the 2 compartment model (figure 4) run the .m file.