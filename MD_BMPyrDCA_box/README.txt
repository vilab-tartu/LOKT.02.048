Simple demonstration of a molecular dynamics simulation

25 BMPyrDCA ionic pairs in a box

All files are prepared for the simulation. The force fields are stored at github: github.com/vladislavivanistsev/RTIL-FF. References are given within the files.


Let's execute some commands to tell the computer how to run the simulation. However, don't forget to allow the computer to read and execute all files. To do so, in the terminal enter "chmod 755 *".

First, we pack 50 BMImDCA into a predefined box. For packmol we should thank J. M. Martínez and L. Martínez. Visit m3g.iqm.unicamp.br/packmol/ to see details.

./packmol < packmol.inp

gmx editconf -f packmol.pdb -o packmol.gro

Now let's create an index of the ionic pairs. Enter q.

gmx make_ndx -f packmol.gro -o index.ndx

Prepare an executable file for the simulation.

gmx grompp -f STEEP.mdp -c packmol.gro -p topol.top -n index.ndx -o STEEP

Execute the first simulation step.

gmx mdrun -deffnm STEEP

Prepare another executable file for the simulation.

gmx grompp -f RUN.mdp -c STEEP.gro -p topol.top -n index.ndx -o NVT

Run the simulation, so-called production run. Note, it should use all available CPUs.

gmx mdrun -deffnm NVT

Convert the trajectory for the whole system. Enter 0.

gmx trjconv -f NVT.xtc -s NVT.tpr -o NVT.trr -ur compact -pbc mol

To see the equilibrated system, open a terminal and type "vmd NVT.gro". In the "VMD main" window select the single line "0 T A D F NVT.gro ...", then from File menu choose "load data into molecule", select "NVT.trr" file.

vmd NVT.gro

Examine the spatial distribution. For example select group 2 and then 3 to see the destribution of 3 (DCA anion) around 2 (Pyr cation).

gmx spatial -f NVT.xtc -s NVT.tpr -n index.ndx

Visualize the grid.cude in VMD. Chech the representations and drawing method: isosurface.

vmd grid.cube
