Simple demonstration of a molecular dynamics simulation

408 BMPyrDCA ionic pairs between two graphene walls

All files are prepared for the simulation. The force fields are stored at github: github.com/vladislavivanistsev/RTIL-FF. References are given within the files.


Let's execute some commands to tell the computer how to run the simulation. However, don't forget to allow the computer to read and execute all files. To do so, in the terminal enter "chmod 755 *".

First we pack 408 BMImDCA ionic pairs into a predefined box.

./packmol < packmol.inp

gmx editconf -f packmol.pdb -o packmol.gro

In the packmol.gro set the box length to 7.0 nm.

Now let's creare an index of the system. Define the groups.

Now let's create an index of the ionic pairs. Enter q.

Prepare an executable file for the simulation.

gmx grompp -f STEEP.mdp -c packmol.gro -p topol.top -n index.ndx -o STEEP

Execute the first simulation step.

gmx mdrun -deffnm STEEP

Prepare an other executable file for the simulation.

gmx grompp -f RUN.mdp -c STEEP.gro -p topol.top -n index.ndx -o NVT

Run the simulation, so called production run. Note, it should use all available CPUs. Caution! The simulation requires a lot of time.

gmx mdrun -deffnm NVT

Check the potential, field, and charge profiles. Plot the .xvg-resulting files using xmgrace (like xmgrace potential.xvg).

gmx potential -f NVT.trr -n index.ndx -s RUN.tpr -sl 1000 -o potential_0.xvg -oc charge_0.xvg -of field_0.xvg

Let's apply field to induce the ionic layering at the modelled interface.

gmx grompp -f FIELD.mdp -c NVT.gro -p topol.top -n index.ndx -o FIELD

Run the simulation with the field. Caution! The simulation requires a lot of time.

gmx mdrun -deffnm FIELD

Check the potential, field, and charge profiles. Plot the .xvg-resulting files using xmgrace (like xmgrace potential.xvg). Compare them with the ones obtained for the simulation without any applied field.

gmx potential -f FIELD.trr -n index.ndx -s FIELD.tpr -sl 1000 -o potential_10.xvg -oc charge_10.xvg -of field_10.xvg
