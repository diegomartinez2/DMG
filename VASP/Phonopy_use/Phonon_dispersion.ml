There are two ways to calculate the phonon frequencies. Here, I discussed about Density Functional Perturbation Theory (DFPT) using Phonopy code. It's already available in the website
https://atztogo.github.io/phonopy/vasp-dfpt.html#vasp-dfpt-interface
Please make sure about that U installed the Phonopy code in clean. Create one folder and make the folder path in Ur terminal. Keep Ur POSCAR in the folder by name of POSCAR-unitcell
Step 1
First we are in need to create a super cell of our system. Lattice constant of the super cell will be around 1-1.5 nm minimum. So create super cell like 3x3x3 or 6x6x6 depends upon Ur system using the comment
phonopy -d --dim="2 2 2" -c POSCAR-unitcell
(222 is a supercell size, POSCAR-unitcell - U want to rename the optimized POSCAR as POSCAR-unitcell (because it's a POSCAR of unitcell)
Now one new SPOSCAR file was generated in the folder. Rename the SPOSCAR to POSCAR.
Step 2
Take the POSCAR create POTCAR corresponds to it. If need reduce the KPOINTS. Run VASP calculation using the following INCAR
PREC = Accurate ENCUT = 500 IBRION = 8 EDIFF = 1.0e-08 IALGO = 38 ISMEAR = 0; SIGMA = 0.1 LREAL = .FALSE. ADDGRID = .TRUE. LWAVE = .FALSE. LCHARG = .FALSE.
Step 3
After the calculation successfully ended copy the vasprun.xml file in the same folder, where U keep all files.
Create the force constant by (Before that make the folder path in Ur terminal)
phonopy --fc vasprun.xml
Step 4
Copy the band.conf file which is available in the Phonopy-Examples-Phonon folder (or search band.conf in folder of Phonopy where U installed in Ur computer)
Check the line FORCE_CONSTANTS = READ in band.conf is enabled or otherwise U just included the line in band.conf. Now run the comment
phonopy --dim="2 2 2" -c POSCAR-unitcell band.conf.
U will get the Phonon spectrum.

---

I just want to add few points.
If you want to use DFPT:
1) Before starting using Phonopy, please make sure the you have optimized your system very accurately. For instance, add these two flags to your optimization INCAR: EDIFF = 1.0e-08, EDIFFG = -1.0e-8
2) make sure that the ENCUT in the INCAR (both OPT and Phonopy) is at least 30% higher than ENMAX for each element in the POTCAR.
3) In Phonopy you have just ONE iteration, it means there is no optimization, Phonopy is just calculating the force on each atom. So, add NSW=1 to your Phonopy INCAR.
4) Use IBRION = 7 or 8 in your Phonopy's INCAR (see the VASP manual).
5) Adjust Kpoints. It means if you are making a supercell and increasing the size of the unit cell by 2*2*2, you need less kpoints than the kpoints that you have used in the optimization step (reciprocal lattice concept)
6) if you need thermodynamic properties, you can use:
phonopy -s -t -c POSCAR-unitcell thermo.conf > THERMO.data
7) At the beginning, when you use the below command to creat the supercell:
phonopy -d --dim="2 2 2" -c POSCAR-unitcell (you can change the size of the supercell, I just put "2 2 2" as an example)
you will get bunch of files like ...-00* . If you want to use DFPT, you caln delete these files. They can be used for "Frozen phonon" method.
8) To review the DFPT theory you can follow Stefano Baroni's papers in this matter.
