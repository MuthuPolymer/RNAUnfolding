#!/bin/sh
#SBATCH -p short
#SBATCH -J 120A
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --time=00-22
#SBATCH --output=energy.out

# Run your executable
echo Run submitted!
echo Time is `date`
./lmp_serial<efield.in> efield.out
echo Run completed!
echo Time is `date`

