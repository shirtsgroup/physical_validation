Examples:

Analysis of gromacs files for vrescale thermostat
python analyze-gromacs.py -d example_data -f vrescale.argon.down.xvg vrescale.argon.up.xvg -t 132.915071475571 137.138128524429 -b 200 -i 40 -v -g vrescale -e total -s 12345 > vrescale.txt
Produces figures vrescale_linear.pdf and vrescale_nonlinear.pdf

Analysis of gromacs files for berendsen thermostat
python analyze-gromacs.py -d example_data -f berendsen.argon.down.xvg berendsen.argon.up.xvg -t 132.915071475571 137.138128524429 -b 200 -i 40 -v -g berendsen -e total -s 12345 > berendsen.txt
Produces figures berendsen_linear.pdf and berendsen_nonlinear.pdf

Analysis of harmonic oscillators:
python harmonic.py -t 0.5 1.5 -g harmonic_energy -s 12345 > harmonic_energy.txt
Produces figures harmonic_energy_linear.pdf and harmonic_energy_nonlinear.pdf

harmonic oscillators with some noise, creating deviations from the fit.
python harmonic.py -t 0.5 1.5 -o 0.1 -g harmonic_energy_noise > harmonic_energy_noise.txt
Produces figures harmonic_energy_noise_linear.pdf and harmonic_energy_noise_nonlinear.pdf

Analysis of harmonic oscillators with work functions:
python harmonic_pressure.py -t 0.5 1.5 -p 1.0 1.0 -e enthalpy -g harmonic_enthalpy -s 12345 > harmonic_enthalpy.txt
Produces figures harmonic_enthalpy_linear.pdf and harmonic_enthalpy_nonlinear.pdf

python harmonic_pressure.py -t 1.0 1.0 -p 0.7 1.3 -e volume -g harmonic_volume -s 12345 > harmonic_volume.txt
Produces figures harmonic_volume_linear.pdf and harmonic_volume_nonlinear.pdf

python harmonic_pressure.py -t 1.7 1.3 -p 0.8 1.2 -e jointEV -s 12345 > harmonic_jointEV.txt
Produces figures harmonic_jointEV_linear.pdf and harmonic_jointEV_nonlinear.pdf
