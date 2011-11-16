Examples:

Analysis of gromacs files for vrescale thermostat
python analyze-gromacs.py -d example_data -f vrescale.argon.down.xvg vrescale.argon.up.xvg -t 132.915071475571 137.138128524429 -b 200 -i 40 -v -g vrescale -e total > vrescale.txt
produces figures vrescale_linear.pdf and vrescale_nonlinear.pdf

Analysis of gromacs files for berendsen thermostat
python analyze-gromacs.py -d example_data -f berendsen.argon.down.xvg berendsen.argon.up.xvg -t 132.915071475571 137.138128524429 -b 200 -i 40 -v -g berendsen -e total > berendsen.txt
produces figures berendsen_linear.pdf and berendsen_nonlinear.pdf

Analysis of harmonic oscillators:
python harmonic.py -t 0.5 1.5 -N 100000 100000 -g harmonic_energy > harmonic_energy.txt

harmonic oscillators with some noise, creating deviations from the fit.
python harmonic.py -t 0.5 1.5 -N 100000 100000 -o 0.1 -g harmonic_energy_noise > harmonic_energy_noise.txt

Analysis of harmonic oscillators with work functions:
python harmonic_pressure.py -t 0.5 1.5 -p 1.0 1.0 -e enthalpy -g harmonic_enthalpy > harmonic_enthalpy.txt
python harmonic_pressure.py -t 1.0 1.0 -p 0.7 1.3 -e volume -g harmonic_volume > harmonic_volume.txt
python harmonic_pressure.py -t 1.7 1.3 -p 0.8 1.2 -e jointEV > harmonic_jointEV.txt
