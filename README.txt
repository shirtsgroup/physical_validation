Examples:

Analysis of gromacs files for vrescale thermostat
python analyze-md.py -d example_data -f vrescale.argon.down.xvg vrescale.argon.up.xvg -t 132.915071475571 137.138128524429 -b 200 -i 40 -v -g vrescale -e total -s 12345 > vrescale.txt
Produces figures vrescale_linear.pdf and vrescale_nonlinear.pdf

Analysis of gromacs files for berendsen thermostat
python analyze-md.py -d example_data -f berendsen.argon.down.xvg berendsen.argon.up.xvg -t 132.915071475571 137.138128524429 -b 200 -i 40 -v -g berendsen -e total -s 12345 > berendsen.txt
Produces figures berendsen_linear.pdf and berendsen_nonlinear.pdf

Analysis of gromacs files for Parrinell-Rahman barostat - tests enthaply distribution
python analyze-md.py -d example_data -f prH.argon.down.xvg prH.argon.up.xvg -t 121.43116059 128.56883941 -p 90 90 -b 200 -i 40 -v -g enthalpy -e enthalpy -s 12345 > enthalpy.txt
Produces figures enthalpy_linear.pdf and enthalpy_nonlinear.pdf

Analysis of gromacs files for Parrinell-Rahman barostat - tests volume distribution
python analyze-md.py -d example_data -f prV.argon.down.xvg prV.argon.up.xvg -t 125 125 -p 30 150 -b 200 -i 40 -v -g volume -e volume -s 12345 > volume.txt
Produces figures volume_linear.pdf and volume_nonlinear.pdf

Analysis of gromacs files for Parrinell-Rahman barostat - test joint energy/volume distribution
python analyze-md.py -d example_data -f prEV.argon.down.xvg prEV.argon.up.xvg -t 121.43116059 128.56883941  -p 30 150 -b 200 -i 40 -v -g energyvolume -e jointEV -s 12345 > energyvolume.txt
Produces figures energyvolume_linear.pdf and energyvolume_nonlinear.pdf

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
