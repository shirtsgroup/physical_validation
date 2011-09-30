Call:

python analyze-gromacs.py -d example_data -f vrescale.argon.down.xvg vrescale.argon.up.xvg -t 132.915071475571 137.138128524429 -b 100 -i 40 -v -p vrescale -e total > vrescale.txt

python analyze-gromacs.py -d example_data -f berendsen.argon.down.xvg berendsen.argon.up.xvg -t 132.915071475571 137.138128524429 -b 100 -i 40 -v -p berendsen -e total > berendsen.txt
