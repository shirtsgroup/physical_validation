import physical_validation as pv
import os

systems = ['none', 'shift', 'switch']
n_dts = 6

# change this to fit to your GROMACS installation
parser = pv.data.GromacsParser(exe='~/bin/gromacs/bin/gmx_d',
                               includepath='~/bin/gromacs/share/gromacs/top')

for sys in systems:
    print('### Analyzing system ' + sys)
    print('## Reading results')
    res = []
    # base simulation
    dir = os.path.join(sys, 'base')
    res.append(parser.get_simulation_data(
        mdp=os.path.join(dir, 'mdout.mdp'),
        top=os.path.join(dir, 'system.top'),
        gro=os.path.join(dir, 'system.gro'),
        edr=os.path.join(dir, 'system.edr')
    ))

    for n in range(1, n_dts-1):
        dir = os.path.join(sys, 'integrator_' + str(n))
        res.append(parser.get_simulation_data(
            mdp=os.path.join(dir, 'mdout.mdp'),
            top=os.path.join(dir, 'system.top'),
            gro=os.path.join(dir, 'system.gro'),
            edr=os.path.join(dir, 'system.edr')
        ))

    # make plot directory
    if not os.path.exists('ana_argon_plots'):
        os.makedirs('ana_argon_plots')
    sysplot = os.path.join('ana_argon_plots', sys)
    
    print('## Validating integrator convergence')
    pv.integrator.convergence(res, verbose=True,
                              filename=sysplot)
    print()
