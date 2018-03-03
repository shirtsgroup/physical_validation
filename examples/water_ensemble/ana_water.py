import physical_validation as pv
import os

systems = ['vr', 'be', 'vr_pr', 'be_pr']

# change this to fit to your GROMACS installation
parser = pv.data.GromacsParser(exe='~/bin/gromacs/bin/gmx',
                               includepath='~/bin/gromacs/share/gromacs/top')

for sys in systems:
    print('### Analyzing system ' + sys)
    print('## Reading lower temperature result')
    dir_low = os.path.join(sys, 'base')
    res_low = parser.get_simulation_data(
        mdp=os.path.join(dir_low, 'mdout.mdp'),
        top=os.path.join(dir_low, 'system.top'),
        gro=os.path.join(dir_low, 'system.gro'),
        edr=os.path.join(dir_low, 'system.edr')
    )
    print('## Reading high temperature result')
    dir_high = os.path.join(sys, 'ensemble_1')
    res_high = parser.get_simulation_data(
        mdp=os.path.join(dir_high, 'mdout.mdp'),
        top=os.path.join(dir_high, 'system.top'),
        gro=os.path.join(dir_high, 'system.gro'),
        edr=os.path.join(dir_high, 'system.edr')
    )

    if not os.path.exists('ana_water_plots'):
        os.makedirs('ana_water_plots')
    sysplot = os.path.join('ana_water_plots', sys)
    
    print('\n## Validating kinetic energy distribution (strict)')
    print('# Low T:')
    pv.kinetic_energy.distribution(res_low, verbosity=2, strict=True,
                                   filename=sysplot + '_low_mb')
    print('# High T:')
    pv.kinetic_energy.distribution(res_high, verbosity=2, strict=True,
                                   filename=sysplot + '_high_mb')

    print('\n## Validating kinetic energy distribution (non-strict)')
    print('# Low T:')
    pv.kinetic_energy.distribution(res_low, verbosity=2, strict=False)
    print('# High T:')
    pv.kinetic_energy.distribution(res_high, verbosity=2, strict=False)
    
    print('\n## Validating ensemble')
    if 'pr' in sys:
        # can't plot in 2d...
        quantiles = pv.ensemble.check(res_low, res_high, verbosity=2)
    else:
        quantiles = pv.ensemble.check(res_low, res_high, verbosity=2,
                                      filename=sysplot + '_ensemble')
    if len(quantiles) == 1:
        q_str = '{:.1f}'.format(quantiles[0])
    else:
        q_str = '('
        for q in quantiles:
            q_str += '{:.1f}, '.format(q)
        q_str = q_str[:-2] + ')'
    print('Calculated slope is ' + q_str +
          ' quantiles from the true slope')
    print('\n')
