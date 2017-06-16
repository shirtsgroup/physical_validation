###########################################################################
#                                                                         #
#    physical_validation,                                                 #
#    a python package to test the physical validity of MD results         #
#                                                                         #
#    Written by Michael R. Shirts <michael.shirts@colorado.edu>           #
#               Pascal T. Merz <pascal.merz@colorado.edu>                 #
#                                                                         #
#    Copyright (C) 2012 University of Virginia                            #
#              (C) 2017 University of Colorado Boulder                    #
#                                                                         #
#    This library is free software; you can redistribute it and/or        #
#    modify it under the terms of the GNU Lesser General Public           #
#    License as published by the Free Software Foundation; either         #
#    version 2.1 of the License, or (at your option) any later version.   #
#                                                                         #
#    This library is distributed in the hope that it will be useful,      #
#    but WITHOUT ANY WARRANTY; without even the implied warranty of       #
#    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU    #
#    Lesser General Public License for more details.                      #
#                                                                         #
#    You should have received a copy of the GNU Lesser General Public     #
#    License along with this library; if not, write to the                #
#    Free Software Foundation, Inc., 51 Franklin Street, Fifth Floor,     #
#    Boston, MA 02110-1301 USA                                            #
#                                                                         #
###########################################################################
r"""
GROMACS python interface.

.. warning:: This is a mere place holder, as an official python API is 
   currently being developed by the gromacs development team. It is 
   probably neither especially elegant nor especially safe. Use of this
   module in any remotely critical application is strongly discouraged.   
"""
import os
import sys
import subprocess
import re
import numpy as np


class GromacsInterface(object):
    def __init__(self, exe=None, dp=None):

        self._exe = None
        self._dp = False

        if dp is not None:
            self.dp = dp

        if exe is None:
            # check whether 'gmx' / 'gmx_d' is in the path
            if self.dp:
                if self._check_exe(quiet=True, exe='gmx_d'):
                    self.exe = 'gmx_d'
                else:
                    print('WARNING: gmx executable not found. Set before attempting to run!')
            else:
                if self._check_exe(quiet=True, exe='gmx'):
                    self.exe = 'gmx'
                else:
                    print('WARNING: gmx executable not found. Set before attempting to run!')
        else:
            self.exe = exe

    @property
    def exe(self):
        """exe is a string pointing to the gmx executable."""
        return self._exe

    @exe.setter
    def exe(self, exe):
        if self._check_exe(exe=exe):
            if os.path.dirname(exe):
                exe = os.path.abspath(exe)
            self._exe = exe

    @property
    def double(self):
        """double is a bool defining whether the simulation was ran at double precision"""
        return self._dp

    @double.setter
    def double(self, dp):
        assert isinstance(dp, bool)
        self._dp = dp

    def get_quantities(self, edr, quantities, cwd=None,
                       begin=None, end=None, args=None):

        if args is None:
            args = []

        tmp_xvg = 'gmxpy_' + os.path.basename(edr).replace('.edr', '') + '.xvg'
        if cwd is not None:
            tmp_xvg = os.path.join(cwd, tmp_xvg)

        q_dict = {}

        for q in quantities:
            not_found = self._create_xvg(edr, tmp_xvg, [q], cwd=cwd,
                                         begin=begin, end=end, args=args)[1]
            if q in not_found:
                q_dict[q] = None
                continue

            skip_line = re.compile("^[#,@]")
            values = []
            times = []
            with open(tmp_xvg, 'r') as xvg:
                for line in xvg:
                    if skip_line.match(line):
                        continue
                    times.append(float(line.split()[0]))
                    values.append(float(line.split()[1]))

            if 'time' in q_dict:
                if not np.array_equal(np.array(times), q_dict['time']):
                    print('WARNING: Time discrepancy in ' + edr)
            else:
                q_dict['time'] = np.array(times)
            q_dict[q] = np.array(values)

            os.remove(tmp_xvg)

        return q_dict

    def read_trr(self, trr):
        tmp_dump = 'gmxpy_' + os.path.basename(trr).replace('.trr', '') + '.dump'
        with open(tmp_dump, 'w') as dump_file:
            self._run('dump', ['-f', trr], stdout=dump_file)

        position = []
        velocity = []
        force = []
        box = []
        with open(tmp_dump) as dump:
            x = []
            v = []
            f = []
            b = []
            for line in dump:
                if 'frame' in line:
                    # new frame
                    if len(x) > 0:
                        # not the first frame - nothing to save there
                        position.append(np.array(x))
                        velocity.append(np.array(v))
                        force.append(np.array(f))
                        box.append(np.array(b))
                    x = []
                    v = []
                    f = []
                    b = []
                    continue
                if 'x[' in line:
                    x.append([float(l.strip()) for l in
                              line.split('{', 1)[1].split('}')[0].split(',')])
                if 'v[' in line:
                    v.append([float(l.strip()) for l in
                              line.split('{', 1)[1].split('}')[0].split(',')])
                if 'f[' in line:
                    f.append([float(l.strip()) for l in
                              line.split('{', 1)[1].split('}')[0].split(',')])
                if 'box[' in line:
                    b.append([float(l.strip()) for l in
                              line.split('{', 1)[1].split('}')[0].split(',')])
            # end loop over file - save last arrays
            position.append(np.array(x))
            velocity.append(np.array(v))
            force.append(np.array(f))
            box.append(np.array(b))

        result = {}
        for key, vector in zip(['position', 'velocity', 'force', 'box'],
                               [position, velocity, force, box]):
            vector = np.array(vector)
            if vector.size > 0:
                result[key] = vector
            else:
                result[key] = None
        return result

    @staticmethod
    def read_gro(gro):
        position = []
        velocity = []
        box = []
        with open(gro) as conf:
            x = []
            v = []
            b = []
            title = conf.readline()
            while title:
                natoms = int(conf.readline().strip())
                for n in range(natoms):
                    line = conf.readline()
                    line = line.split()
                    x.append([float(x) for x in line[3:6]])
                    v.append([float(v) for v in line[6:9]])

                line = conf.readline()
                line = line.split()
                b.append([float(v) for v in line[0:3]])
                title = conf.readline()

        result = {}
        for key, vector in zip(['position', 'velocity', 'force', 'box'],
                               [position, velocity, [], box]):
            vector = np.array(vector)
            if vector.size > 0:
                result[key] = vector
            else:
                result[key] = None
        return result

    def _check_exe(self, quiet=False, exe=None):
        # could we also check that it is actually Gromacs, not just any existing executable?
        if exe is None:
            exe = self._exe
        try:
            devnull = open(os.devnull)
            subprocess.call([exe], stdout=devnull, stderr=devnull)
        except OSError as e:
            if e.errno == os.errno.ENOENT:
                # file not found error.
                if not quiet:
                    print('ERROR: gmx executable not found')
                    print(exe)
                return False
        return True

    def _run(self, cmd, args, cwd=None, stdin=None, stdout=None, stderr=None):
        if self.exe is None:
            print('ERROR: No gmx executable defined. Set before attempting to run!')
        command = [self.exe, cmd]
        command.extend(args)
        return subprocess.Popen(command, cwd=cwd,
                                stdin=stdin, stdout=stdout, stderr=stderr)

    def _create_xvg(self, edr, xvg, quantities, cwd=None,
                    begin=None, end=None, args=None):
        assert os.path.exists(edr)
        assert os.path.exists(os.path.abspath(os.path.dirname(xvg)))

        if args is None:
            args = []

        if self._dp:
            args.append('-dp')
        if begin is not None:
            args.extend(['-b', str(begin)])
        if end is not None:
            args.extend(['-e', str(end)])

        quants = ''
        for q in quantities:
            quants += str(q) + '\n'

        args = ['-f', edr, '-o', xvg] + args
        proc = self._run('energy', args, cwd=cwd,
                         stdin=subprocess.PIPE, stdout=subprocess.PIPE, stderr=subprocess.PIPE)

        err = proc.communicate(quants.encode(sys.stdin.encoding))[1]

        encoding = sys.stderr.encoding
        if encoding is None:
            encoding = 'UTF-8'

        not_found = []
        if 'does not match anything' in err.decode(encoding):
            for q in quantities:
                if "String '" + q + "' does not match anything" in err.decode(encoding):
                    not_found.append(q)

        return proc.wait(), not_found
