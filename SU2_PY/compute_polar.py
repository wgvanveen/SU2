#!/usr/bin/env python

## \file shape_optimization.py
#  \brief Python script for performing the shape optimization.
#  \author Francisco Palacios, Trent Lukaczyk, Aerospace Design Laboratory (Stanford University) <http://su2.stanford.edu>.
#  \version 3.2.3 "eagle"
#
# SU2, Copyright (C) 2012-2013 Aerospace Design Laboratory (ADL).
#
# SU2 is free software; you can redistribute it and/or
# modify it under the terms of the GNU Lesser General Public
# License as published by the Free Software Foundation; either
# version 2.1 of the License, or (at your option) any later version.
#
# SU2 is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the GNU
# Lesser General Public License for more details.
#
# You should have received a copy of the GNU Lesser General Public
# License along with SU2. If not, see <http://www.gnu.org/licenses/>.

# imports
import numpy as np
from optparse import OptionParser
import os, sys, shutil, copy
sys.path.append(os.environ['SU2_RUN'])
import SU2

# Command Line Options
parser = OptionParser()
parser.add_option("-f", "--file", dest="filename",
                  help="read config from FILE", metavar="FILE")
parser.add_option("-n", "--partitions", dest="partitions", default=2,
                  help="number of PARTITIONS", metavar="PARTITIONS")
parser.add_option("-p", "--oldpartitions", dest="oldpartitions", default="oldpartitions",
                  help="old number of PARTITIONS (use -n instead)", metavar="OLDPARTITIONS")
parser.add_option("-i", "--iterations", dest="iterations", default=99999,
                  help="number of ITERATIONS", metavar="ITERATIONS")

(options, args)=parser.parse_args()
options.partitions = int( options.partitions )
options.iterations = int( options.iterations )

if options.oldpartitions != "oldpartitions":
  print ("\n IMPORTANT: -p is no longer available in SU2 v3.2.3, use -n flag instead \n")
  sys.exit()

# load config, start state
config = SU2.io.Config(options.filename)
state  = SU2.io.State()

# prepare config
config.NUMBER_PART = options.partitions
config.EXT_ITER    = options.iterations

# find solution files if they exist
state.find_files(config)

# angles to run
angles = np.linspace(0.,14.,8)

# mach numbers to run
mach = np.linspace(0.15,0.95,10)

# start results data
results = SU2.util.bunch()
results.AoA = angles
results.MACH = mach
results.DRAG = []
results.LIFT = []
results.MOMENT_Z = []

# iterate mach
for MachNumber in mach:

  f = open('Polar_M' + str(MachNumber) + '.dat', 'w')
  f.write('VARIABLES = "AoA", "C<sub>L</sub>", "C<sub>D</sub>", "C<sub>Mz</sub>" \n')
  
  # iterate angles
  for AngleAttack in angles:
    
    # local config and state
    konfig = copy.deepcopy(config)
    ztate  = copy.deepcopy(state)
    
    # set angle of attack
    konfig.AoA = AngleAttack
    konfig.MACH_NUMBER = MachNumber
    print 'Mach = ' , konfig.MACH_NUMBER , 'AoA = ' , konfig.AoA
    
    # run su2
    drag = SU2.eval.func('DRAG',konfig,ztate)
    lift = SU2.eval.func('LIFT',konfig,ztate)
    moment = SU2.eval.func('MOMENT_Z',konfig,ztate)

    # append results
    results.DRAG.append(drag)
    results.LIFT.append(lift)
    results.MOMENT_Z.append(moment)

    output = str(AngleAttack) + ", " + str(lift) + ", " + str(drag) + ", " + str(moment) + "\n"

    f.write(output)

  f.close()

# Close open file

#: for each angle

# plotting
# plt.figure()
# plt.plot( results.MACH_NUMBER, results.AoA , results.LIFT , results.DRAG )
# plt.show()

# save data
SU2.io.save_data('results.pkl',results)
