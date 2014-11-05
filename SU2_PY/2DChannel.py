#!/usr/bin/env python 

## \file 2DChannel.py
#  \brief Python script for creating grid for freesurface channels.
#  \author F. Palacios
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

from optparse import OptionParser

parser=OptionParser()
parser.add_option("-f", "--file", dest="filename", default="Channel.su2",
                  help="write mesh to FILE", metavar="FILE")
parser.add_option("-n", "--nNode", dest="nNode", default=125,
                  help="use this NNODE in x direction", metavar="NNODE")
parser.add_option("-m", "--mNode", dest="mNode", default=100,
                  help="use this MNODE in y direction", metavar="MNODE")
parser.add_option("-x", "--xLength", dest="xLength", default=7.0,
                  help="use this XLENGTH", metavar="XLENGTH")
parser.add_option("--offsetx", dest="offsetx", default=1.0,
                  help="use this OFFSETX", metavar="OFFSETX")
parser.add_option("--xAdapt", dest="xAdapt", default="True",
                  help="Adapt the grid XADAPT", metavar="XADAPT")
parser.add_option("-y", "--yLength", dest="yLength", default=1.0,
                  help="use this YLENGTH", metavar="YLENGTH")
parser.add_option("--offsety", dest="offsety", default=0.5,
                  help="use this OFFSETY", metavar="OFFSETY")
parser.add_option("--yAdapt", dest="yAdapt", default="True",
                  help="Adapt the grid YADAPT", metavar="YADAPT")
(options, args)=parser.parse_args()

nNode = int(options.nNode)
mNode = int(options.mNode)
xLength = float(options.xLength)
yLength = float(options.yLength)
xAdapt = options.xAdapt
yAdapt = options.yAdapt

Mesh_File = open(options.filename,"w")

Mesh_File.write( "%\n" )
Mesh_File.write( "% Problem dimension\n" )
Mesh_File.write( "%\n" )
Mesh_File.write( "NDIME=2\n" )
Mesh_File.write( "%\n" )
Mesh_File.write( "% Inner elements\n" )
Mesh_File.write( "%\n" )
Mesh_File.write( "NELEM=%s\n" % ((nNode-1)*(mNode-1)))


iElem = 0
for jNode in range(mNode-1):
    for iNode in range(nNode-1):
        iPoint = jNode*nNode + iNode
        jPoint = jNode*nNode + iNode + 1
        kPoint = (jNode + 1)*nNode + iNode
        mPoint = (jNode + 1)*nNode + (iNode + 1)
        Mesh_File.write( "9 \t %s \t %s \t %s \t %s \t %s\n" % (iPoint, jPoint, mPoint, kPoint, iElem) )
        iElem = iElem + 1

nPoint = (nNode)*(mNode)
Mesh_File.write( "%\n" )
Mesh_File.write( "NPOIN=%s\n" % ((nNode)*(mNode)) )
iPoint = 0

for jNode in range(mNode):
    for iNode in range(nNode):
        xCoord = float(iNode)/float(nNode-1)
        yCoord = float(jNode)/float(mNode-1)

        if xAdapt == "True":
          
          x0 = 0.0; x1 = 0.12; x2 = 0.30; x3 = 1.0;
          fp3 = 10.0;  fp0 = fp3; fpp2 = 0.0; fpp1 = -fpp2;
          if (xCoord >= x0) and (xCoord < x1):
              incr = x1-x0; incr2 = x1*x1-x0*x0; incr3 = x1*x1*x1-x0*x0*x0
              d = ((1.0 - fp0) - fpp1*incr)/ (3.0*incr2-6.0*x1*incr); c =  (fpp1/2.0) -3.0*d*x1
              b = 1.0 - 2.0*c*x1 - 3.0*d*x1*x1; a = x1 - b*x1-c*x1*x1-d*x1*x1*x1
              x = a + b*xCoord + c*xCoord*xCoord + d*xCoord*xCoord*xCoord
          if (xCoord >= x1) and (xCoord <= x2):
              x = xCoord
          if (xCoord > x2) and (xCoord <= x3):
              incr = x2-x3; incr2 = x2*x2-x3*x3; incr3 = x2*x2*x2-x3*x3*x3
              d = ((1.0 - fp3) - fpp2*incr)/ (3.0*incr2-6.0*x2*incr); c = (fpp2/2.0) - 3.0*d*x2
              b = 1.0 - 2.0*c*x2 - 3.0*d*x2*x2; a = x2 - b*x2-c*x2*x2-d*x2*x2*x2
              x = a + b*xCoord + c*xCoord*xCoord + d*xCoord*xCoord*xCoord

          incr = x1-x0; incr2 = x1*x1-x0*x0; incr3 = x1*x1*x1-x0*x0*x0
          d = ((1.0 - fp0) - fpp1*incr)/ (3.0*incr2-6.0*x1*incr); c =  (fpp1/2.0) -3.0*d*x1
          b = 1.0 - 2.0*c*x1 - 3.0*d*x1*x1; a = x1 - b*x1-c*x1*x1-d*x1*x1*x1
          x_min = a + b*x0 + c*x0*x0 + d*x0*x0*x0
          incr = x2-x3; incr2 = x2*x2-x3*x3; incr3 = x2*x2*x2-x3*x3*x3
          d = ((1.0 - fp3) - fpp2*incr)/ (3.0*incr2-6.0*x2*incr); c = (fpp2/2.0) - 3.0*d*x2
          b = 1.0 - 2.0*c*x2 - 3.0*d*x2*x2; a = x2 - b*x2-c*x2*x2-d*x2*x2*x2
          x_max = a + b*x3 + c*x3*x3 + d*x3*x3*x3
            
          xCoord_new =  xLength*((x-x_min)/(x_max-x_min)) - float(options.offsetx)

        else:
          
          xCoord_new =  xLength*xCoord - float(options.offsetx)

        if yAdapt == "True":
          
          y0 = 0.0; y1 = 0.000001; y2 = 0.001; y3 = 1.0;
          fp3 = 10.0;  fp0 = fp3; fpp2 = 0.0; fpp1 = -fpp2;
          if (yCoord >= y0) and (yCoord < y1):
              incr = y1-y0; incr2 = y1*y1-y0*y0; incr3 = y1*y1*y1-y0*y0*y0
              d = ((1.0 - fp0) - fpp1*incr)/ (3.0*incr2-6.0*y1*incr); c =  (fpp1/2.0) -3.0*d*y1
              b = 1.0 - 2.0*c*y1 - 3.0*d*y1*y1; a = y1 - b*y1-c*y1*y1-d*y1*y1*y1
              y = a + b*yCoord + c*yCoord*yCoord + d*yCoord*yCoord*yCoord
          if (yCoord >= y1) and (yCoord <= y2):
              y = yCoord
          if (yCoord > y2) and (yCoord <= y3):
              incr = y2-y3; incr2 = y2*y2-y3*y3; incr3 = y2*y2*y2-y3*y3*y3
              d = ((1.0 - fp3) - fpp2*incr)/ (3.0*incr2-6.0*y2*incr); c = (fpp2/2.0) - 3.0*d*y2
              b = 1.0 - 2.0*c*y2 - 3.0*d*y2*y2; a = y2 - b*y2-c*y2*y2-d*y2*y2*y2
              y = a + b*yCoord + c*yCoord*yCoord + d*yCoord*yCoord*yCoord

          incr = y1-y0; incr2 = y1*y1-y0*y0; incr3 = y1*y1*y1-y0*y0*y0
          d = ((1.0 - fp0) - fpp1*incr)/ (3.0*incr2-6.0*y1*incr); c =  (fpp1/2.0) -3.0*d*y1
          b = 1.0 - 2.0*c*y1 - 3.0*d*y1*y1; a = y1 - b*y1-c*y1*y1-d*y1*y1*y1
          y_min = a + b*y0 + c*y0*y0 + d*y0*y0*y0

          incr = y2-y3; incr2 = y2*y2-y3*y3; incr3 = y2*y2*y2-y3*y3*y3
          d = ((1.0 - fp3) - fpp2*incr)/ (3.0*incr2-6.0*y2*incr); c = (fpp2/2.0) - 3.0*d*y2
          b = 1.0 - 2.0*c*y2 - 3.0*d*y2*y2; a = y2 - b*y2-c*y2*y2-d*y2*y2*y2
          y_max = a + b*y3 + c*y3*y3 + d*y3*y3*y3
          
          yCoord_new =  yLength*((y-y_min)/(y_max-y_min)) - float(options.offsety)

        else:
  
          yCoord_new =  yLength*yCoord - float(options.offsety)

        Mesh_File.write( "%15.14f \t %15.14f \t %s\n" % (xCoord_new, yCoord_new, iPoint) )
        iPoint = iPoint + 1

Mesh_File.write( "%\n" )
Mesh_File.write( "% Boundary elements\n" )
Mesh_File.write( "%\n" )
Mesh_File.write( "NMARK=4\n" )
Mesh_File.write( "MARKER_TAG= lower\n" )
Mesh_File.write( "MARKER_ELEMS=%s\n" % (nNode-1))
for iNode in range(nNode-1):
    Mesh_File.write( "3 \t %s \t %s\n" % (iNode, iNode + 1) )

Mesh_File.write( "MARKER_TAG= outlet\n" )
Mesh_File.write( "MARKER_ELEMS=%s\n" % (mNode-1))
for jNode in range(mNode-1):
    Mesh_File.write( "3 \t %s \t %s\n" % (jNode*nNode + (nNode - 1),  (jNode + 1)*nNode + (nNode - 1) ) )
Mesh_File.write( "MARKER_TAG= upper\n" )
Mesh_File.write( "MARKER_ELEMS=%s\n" % (nNode-1))
for iNode in range(nNode-1):
    Mesh_File.write( "3 \t %s \t %s\n" % ((nNode*mNode - 1) - iNode, (nNode*mNode - 1) - (iNode + 1)) )
Mesh_File.write( "MARKER_TAG= inlet\n" )
Mesh_File.write( "MARKER_ELEMS=%s\n" % (mNode-1))
for jNode in range(mNode-2, -1, -1):
    Mesh_File.write( "3 \t %s \t %s\n" % ((jNode + 1)*nNode, jNode*nNode ) )
    
Mesh_File.close()
