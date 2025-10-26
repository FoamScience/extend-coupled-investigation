# trace generated using paraview version 5.11.2
#import paraview
#paraview.compatibility.major = 5
#paraview.compatibility.minor = 11

import sys, os
#### import the simple module from the paraview
from paraview.simple import *
#### disable automatic camera reset on 'Show'
paraview.simple._DisableFirstRenderCameraReset()

case = os.environ['PWD']
if len(sys.argv) > 1:
    case = sys.argv[1]

# create a new 'OpenFOAMReader'
casefoam = OpenFOAMReader(registrationName='case.foam', FileName=f'{case}/case.foam')
casefoam.MeshRegions = ['internalMesh']
casefoam.CellArrays = ['U', 'Urel', 'p']
casefoam.Decomposepolyhedra = 0

# get animation scene
animationScene1 = GetAnimationScene()

# get the time-keeper
timeKeeper1 = GetTimeKeeper()

# update animation scene based on data timesteps
animationScene1.UpdateAnimationUsingDataTimeSteps()

UpdatePipeline(time=4500.0, proxy=casefoam)

animationScene1.GoToLast()

# create a new 'Slice'
slice1 = Slice(registrationName='Slice1', Input=casefoam)
slice1.SliceType = 'Plane'
slice1.HyperTreeGridSlicer = 'Plane'
slice1.SliceOffsetValues = [0.0]

# init the 'Plane' selected for 'SliceType'
slice1.SliceType.Origin = [0.0, 0.0, 0.10000000149011612]

# init the 'Plane' selected for 'HyperTreeGridSlicer'
slice1.HyperTreeGridSlicer.Origin = [0.0, 0.0, 0.10000000149011612]

# Properties modified on slice1.SliceType
slice1.SliceType.Normal = [0.0, 0.0, 1.0]

UpdatePipeline(time=5000.0, proxy=slice1)

# create a new 'Calculator'
calculator1 = Calculator(registrationName='Calculator1', Input=slice1)
calculator1.Function = ''

# rename source object
RenameSource('r', calculator1)

# Properties modified on calculator1
calculator1.ResultArrayName = 'r'
calculator1.Function = 'sqrt(coordsX^2 + coordsY^2)'

UpdatePipeline(time=5000.0, proxy=calculator1)

# create a new 'Calculator'
calculator1_1 = Calculator(registrationName='Calculator1', Input=calculator1)
calculator1_1.Function = ''

# rename source object
RenameSource('theta', calculator1_1)

# Properties modified on calculator1_1
calculator1_1.ResultArrayName = 'theta'
calculator1_1.Function = 'atan2(coordsY, coordsX)'

UpdatePipeline(time=5000.0, proxy=calculator1_1)

# create a new 'Calculator'
calculator1_2 = Calculator(registrationName='Calculator1', Input=calculator1_1)
calculator1_2.Function = ''

# rename source object
RenameSource('Ur', calculator1_2)

# Properties modified on calculator1_2
calculator1_2.ResultArrayName = 'Ur'
calculator1_2.Function = 'U_X*cos(theta) + U_Y*sin(theta)'

UpdatePipeline(time=5000.0, proxy=calculator1_2)

# create a new 'Calculator'
calculator1_3 = Calculator(registrationName='Calculator1', Input=calculator1_2)
calculator1_3.Function = ''

# rename source object
RenameSource('Utheta', calculator1_3)

# Properties modified on calculator1_3
calculator1_3.ResultArrayName = 'Utheta'
calculator1_3.Function = '-U_X*sin(theta) + U_Y*cos(theta)'

UpdatePipeline(time=5000.0, proxy=calculator1_3)

# create a new 'Plot Over Line'
plotOverLine1 = PlotOverLine(registrationName='PlotOverLine1', Input=calculator1_3)
plotOverLine1.Point1 = [-0.05000000074505806, -0.05000000074505806, 0.10000000149011612]
plotOverLine1.Point2 = [0.05000000074505806, 0.05000000074505806, 0.10000000149011612]

# Properties modified on plotOverLine1
plotOverLine1.Point1 = [0.025, 0.0, 0.10000000149011612]
plotOverLine1.Point2 = [0.05, 0.0, 0.10000000149011612]

UpdatePipeline(time=5000.0, proxy=plotOverLine1)

# save data
SaveData(f'{case}/data_along_radial.csv', proxy=plotOverLine1, PointDataArrays=['U', 'Ur', 'Urel', 'Utheta', 'arc_length', 'p', 'r', 'theta', 'vtkValidPointMask'])
