# script-version: 2.0
# Catalyst state generated using paraview version 5.9.1

#### import the simple module from the paraview
from paraview.simple import *
#### disable automatic camera reset on 'Show'
paraview.simple._DisableFirstRenderCameraReset()

print("into coprocess")

# ----------------------------------------------------------------
# setup views used in the visualization
# ----------------------------------------------------------------

# Create a new 'Render View'
renderView1 = CreateView('RenderView')
renderView1.UseColorPaletteForBackground = 0
renderView1.ViewSize = [1024, 768]
renderView1.AxesGrid = 'GridAxes3DActor'
renderView1.CenterOfRotation = [32.0, 32.0, 32.0]
renderView1.CameraPosition = [32.0, 60.0, 96.0]
renderView1.CameraFocalPoint = [32.0, 32.0, 32.0]
renderView1.CameraViewUp = [1.0, 0.0, 0.0]
renderView1.EnableRayTracing = 1
renderView1.BackEnd = 'OSPRay raycaster'
renderView1.SamplesPerPixel = 1
renderView1.Background = [0.0, 0.0, 0.0]

SetActiveView(None)

# ----------------------------------------------------------------
# setup view layouts
# ----------------------------------------------------------------

# create new layout object 'Layout #1'
layout1 = CreateLayout(name='Layout #1')
layout1.AssignView(0, renderView1)
layout1.SetSize(1024, 768)

# ----------------------------------------------------------------
# restore active view
SetActiveView(renderView1)
# ----------------------------------------------------------------

# ----------------------------------------------------------------
# setup the data processing pipelines
# ----------------------------------------------------------------

# create a new 'VisItChomboReader'
input = VisItChomboReader(registrationName='input', FileName=['/home/miren/NR/GRChombo-public/Examples/BinaryBH/hdf5/BinaryBHp_000000.3d.hdf5', '/home/miren/NR/GRChombo-public/Examples/BinaryBH/hdf5/BinaryBHp_000001.3d.hdf5', '/home/miren/NR/GRChombo-public/Examples/BinaryBH/hdf5/BinaryBHp_000002.3d.hdf5', '/home/miren/NR/GRChombo-public/Examples/BinaryBH/hdf5/BinaryBHp_000003.3d.hdf5'])
input.MeshStatus = ['Mesh']
input.CellArrayStatus = ['chi']

# ----------------------------------------------------------------
# setup the visualization in view 'renderView1'
# ----------------------------------------------------------------

# show data from resampleToImage1
display = Show(input, renderView1, 'AMRRepresentation')

# get opacity transfer function/opacity map for 'chi'
chiPWF = GetOpacityTransferFunction('chi')

chiPWF.Points = [0.0, 0.8, 0.5, 0.0,
                 0.05, 0.8, 0.5, 0.0,
                 0.15, 0.4, 0.5, 0.0,
                 0.4, 0.1, 0.5, 0.0,
                 0.8, 0.001, 0.5, 0.0,
                 0.9673850536346436, 0.0, 0.5, 0.0]
"""
chiPWF.Points = [0.0, 0.0, 0.5, 0.0,
                 0.5, 0.3, 0.5, 0.0,
                 1.0, 0.6, 0.5, 0.0]
"""

# trace defaults for the display properties.
display.SetRepresentationType('Volume')
display.ColorArrayName = ['CELLS', 'chi']
display.ScalarOpacityUnitDistance = 1.0

display.RescaleTransferFunctionToDataRange(True, True)

# ----------------------------------------------------------------
# setup color maps and opacity mapes used in the visualization
# note: the Get..() functions create a new object, if needed
# ----------------------------------------------------------------

# ----------------------------------------------------------------
# setup extractors
# ----------------------------------------------------------------

# create extractor
pNG1 = CreateExtractor('PNG', renderView1, registrationName='PNG1')
# trace defaults for the extractor.
# init the 'PNG' selected for 'Writer'
pNG1.Writer.FileName = 'OSPRayVolumeChi_{timestep:06d}{camera}.png'
pNG1.Writer.ImageResolution = [1024, 768]
pNG1.Writer.Format = 'PNG'

import math
import numpy

cam = GetActiveCamera()
cam.SetRoll(90)
pos = numpy.array(cam.GetPosition())
at = numpy.array(cam.GetFocalPoint())
up = numpy.array(cam.GetViewUp())

d = numpy.linalg.norm(pos-at)
v1 = (pos-at)/d
v2 = numpy.cross(v1,up)

numsteps = 7
for x in range(0,numsteps):
    print("T0 CAMERA VIEW ", x)
    theta = x*2*math.pi/numsteps
    pnow = at+d*((v1*math.cos(theta))+(v2*math.sin(theta)))
    cam.SetPosition(pnow)
    # SaveScreenshot("foo_"+str(x)+".png") #not respecting viewsize for some reason
    WriteImage("foo_"+str(x)+".png")


# ----------------------------------------------------------------
# restore active source
SetActiveSource(input)
# ----------------------------------------------------------------

# ------------------------------------------------------------------------------
# Catalyst options
from paraview import catalyst
options = catalyst.Options()
options.GlobalTrigger = 'TimeStep'
options.EnableCatalystLive = 0
#options.GlobalTrigger.UseStartTimeStep = 1
#options.GlobalTrigger.StartTimeStep = 1
#options.GlobalTrigger.Frequency = 10

# ------------------------------------------------------------------------------
if __name__ == '__main__':
    from paraview.simple import SaveExtractsUsingCatalystOptions
    # Code for non in-situ environments; if executing in post-processing
    # i.e. non-Catalyst mode, let's generate extracts using Catalyst options
    SaveExtractsUsingCatalystOptions(options)
