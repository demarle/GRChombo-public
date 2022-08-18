# script-version: 2.0
# Catalyst state generated using paraview version 5.9.1

#### import the simple module from the paraview
from paraview.simple import *
#### disable automatic camera reset on 'Show'
paraview.simple._DisableFirstRenderCameraReset()

print("into coprocess")


def rescaleTF(pts, newRange):
  vals = pts[::4] # get values, ignore R,G,B
  npts = len(vals)
  mV = min(vals)
  MV = max(vals)
  rV = MV-mV
  mnR = newRange[0]
  MnR = newRange[1]
  rnR = newRange[1]-newRange[0]
  for x in range(0, npts):
      pts[x*4+0] = (pts[x*4+0]-mV)/rV*rnR+mnR
  return pts

# ----------------------------------------------------------------
# setup views used in the visualization
# ----------------------------------------------------------------

# Create a new 'Render View'
renderView1 = CreateView('RenderView')
renderView1.UseColorPaletteForBackground = 0
renderView1.ViewSize = [1024, 768]
renderView1.AxesGrid = 'GridAxes3DActor'
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

pythonCalculator1 = PythonCalculator(registrationName='PythonCalculator1', Input=input)
pythonCalculator1.ArrayAssociation = 'Cell Data'
pythonCalculator1.Expression = "chi"

# ----------------------------------------------------------------
# setup the visualization in view 'renderView1'
# ----------------------------------------------------------------

# show data from resampleToImage1
display = Show(input, renderView1, 'AMRRepresentation')

# trace defaults for the display properties.
display.SetRepresentationType('Volume')
display.ColorArrayName = ['CELLS', 'chi']
display.ScalarOpacityUnitDistance = 1.0

chiCTF = GetColorTransferFunction('chi')
chiCTF.ScalarRangeInitialized = 1
chiCTF.AutomaticRescaleRangeMode = "Never"
pts = chiCTF.RGBPoints
pts = rescaleTF(pts, [0,1])
chiCTF.RGBPoints = pts;

chiPWF = GetOpacityTransferFunction('chi')
pts = [0.0, 1.0, 0.5, 0.0,
       0.03, 0.0, 0.5, 0.0,
       0.1, 0.5, 0.5, 0.0,
       0.15, 0.0, 0.5, 0.0,
       0.2, 0.1, 0.5, 0.0,
       0.3, 0.0, 0.5, 0.0,
       0.4, 0.01, 0.5, 0.0,
       0.3, 0.0, 0.5, 0.0,
       0.5, 0.01, 0.5, 0.0,
       0.6, 0.0, 0.5, 0.0,
       0.7, 0.01, 0.5, 0.0,
       0.8, 0.0, 0.5, 0.0]
pts = rescaleTF(pts, [0,0.8])
chiPWF.Points = pts


#display.RescaleTransferFunctionToDataRange(False, True) #

# ----------------------------------------------------------------
# setup color maps and opacity mapes used in the visualization
# note: the Get..() functions create a new object, if needed
# ----------------------------------------------------------------

renderView1.ResetCamera()

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
pNG1.Writer.CameraMode = 'Phi-Theta'
pNG1.Writer.PhiResolution = 1
pNG1.Writer.ThetaResolution = 4

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
