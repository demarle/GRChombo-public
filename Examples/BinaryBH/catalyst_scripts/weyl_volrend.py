# script-version: 2.0
# Catalyst state generated using paraview version 5.10.0

#### import the simple module from the paraview
from paraview.simple import *
#### disable automatic camera reset on 'Show'
paraview.simple._DisableFirstRenderCameraReset()


# ----------------------------------------------------------------
# setup views used in the visualization
# ----------------------------------------------------------------

# get the material library
materialLibrary1 = GetMaterialLibrary()

# Create a new 'Render View'
renderView1 = CreateView('RenderView')
renderView1.ViewSize = [1920, 1080]
renderView1.AxesGrid = 'GridAxes3DActor'
renderView1.OrientationAxesVisibility = 1
renderView1.StereoType = 'Crystal Eyes'
renderView1.CenterOfRotation = [192, 192, 192]
renderView1.CameraFocalPoint = [192, 192, 192]
renderView1.CameraPosition = [192, 192, -1500]
renderView1.CameraViewUp = [0.0, 1.0, 0.0]
renderView1.CameraFocalDisk = 1.0
renderView1.CameraParallelScale = 48.174170502068606
renderView1.Background = [0.02, 0.00, 0.2]
renderView1.UseColorPaletteForBackground = 0
renderView1.EnableRayTracing = 1
renderView1.BackEnd = 'OSPRay raycaster'
renderView1.AmbientSamples = 1
renderView1.OSPRayMaterialLibrary = materialLibrary1

SetActiveView(None)

# ----------------------------------------------------------------
# setup view layouts
# ----------------------------------------------------------------

# create new layout object 'Layout #1'
layout1 = CreateLayout(name='Layout #1')
layout1.AssignView(0, renderView1)
layout1.SetSize(1920, 1080)

# ----------------------------------------------------------------
# restore active view
SetActiveView(renderView1)
# ----------------------------------------------------------------

# ----------------------------------------------------------------
# setup the data processing pipelines
# ----------------------------------------------------------------

# create a new 'VisItChomboReader'
input = VisItChomboReader(registrationName='input', FileName=['/home/mr618/rds/rds-intel-vis-0hKvbqKpczQ/mr618/BinaryBH/CatalystTests/q2-insitu-off/hdf5/BinaryBHp_000000.3d.hdf5'])
input.MeshStatus = ['Mesh']
input.CellArrayStatus = ['Weyl4_Re']

# ----------------------------------------------------------------
# setup the visualization in view 'renderView1'
# ----------------------------------------------------------------

# show data from resampleToImage1
#resampleToImage1Display = Show(resampleToImage1, renderView1, 'UniformGridRepresentation')
resampleToImage1Display = Show(input, renderView1, 'AMRRepresentation')

# get color transfer function/color map for 'Weyl4_Re'
weyl4_ReLUT = GetColorTransferFunction('Weyl4_Re')
lutmodes = ["Never",
            "Grow and update on 'Apply'",
            "Grow and update every timestep",
            "Update on 'Apply'",
            "Clamp and update every timestep"]
weyl4_ReLUT.RGBPoints = [-1.00, 1.0, 0.0, 0.0,
                         -0.01, 1.0, 0.0, 0.0,
                          0.00, 0.1, 0.0, 0.1,
                          0.01, 0.0, 0.0, 1.0,
                          1.00, 0.0, 0.0, 1.0,]
weyl4_ReLUT.ColorSpace = 'RGB'
weyl4_ReLUT.ScalarRangeInitialized = 1.0

# get opacity transfer function/opacity map for 'Weyl4_Re'
weyl4_RePWF = GetOpacityTransferFunction('Weyl4_Re')
weyl4_RePWF.Points = [-1.0000, 0.4, 0.0, 0.0,
                      -0.0101, 0.4, 0.0, 0.0,
                      -0.0100, 0.0, 0.0, 0.0,
                       0.0100, 0.0, 0.0, 0.0,
                       0.0101, 0.4, 0.0, 0.0,
                       1.0000, 0.4, 0.0, 0.0]
weyl4_RePWF.ScalarRangeInitialized = 1
weyl4_ReLUT.AutomaticRescaleRangeMode = lutmodes[0]
weyl4_ReLUT.Build()

# trace defaults for the display properties.
resampleToImage1Display.SetRepresentationType('Volume')
resampleToImage1Display.ColorArrayName = ['POINTS', 'Weyl4_Re']
resampleToImage1Display.LookupTable = weyl4_ReLUT
resampleToImage1Display.ScalarOpacityUnitDistance = 1.0
resampleToImage1Display.ScalarOpacityFunction = weyl4_RePWF

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
pNG1.Trigger = 'TimeStep'

# init the 'PNG' selected for 'Writer'
pNG1.Writer.FileName = 'ChiContour_GWs_{timestep:06d}{camera}.png'
pNG1.Writer.ImageResolution = [1920, 1080]
pNG1.Writer.OverrideColorPalette = 'BlackBackground'
pNG1.Writer.Format = 'PNG'
pNG1.Writer.CameraMode = 'Phi-Theta'
pNG1.Writer.PhiResolution = 1
pNG1.Writer.ThetaResolution = 1

# ----------------------------------------------------------------
# restore active source
SetActiveSource(input)
# ----------------------------------------------------------------

# ------------------------------------------------------------------------------
# Catalyst options
from paraview import catalyst
options = catalyst.Options()
options.GlobalTrigger = 'TimeStep'
options.CatalystLiveTrigger = 'TimeStep'

# ------------------------------------------------------------------------------
if __name__ == '__main__':
    from paraview.simple import SaveExtractsUsingCatalystOptions
    # Code for non in-situ environments; if executing in post-processing
    # i.e. non-Catalyst mode, let's generate extracts using Catalyst options
    SaveExtractsUsingCatalystOptions(options)
