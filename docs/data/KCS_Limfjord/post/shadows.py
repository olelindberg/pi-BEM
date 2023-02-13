#!/usr/bin/env python
import vtkmodules.all as vtk

from matplotlib import cm



# noinspection PyUnresolvedReferences
import vtkmodules.vtkInteractionStyle
# noinspection PyUnresolvedReferences
import vtkmodules.vtkRenderingOpenGL2
from vtkmodules.vtkCommonColor import vtkNamedColors
from vtkmodules.vtkFiltersSources import (
    vtkCubeSource,
    vtkSphereSource
)
from vtkmodules.vtkIOGeometry import (
    vtkBYUReader,
    vtkOBJReader,
    vtkSTLReader
)
from vtkmodules.vtkIOPLY import vtkPLYReader
from vtkmodules.vtkRenderingCore import (
    vtkActor,
    vtkLight,
    vtkPolyDataMapper,
    vtkRenderWindow,
    vtkRenderWindowInteractor,
    vtkRenderer
)
from vtkmodules.vtkRenderingOpenGL2 import (
    vtkCameraPass,
    vtkRenderPassCollection,
    vtkSequencePass,
    vtkShadowMapPass
)


def get_program_parameters():
    import argparse
    description = 'Read a polydata file of a surface and display it with shadows.'
    epilogue = '''
If no file is entered a sphere is used.
   '''
    parser = argparse.ArgumentParser(description=description, epilog=epilogue,
                                     formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument('filename', default=None, nargs='?', help='Enter a polydata file e.g cow.g.')
    args = parser.parse_args()
    return args.filename


def ReadPolyData(file_name):
    import os
    path, extension = os.path.splitext(file_name)
    extension = extension.lower()
    if extension == '.ply':
        reader = vtkPLYReader()
        reader.SetFileName(file_name)
        reader.Update()
        poly_data = reader.GetOutput()
    elif extension == '.vtp':
        reader = vtkXMLpoly_dataReader()
        reader.SetFileName(file_name)
        reader.Update()
        poly_data = reader.GetOutput()
    elif extension == '.obj':
        reader = vtkOBJReader()
        reader.SetFileName(file_name)
        reader.Update()
        poly_data = reader.GetOutput()
    elif extension == '.stl':
        reader = vtkSTLReader()
        reader.SetFileName(file_name)
        reader.Update()
        poly_data = reader.GetOutput()
    elif extension == '.vtk':
        reader = vtkpoly_dataReader()
        reader.SetFileName(file_name)
        reader.Update()
        poly_data = reader.GetOutput()
    elif extension == '.g':
        reader = vtkBYUReader()
        reader.SetGeometryFileName(file_name)
        reader.Update()
        poly_data = reader.GetOutput()
    else:
        # Return a None if the extension is unknown.
        poly_data = None
    return poly_data

def makeLookupTable(xmin,xmax):

    cmap    = cm.get_cmap('jet')
    ctf     = vtk.vtkColorTransferFunction()
    for i in range(256):
        x   = i/255.0
        rgb = cmap(x)[:3]
        x   = xmin+(xmax-xmin)*x
        ctf.AddRGBPoint(x,rgb[0],rgb[1],rgb[2])

    return ctf

def main():

    colors  = vtk.vtkNamedColors()


    fn = get_program_parameters()
    if fn:
        polyData = ReadPolyData(fn)
    else:
        # Use a sphere
        source = vtkSphereSource()
        source.SetThetaResolution(100)
        source.SetPhiResolution(100)
        source.Update()
        polyData = source.GetOutput()


    reader = vtk.vtkXMLUnstructuredGridReader()
    reader.SetFileName('/home/ole/Projects/pi-BEM/docs/data/KCS_Limfjord/output/result_scalar_results_0.vtu')
    reader.Update()  # Needed because of GetScalarRange
    output = reader.GetOutput()

    # Create the mapper that corresponds the objects of the vtk.vtk file
    # into graphics elements
    mapper = vtk.vtkDataSetMapper()
    mapper.SetLookupTable(makeLookupTable(-1000,1600))
    mapper.SetInputData(output)
    mapper.SetScalarRange(-1000,1600)
    # mapper.ScalarVisibilityOff()
    mapper.SelectColorArray("hydrodynamic_pressure")
    mapper.SetScalarModeToUsePointFieldData()
    mapper.InterpolateScalarsBeforeMappingOn()


    # Create the Actor
    actor = vtk.vtkActor()
    actor.SetMapper(mapper)
    actor.GetProperty().SetInterpolationToGouraud()
    #actor.GetProperty().EdgeVisibilityOn()
    #actor.GetProperty().SetLineWidth(2.0)
#    actor.GetProperty().SetColor(colors.GetColor3d("MistyRose"))
    actor.GetProperty().SetSpecular(1.0)
    actor.GetProperty().SetDiffuse(0.4)
    actor.GetProperty().SetAmbient(0.4)
    actor.GetProperty().SetSpecularPower(30.0)
    actor.GetProperty().SetOpacity(1.0)

    backface = vtk.vtkProperty()
    backface.SetColor(colors.GetColor3d('Tomato'))
    actor.SetBackfaceProperty(backface)



    colors = vtkNamedColors()
    colors.SetColor('HighNoonSun', [255, 255, 251, 255])  # Color temp. 5400°K
    colors.SetColor('100W Tungsten', [255, 214, 170, 255])  # Color temp. 2850°K

    renderer = vtkRenderer()
    renderer.AddActor(actor)

    renderer.SetBackground(colors.GetColor3d('Silver'))

    renderWindow = vtkRenderWindow()
    renderWindow.SetSize(640, 480)
    renderWindow.AddRenderer(renderer)

    interactor = vtkRenderWindowInteractor()
    interactor.SetRenderWindow(renderWindow)

    light1 = vtkLight()
    light1.SetFocalPoint(0, 0, 0)
    light1.SetPosition(0, 1, 0.2)
    light1.SetColor(colors.GetColor3d('HighNoonSun'))
    light1.SetIntensity(0.5)
    renderer.AddLight(light1)

    light2 = vtkLight()
    light2.SetFocalPoint(0, 0, 0)
    light2.SetPosition(1.0, 1.0, 1.0)
    light2.SetColor(colors.GetColor3d('100W Tungsten'))
    light2.SetIntensity(0.5)
    renderer.AddLight(light2)

    mapper = vtkPolyDataMapper()
    mapper.SetInputData(polyData)

    actor = vtkActor()
    actor.SetMapper(mapper)
    actor.GetProperty().SetAmbientColor(colors.GetColor3d('SaddleBrown'))
    actor.GetProperty().SetDiffuseColor(colors.GetColor3d('Sienna'))
    actor.GetProperty().SetSpecularColor(colors.GetColor3d('White'))
    actor.GetProperty().SetSpecular(0.51)
    actor.GetProperty().SetDiffuse(0.7)
    actor.GetProperty().SetAmbient(0.7)
    actor.GetProperty().SetSpecularPower(10.0)
    actor.GetProperty().SetOpacity(1.0)
    renderer.AddActor(actor)

    # Add a plane
    bounds = polyData.GetBounds()

    rnge = [0] * 3
    rnge[0] = bounds[1] - bounds[0]
    rnge[1] = bounds[3] - bounds[2]
    rnge[2] = bounds[5] - bounds[4]
    print('range: ', ', '.join(['{0:0.6f}'.format(i) for i in rnge]))
    expand = 1.0
    thickness = rnge[2] * 0.1
    plane = vtkCubeSource()
    plane.SetCenter((bounds[1] + bounds[0]) / 2.0,
                    bounds[2] - thickness / 2.0,
                    (bounds[5] + bounds[4]) / 2.0)
    plane.SetXLength(bounds[1] - bounds[0] + (rnge[0] * expand))
    plane.SetYLength(thickness)
    plane.SetZLength(bounds[5] - bounds[4] + (rnge[2] * expand))

    planeMapper = vtkPolyDataMapper()
    planeMapper.SetInputConnection(plane.GetOutputPort())

    planeActor = vtkActor()
    planeActor.SetMapper(planeMapper)
    renderer.AddActor(planeActor)

    renderWindow.SetMultiSamples(0)

    shadows = vtkShadowMapPass()

    seq = vtkSequencePass()

    passes = vtkRenderPassCollection()
    passes.AddItem(shadows.GetShadowMapBakerPass())
    passes.AddItem(shadows)
    seq.SetPasses(passes)

    cameraP = vtkCameraPass()
    cameraP.SetDelegatePass(seq)

    # Tell the renderer to use our render pass pipeline
    glrenderer = renderer
    glrenderer.SetPass(cameraP)

    #renderer.GetActiveCamera().SetPosition(-0.2, 0.2, 1)
    #renderer.GetActiveCamera().SetFocalPoint(0, 0, 0)
    #renderer.GetActiveCamera().SetViewUp(0, 0, 1)
    #renderer.GetActiveCamera().Azimuth(60)
    #renderer.GetActiveCamera().Elevation(-60)
    camera = renderer.GetActiveCamera()
    camera.Roll(120)
    #camera.Azimuth(40.0)
    camera.Elevation(-70.0)

    renderer.ResetCamera()
    renderer.GetActiveCamera().Dolly(8)
    renderer.ResetCameraClippingRange()
    renderWindow.SetWindowName('Shadows')
    renderWindow.Render()
    renderWindow.SetWindowName('Shadows')

    interactor.Start()


if __name__ == '__main__':
    main()