
from __future__ import print_function

from OCC.Display.SimpleGui import init_display
from OCC.Extend.DataExchange import read_iges_file
from OCC.Core.TopAbs import TopAbs_FACE
from OCC.Core.TopExp import TopExp_Explorer
from OCC.Core.TopoDS import topods_Face
from OCC.Core.BRep import BRep_Tool

import numpy as np
import vtkmodules.all as vtk
from matplotlib import cm

import pyoccvtk as occvtk

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

def sphereActor(center,radius):

    sphere = vtkSphereSource()
    sphere.SetCenter(center)
    sphere.SetRadius(radius)

    mapper = vtkPolyDataMapper()
    mapper.SetInputConnection(sphere.GetOutputPort())

    actor = vtkActor()
    actor.SetMapper(mapper)

    return actor

def arrowActor(origin,direction):

    arrow = vtk.vtkArrowSource()
    arrow.SetTipResolution(64)
    arrow.SetShaftResolution(64)


    # The X axis is a vector from start to end
    normalizedX = direction
    length = vtk.vtkMath.Norm(normalizedX)
    vtk.vtkMath.Normalize(normalizedX)

    # The Z axis is an arbitrary vector cross X
    normalizedZ = [0] * 3
    rng = vtk.vtkMinimalStandardRandomSequence()
    rng.SetSeed(8775070)  # For testing.    
    arbitrary = [0] * 3
    for i in range(0, 3):
        rng.Next()
        arbitrary[i] = rng.GetRangeValue(-10, 10)
    vtk.vtkMath.Cross(normalizedX, arbitrary, normalizedZ)
    vtk.vtkMath.Normalize(normalizedZ)

    # The Y axis is Z cross X
    normalizedY = [0] * 3
    vtk.vtkMath.Cross(normalizedZ, normalizedX, normalizedY)
    matrix = vtk.vtkMatrix4x4()

    # Create the direction cosine matrix
    matrix.Identity()
    for i in range(0, 3):
        matrix.SetElement(i, 0, normalizedX[i])
        matrix.SetElement(i, 1, normalizedY[i])
        matrix.SetElement(i, 2, normalizedZ[i])

    # Apply the transforms
    transform = vtk.vtkTransform()
    transform.Translate(origin)
    transform.Concatenate(matrix)
    transform.Scale(length, length, length)

    mapper = vtkPolyDataMapper()
    mapper.SetInputConnection(arrow.GetOutputPort())

    actor = vtkActor()
    actor.SetUserMatrix(transform.GetMatrix())
    actor.SetMapper(mapper)
    return actor

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

    shape = read_iges_file("/home/ole/dev/projects/pi-BEM/docs/data/KCS_Limfjord/limfjord_seabed.iges")
    seabed_actors = []
    topExp = TopExp_Explorer()
    topExp.Init(shape, TopAbs_FACE)
    while topExp.More():
        face = topods_Face(topExp.Current())
        surface = BRep_Tool().Surface(face)
        mapper = occvtk.OccVtk_makeMapper(surface)
        actor = vtk.vtkActor()
        actor.SetMapper(mapper)
        seabed_actors.append(actor)
        topExp.Next()

    # display, start_display, add_menu, add_function_to_menu = init_display()
    # display.DisplayShape(shapes, update=True)
    # start_display()




    for step in range(1):

        colors  = vtk.vtkNamedColors()

        reader = vtk.vtkXMLUnstructuredGridReader()
        reader.SetFileName('./../output/result_scalar_results_' + str(step) + '.vtu')
        reader.Update()  # Needed because of GetScalarRange
        output = reader.GetOutput()

        static_force = np.genfromtxt("./../output/hydrostatic_force.csv",delimiter=",")
        dynamic_force = np.genfromtxt("./../output/hydrodynamic_force.csv",delimiter=",")

        static_cp  = static_force[step,:3]
        dynamic_cp = dynamic_force[step,:3]
        hydrodynamic_force = dynamic_force[step,3:6]

        Lpp = 230.0
        Tm  = 10.8
        area = Lpp*Tm

        LCB_sphere = sphereActor([Lpp/2-0.0148*Lpp,0,Tm],4)
        LCB_sphere.GetProperty().SetColor(colors.GetColor3d("Black"))

        static_cpsphere = sphereActor(static_cp,4)
        static_cpsphere.GetProperty().SetColor(colors.GetColor3d("Red"))
        
        dynamic_cpsphere = sphereActor(dynamic_cp,4)
        dynamic_cpsphere.GetProperty().SetColor(colors.GetColor3d("Green"))

        b       = 16
        rho     = 1000
        u       = 1.799786422884671
        Awp     = 6227.87
        Vol     = 52030.0
        sinkage = hydrodynamic_force[2]/(rho*Awp)
        print(sinkage)
        pressure_scale = 1/2*rho*u**2
        force_scale = pressure_scale*area

        xforcearrow = arrowActor([0,0,0],[hydrodynamic_force[0]/force_scale,0,0])
        yforcearrow = arrowActor(dynamic_cp,[0,hydrodynamic_force[1]/force_scale*20000,0])
        zforcearrow = arrowActor([0,0,0],[0,0,hydrodynamic_force[2]])

        # Create the mapper that corresponds the objects of the vtk.vtk file
        # into graphics elements
        mapper = vtk.vtkDataSetMapper()
        mapper.SetLookupTable(makeLookupTable(-pressure_scale,pressure_scale))
        mapper.SetInputData(output)
        mapper.SetScalarRange(-pressure_scale,pressure_scale)
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
        #actor.GetProperty().SetSpecular(1.0)
        #actor.GetProperty().SetDiffuse(0.4)
        #actor.GetProperty().SetAmbient(0.4)
        #actor.GetProperty().SetSpecularPower(30.0)
        actor.GetProperty().SetOpacity(0.8)

        backface = vtk.vtkProperty()
        backface.SetColor(colors.GetColor3d('Tomato'))
        actor.SetBackfaceProperty(backface)

        colors.SetColor('HighNoonSun', [255, 255, 251, 255])  # Color temp. 5400°K
        colors.SetColor('100W Tungsten', [255, 214, 170, 255])  # Color temp. 2850°K

        renderer = vtkRenderer()
        # renderer.AddActor(actor)
        # renderer.AddActor(xforcearrow)
        # renderer.AddActor(yforcearrow)
        # renderer.AddActor(static_cpsphere)
        # renderer.AddActor(dynamic_cpsphere)
        # renderer.AddActor(LCB_sphere)
        
        for act in seabed_actors:
            renderer.AddActor(act)


        #renderer.AddActor(zforcearrow)

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
        #renderer.AddLight(light1)

        light2 = vtkLight()
        light2.SetFocalPoint(0, 0, 0)
        light2.SetPosition(1.0, 1.0, 1.0)
        light2.SetColor(colors.GetColor3d('100W Tungsten'))
        light2.SetIntensity(0.5)
        #renderer.AddLight(light2)

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
    #    glrenderer.SetPass(cameraP)

        #renderer.GetActiveCamera().SetPosition(-0.2, 0.2, 1)
        #renderer.GetActiveCamera().SetFocalPoint(0, 0, 0)
        #renderer.GetActiveCamera().SetViewUp(0, 0, 1)
        #renderer.GetActiveCamera().Azimuth(60)
        #renderer.GetActiveCamera().Elevation(-60)
        #camera = renderer.GetActiveCamera()
        #camera.Roll(120)
        #camera.Azimuth(40.0)
        #camera.Elevation(-70.0)

        renderer.ResetCamera()
        renderer.GetActiveCamera().Dolly(8)
        renderer.ResetCameraClippingRange()
        renderWindow.SetWindowName('Shadows')
        renderWindow.Render()
        renderWindow.SetWindowName('Shadows')


        # # screenshot code:
        # w2if = vtk.vtkWindowToImageFilter()
        # w2if.SetInput(renderWindow)
        # w2if.SetInputBufferTypeToRGB()
        # w2if.ReadFrontBufferOff()
        # w2if.Update()

        # writer = vtk.vtkPNGWriter()
        # writer.SetFileName('pictures/KCS_Limfjord_' + str(step) + '.png')
        # writer.SetInputConnection(w2if.GetOutputPort())
        # writer.Write()


    interactor.Start()


if __name__ == '__main__':
    main()