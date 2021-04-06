# -*- coding: utf-8 -*-
import pychrono as chrono
import pychrono.fea as fea
import pychrono.irrlicht as chronoirr
import pychrono.pardisomkl as mkl
import numpy as np
import matplotlib.pyplot as plt

chrono.SetChronoDataPath('../../../miniconda3/pkgs/pychrono-6.0.0-py37_0/Library/data/')

# Unit conversion
M_TO_L = 1e2 # convert length unit
KG_TO_W = 1e3 # convert weight unit
S_TO_T = 1e0
PI = np.pi

# Constants
floor_length = 0.03*M_TO_L
floor_thickness = 0.001*M_TO_L

beam_arc_length = 0.01*M_TO_L # Full angle of the curved beam
beam_radius = 0.003*M_TO_L
beam_angle = beam_arc_length/beam_radius 
beam_thickness = 0.0001*M_TO_L
beam_length = 0.05*M_TO_L

g = 0*M_TO_L/S_TO_T**2 # m/s^2

rho = 2700*KG_TO_W/(M_TO_L**3) # kg/m^3
E = 2.0e11*KG_TO_W/M_TO_L/S_TO_T**2 # kg*m/s^2/m^2
nu = 0.3
mb = 0.01

# chrono.ChCollisionModel.SetDefaultSuggestedEnvelope(1e-2) # Adjust if units change
# chrono.ChCollisionModel.SetDefaultSuggestedMargin(1e-3)

# col_mat = chrono.ChMaterialSurfaceSMC()
# Efloor = 2.0e9*KG_TO_W/M_TO_L/S_TO_T**2 # kg*m/s^2/m^2
# col_mat.SetYoungModulus(Efloor)
# col_mat.SetFriction(1.0)

step = 1e-3*S_TO_T
tfinal = 1*S_TO_T

deform = 0.02*M_TO_L
cycle = 1
    
# Helper functions
def cast_node(nb):
    feaNB = fea.CastToChNodeFEAbase(nb)
    nodeFead = fea.CastToChNodeFEAxyzD(feaNB)
    return nodeFead

def experiment():
    # Chrono Simulation
    mysystem = chrono.ChSystemSMC()
    mysystem.Set_G_acc(chrono.ChVectorD(0,-g,0))
    
    mground = chrono.ChBodyEasyBox(floor_thickness, floor_length, floor_length, 0)
    mground.SetPos(chrono.ChVectorD(-floor_thickness/2,0,0))
    mground.SetBodyFixed(True)
    mysystem.Add(mground)
    
    mslider = chrono.ChBodyEasyBox(floor_thickness, floor_length, floor_length, 0)
    mslider.SetPos(chrono.ChVectorD(beam_length+floor_thickness/2,0,0))
    # mslider.SetBodyFixed(True)
    mysystem.Add(mslider)
    
    mmotor = chrono.ChLinkMotorLinearPosition()
    mmotor.Initialize(
        mslider,
        mground,
        chrono.ChFrameD(chrono.ChVectorD(0,0,0),chrono.Q_from_AngZ(PI/2))
    )
    # mmotor_position = MotorPosition()
    mmotor_position = chrono.ChFunction_Sine(0,cycle/tfinal,deform)
    mmotor.SetMotorFunction(mmotor_position)
    mysystem.Add(mmotor)
    
    mmesh = fea.ChMesh()
    
    num_div_x = 25
    num_div_y = 0
    num_div_z = 5
    
    num_node_x = num_div_x+1
    num_node_y = num_div_y+1
    num_node_z = num_div_z+1
    
    num_elements = num_div_x*num_div_z*(num_div_y+1)
    num_nodes = num_node_x*num_node_z*num_node_y
    
    dx = beam_length/num_div_x
    dy = beam_thickness
    dz = beam_angle/num_div_z # rad
    
    # Nodes
    for i in range(num_nodes):
        # Position of node 
        x = i%num_node_x*dx
        t = (i//num_node_x)%num_node_z*dz+(PI-beam_angle)/2 # angle from positive x
        r = (i//(num_node_x*num_node_z)*dy+beam_radius)
        y = r*np.sin(t)
        z = r*np.cos(t)
        
        Dx = 0
        Dy = np.sin(t)
        Dz = np.cos(t)
        
        # If nodes added to element in CCW then -y
        node = fea.ChNodeFEAxyzD(
            chrono.ChVectorD(x,y,z),
            chrono.ChVectorD(Dx,Dy,Dz),
        )
        node.SetMass(0)
    
        mmesh.AddNode(node)
        
    
    # Elements
    mat = fea.ChMaterialShellANCF(rho, E, nu)
    for i in range(num_elements):
        # counter clockwise of a rectangle
        x = i%num_div_x    
        y = i//(num_div_x*num_div_z)
        z = (i//num_div_x)%num_div_z
    
        nodeA = z*num_node_x+x+y*(num_node_x*num_node_z)
        nodeB = z*num_node_x+x+1+y*(num_node_x*num_node_z)
        nodeC = (z+1)*num_node_x+x+1+y*(num_node_x*num_node_z)
        nodeD = (z+1)*num_node_x+x+y*(num_node_x*num_node_z)
    
        element = fea.ChElementShellANCF()
        element.SetNodes(
            cast_node(mmesh.GetNode(nodeA)),
            cast_node(mmesh.GetNode(nodeB)),
            cast_node(mmesh.GetNode(nodeC)),
            cast_node(mmesh.GetNode(nodeD))
        )
        element.SetDimensions(dx, beam_radius*np.sin(dz/2)*2)
        element.AddLayer(dy, 0*chrono.CH_C_DEG_TO_RAD, mat)
        
        element.SetAlphaDamp(mb)
        element.SetGravityOn(False)
        
        mmesh.AddElement(element)
    
    mysystem.Add(mmesh)
    
    for i in range(num_nodes):
        # Position of node 
        x = i%num_node_x
        
        if x == 0: 
            node_root = cast_node(mmesh.GetNode(i))
            node_root.SetFixed(True)
    
        if x == num_node_x-1: 
            node_end = cast_node(mmesh.GetNode(i))
            mjoint = fea.ChLinkPointFrameGeneric(False,True,True)
            mjoint.Initialize(node_end,mslider,node_end.GetPos())
            mysystem.Add(mjoint)
            
    
    # Visuals
    mground_color = chrono.ChColorAsset()
    mground_color.SetColor(chrono.ChColor(0,1,0))
    mground.AddAsset(mground_color)
    
    mslider_color = chrono.ChColorAsset()
    mslider_color.SetColor(chrono.ChColor(0,0,1))
    mslider.AddAsset(mslider_color)
    
    vmeshA = fea.ChVisualizationFEAmesh(mmesh)
    vmeshA.SetFEMdataType(fea.ChVisualizationFEAmesh.E_PLOT_SURFACE)
    vmeshA.SetWireframe(True)
    mmesh.AddAsset(vmeshA)
    
    # vmeshB = fea.ChVisualizationFEAmesh(mmesh)
    # vmeshB.SetFEMdataType(fea.ChVisualizationFEAmesh.E_PLOT_NONE)
    # vmeshB.SetFEMglyphType(fea.ChVisualizationFEAmesh.E_GLYPH_NODE_CSYS)
    # vmeshB.SetSymbolsThickness(0.2)
    # mmesh.AddAsset(vmeshB)
    
    # Solver and stepper
    mkl_solver = mkl.ChSolverPardisoMKL()
    mysystem.SetSolver(mkl_solver)
    
    hht_stepper = chrono.ChTimestepperHHT(mysystem)
    hht_stepper.SetStepControl(False)
    mysystem.SetTimestepper(hht_stepper)
    
    application = chronoirr.ChIrrApp(mysystem, "Curve beam", chronoirr.dimension2du(1024, 768))
    application.AddTypicalSky()
    application.AddTypicalLights()
    application.AddTypicalCamera(chronoirr.vector3df(beam_length/2, 0, -floor_length*1.5),chronoirr.vector3df(beam_length/2, 0, 0))
    application.AssetBindAll()
    application.AssetUpdateAll()
    # application.SetShowInfos(True)
    # application.SetVideoframeSaveInterval(int(1/step/25)) # N frame per unit time
    # application.SetVideoframeSave(True)
    
    def drawSysFrame(scale=0.01*M_TO_L):
        chronoirr.ChIrrTools().drawSegment(
            application.GetVideoDriver(),
            chrono.ChVectorD(0,0,0),
            chrono.ChVectorD(scale,0,0),
            chronoirr.SColor(1,255,0,0)
        )
        chronoirr.ChIrrTools().drawSegment(
            application.GetVideoDriver(),
            chrono.ChVectorD(0,0,0),
            chrono.ChVectorD(0,scale,0),
            chronoirr.SColor(1,0,255,0)
        )
        chronoirr.ChIrrTools().drawSegment(
            application.GetVideoDriver(),
            chrono.ChVectorD(0,0,0),
            chrono.ChVectorD(0,0,scale),
            chronoirr.SColor(1,0,0,255)
        )    
    
    application.SetTimestep(step)
    
    y = []
    f = []
    while application.GetDevice().run():
        y.append(mslider.GetPos().y/M_TO_L*1000)    
        f.append(mmotor.GetMotorForce()/KG_TO_W/M_TO_L*S_TO_T**2)
        
        application.BeginScene()
        application.DrawAll()
        drawSysFrame()    
        application.DoStep()
        application.EndScene()
        
        if mysystem.GetChTime() > tfinal: # in system seconds
              application.GetDevice().closeDevice()   
    
    return y, f

plt.figure()
for beam_arc_length, line in zip([0.01*M_TO_L],['-']):
    for beam_radius, color in zip([0.004*M_TO_L,0.006*M_TO_L,0.008*M_TO_L,0.010*M_TO_L],['r','g','b','k']):
        beam_angle = beam_arc_length/beam_radius 
        
        y, f = experiment()
        plt.plot(y,f,line+color,label='{:.0f},{:.0f}'.format(beam_radius/M_TO_L*1000,beam_arc_length/M_TO_L*1000))

plt.title('Force vs Deformation, Aluminum Curved Beam, ANCF')        
plt.ylabel('Force[N]')
plt.xlabel('Deformation[mm]')
plt.legend(loc='upper left')