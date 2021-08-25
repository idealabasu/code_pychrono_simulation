# -*- coding: utf-8 -*-
import pychrono as chrono
import pychrono.fea as fea
import pychrono.irrlicht as chronoirr
import pychrono.pardisomkl as mkl
import numpy as np
import matplotlib.pyplot as plt

chrono.SetChronoDataPath('./chrono_data/')

# Unit conversion
M_TO_L = 1e2 # convert length unit
KG_TO_W = 1e3 # convert weight unit
S_TO_T = 1e0
PI = np.pi

# Constants
floor_length = 0.15*M_TO_L
floor_width = 0.3*M_TO_L
floor_thickness = 0.001*M_TO_L

body_length = 0.015*M_TO_L
body_width = 0.04*M_TO_L
body_height = 0.012*M_TO_L

beam_arc_length = 0.01*M_TO_L # Full angle of the curved beam
beam_radius = 0.004*M_TO_L
beam_thickness = 0.0001*M_TO_L
beam_length = 0.01*M_TO_L
beam_start_angle = -PI/2
beam_angle = beam_arc_length/beam_radius
beam_chord_length = np.sin(beam_angle/2)*beam_radius*2
beam_chord_center_offset = np.cos(beam_angle/2)*beam_radius

slider_length = 0.01*M_TO_L
slider_width = beam_radius-beam_chord_center_offset
slider_height = beam_chord_length

pusher_length = 0.01*M_TO_L
pusher_width = beam_radius-beam_chord_center_offset
pusher_height = beam_chord_length

peg_radius = 0.001*M_TO_L
peg_height = body_height
peg_x = beam_length+body_length+slider_length+pusher_length-0.002*M_TO_L
peg_dz = 0.02*M_TO_L
peg_num = 13
g = 9.81*M_TO_L/S_TO_T**2 # m/s^2

# Polyester film
rho = 1383*KG_TO_W/(M_TO_L**3) # kg/m^3
E = 1.9e9*KG_TO_W/M_TO_L/S_TO_T**2 # kg*m/s^2/m^2
nu = 0.3
mb = 0.01

# chrono.ChCollisionModel.SetDefaultSuggestedEnvelope(1e-2) # Adjust if units change
# chrono.ChCollisionModel.SetDefaultSuggestedMargin(1e-3)

col_mat = chrono.ChMaterialSurfaceSMC()
col_mat.SetYoungModulus(1.0e7*KG_TO_W/M_TO_L/S_TO_T**2)
col_mat.SetFriction(0.4)
# col_mat.SetRestitution(0.2)
# col_mat.SetAdhesion(0)

count = 10
step = 5e-5*S_TO_T
period = 1.5*S_TO_T
tfinal = period*count

deform = 0.02*M_TO_L

# Helper functions
def cast_node(nb):
    feaNB = fea.CastToChNodeFEAbase(nb)
    nodeFead = fea.CastToChNodeFEAxyzrot(feaNB)
    return nodeFead


def experiment():
    # Chrono Simulation
    mysystem = chrono.ChSystemSMC()
    mysystem.Set_G_acc(chrono.ChVectorD(0,-g,0))

    mfloor = chrono.ChBodyEasyBox(floor_length, floor_thickness, floor_width, 0, True, True, col_mat)
    mfloor.SetPos(chrono.ChVectorD(0,-floor_thickness/2,floor_width/2-body_width))
    mfloor.SetBodyFixed(True)
    mysystem.Add(mfloor)

    rho_body = 0.3*KG_TO_W/body_height/body_length/body_width
    mbody = chrono.ChBodyEasyBox(body_length, body_height, body_width, rho_body, True, True, col_mat)
    mbody.SetPos(chrono.ChVectorD(body_length/2,body_height/2,0))
    # mbody.SetBodyFixed(True)
    mysystem.Add(mbody)

    mslider = chrono.ChBodyEasyBox(slider_length, slider_height, slider_width, 0)
    mslider.SetPos(chrono.ChVectorD(body_length+slider_length/2,body_height/2,slider_width/2))
    mysystem.Add(mslider)

    mpusher = chrono.ChBodyEasyBox(pusher_length, pusher_height, pusher_width, 0, True, True, col_mat)
    mpusher.SetPos(chrono.ChVectorD(body_length+slider_length+beam_length+pusher_length/2,body_height/2,pusher_width/2))
    mysystem.Add(mpusher)

    mpegs = []
    for i in range(-1,int(peg_num-1),1):
        mpeg = chrono.ChBodyEasyCylinder(peg_radius, peg_height, 0, True, True, col_mat)
        mpeg.SetPos(chrono.ChVectorD(peg_x,peg_height/2,peg_dz*i+pusher_width+peg_radius))
        mpeg.SetBodyFixed(True)
        mysystem.Add(mpeg)

        mpegs.append(mpeg)

    # Body can only move in yz plane
    mjointA = chrono.ChLinkMateGeneric(True,False,False,True,True,True)
    mjointA.Initialize(mbody,mfloor,chrono.ChFrameD(chrono.VNULL))
    mysystem.Add(mjointA)

    # Motor between body and slider
    mmotor = chrono.ChLinkMotorLinearPosition()
    mmotor.Initialize(
        mslider,
        mbody,
        chrono.ChFrameD(chrono.ChVectorD(0,0,0),chrono.Q_from_AngY(-PI/2))
    )
    mmotor_position = chrono.ChFunction_Sine(0,1/period,deform)
    mmotor.SetMotorFunction(mmotor_position)
    mysystem.Add(mmotor)

    mmesh = fea.ChMesh()

    num_div_x = 5
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
        # offset to side of body
        x = i%num_node_x*dx+body_length+slider_length
        # angle from positive z with start offset
        t = (i//num_node_x)%num_node_z*dz+(PI-beam_angle)/2+beam_start_angle
        r = i//(num_node_x*num_node_z)*dy+beam_radius
        # offset so that chord center is at body center
        y = r*np.sin(t)+body_height/2
        # offset so that center of chord crosses origin
        z = r*np.cos(t)-beam_chord_center_offset


        # If nodes added to element in CCW then -y
        node = fea.ChNodeFEAxyzrot(chrono.ChFrameD(
            chrono.ChVectorD(x,y,z),
            chrono.QUNIT
        ))
        # node.SetMass(0)

        mmesh.AddNode(node)


    # Elements
    melasticity = fea.ChElasticityReissnerIsothropic(E, nu)
    mdamping = fea.ChDampingReissnerRayleigh(melasticity,mb)
    mat = fea.ChMaterialShellReissner(melasticity, None, mdamping)
    mat.SetDensity(rho)
    for i in range(num_elements):
        # counter clockwise of a rectangle
        x = i%num_div_x
        y = i//(num_div_x*num_div_z)
        z = (i//num_div_x)%num_div_z

        nodeA = z*num_node_x+x+y*(num_node_x*num_node_z)
        nodeB = z*num_node_x+x+1+y*(num_node_x*num_node_z)
        nodeC = (z+1)*num_node_x+x+1+y*(num_node_x*num_node_z)
        nodeD = (z+1)*num_node_x+x+y*(num_node_x*num_node_z)

        element = fea.ChElementShellReissner4()
        element.SetNodes(
            cast_node(mmesh.GetNode(nodeA)),
            cast_node(mmesh.GetNode(nodeB)),
            cast_node(mmesh.GetNode(nodeC)),
            cast_node(mmesh.GetNode(nodeD))
        )
        element.AddLayer(dy, 0*chrono.CH_C_DEG_TO_RAD, mat)

        mmesh.AddElement(element)

    mysystem.Add(mmesh)

    for i in range(num_nodes):
        # Position of node
        x = i%num_node_x

        # Fixe nodes to slider
        if x == 0:
            node_root = cast_node(mmesh.GetNode(i))
            mjoint = chrono.ChLinkMateGeneric(True,True,True,True,True,True)
            mjoint.Initialize(mslider,node_root,node_root.Frame())
            mysystem.Add(mjoint)

        # Fix nodes to pusher
        if x == num_node_x-1:
            node_end = cast_node(mmesh.GetNode(i))
            mjoint = chrono.ChLinkMateGeneric(True,True,True,True,True,True)
            mjoint.Initialize(mpusher,node_end,node_end.Frame())
            mysystem.Add(mjoint)

    # Visuals
    mfloor_color = chrono.ChColorAsset()
    mfloor_color.SetColor(chrono.ChColor(0,1,0))
    mfloor.AddAsset(mfloor_color)

    mbody_color = chrono.ChColorAsset()
    mbody_color.SetColor(chrono.ChColor(0,0,1))
    mbody.AddAsset(mbody_color)

    mslider_color = chrono.ChColorAsset()
    mslider_color.SetColor(chrono.ChColor(1,0,0))
    mslider.AddAsset(mslider_color)

    mpusher_color = chrono.ChColorAsset()
    mpusher_color.SetColor(chrono.ChColor(0,0,1))
    mpusher.AddAsset(mpusher_color)


    for mpeg in mpegs:
        mpeg_color = chrono.ChColorAsset()
        mpeg_color.SetColor(chrono.ChColor(1,0,0))
        mpeg.AddAsset(mpeg_color)

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
    application.AddTypicalCamera(chronoirr.vector3df(0, body_height*12.5, floor_width/2-body_width),chronoirr.vector3df(0, 0, floor_width/2-body_width))
    application.AssetBindAll()
    application.AssetUpdateAll()
    # application.SetShowInfos(True)
    # application.SetVideoframeSaveInterval(int(1/step/120)) # N frames per unit time
    # application.SetVideoframeSave(True)

    def drawSysFrame(s=0.01*M_TO_L):
        chronoirr.drawSegment(application.GetVideoDriver(),chrono.ChVectorD(0,0,0),chrono.ChVectorD(s,0,0),chronoirr.SColor(1,255,0,0))
        chronoirr.drawSegment(application.GetVideoDriver(),chrono.ChVectorD(0,0,0),chrono.ChVectorD(0,s,0),chronoirr.SColor(1,0,255,0))
        chronoirr.drawSegment(application.GetVideoDriver(),chrono.ChVectorD(0,0,0),chrono.ChVectorD(0,0,s),chronoirr.SColor(1,0,0,255))

    application.SetTimestep(step)


    t = []
    z = []
    while application.GetDevice().run():
        t.append(mysystem.GetChTime()/period)
        z.append(mbody.GetPos().z/M_TO_L*1000)

        application.BeginScene()
        application.DrawAll()
        drawSysFrame()
        application.DoStep()
        application.EndScene()

        if mysystem.GetChTime() > tfinal: # in system seconds
              application.GetDevice().closeDevice()

    return t, z


plt.figure()
plt.title('Distance vs Swing Count, Polyester Film Curved Beam')
for period in np.array([0.5])*S_TO_T:
    tfinal = period*count
    t, z = experiment()
    plt.plot(t,z, label='{:.2f}Hz'.format(1/period))
plt.ylabel('Distance[mm]')
plt.xlabel('Count')
plt.legend(loc='upper left')
