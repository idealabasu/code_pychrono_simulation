import pychrono as chrono
import pychrono.fea as fea
import pychrono.irrlicht as chronoirr
import pychrono.pardisomkl as mkl
import numpy as np
import matplotlib.pyplot as plt

chrono.SetChronoDataPath('../../../miniconda3/pkgs/pychrono-6.0.0-py37_223/Library/data/')

# Unit conversion
M_TO_L = 1e0
KG_TO_W = 1e0
S_TO_T = 1e0
PI = np.pi

# Constants
beam_thickness = 0.003*M_TO_L
beam_length = 0.55*M_TO_L
beam_height = 0.04*M_TO_L

g = -9.81*M_TO_L/S_TO_T**2 # m/s^2

rho = 1220*KG_TO_W/(M_TO_L**3) # kg/m^3
E = 26e6*KG_TO_W/M_TO_L/S_TO_T**2 # kg*m/s^2/m^2
nu = 0.38
mb = 0.01

freq = 1.0
amplitude = 0*PI/2/S_TO_T

step = 1e-3*S_TO_T
tfinal = 5*S_TO_T
    
# Helper functions
def cast_node(nb):
    feaNB = fea.CastToChNodeFEAbase(nb)
    nodeFead = fea.CastToChNodeFEAxyzD(feaNB)
    return nodeFead

# Chrono Simulation
mysystem = chrono.ChSystemSMC()
mysystem.Set_G_acc(chrono.ChVectorD(0,0,g))

ground = chrono.ChBodyEasyBox(0.1, 0.1, 0.1, 0, False)
ground.SetPos(chrono.ChVectorD(0,0,0))
ground.SetBodyFixed(True)
mysystem.Add(ground)

slider = chrono.ChBodyEasyBox(0.1, 0.1, 0.1, 0, False)
slider.SetPos(chrono.ChVectorD(0,0,0))
# mslider.SetBodyFixed(True)
mysystem.Add(slider)

class MotorSpeed(chrono.ChFunction) :
  def __init__(self): 
        super().__init__()

  def Get_y(self, x) :
      period = 1/freq
      t = np.mod(x, period)
      
      if t < period/4 or t > period*3/4:
          return amplitude
      else:
          return -amplitude

motor = chrono.ChLinkMotorRotationSpeed()
motor.Initialize(
    slider,
    ground,
    chrono.ChFrameD(chrono.ChVectorD(-beam_length/2,0,0),chrono.QUNIT)
)
motor_speed = chrono.ChFunction_Sine(PI/2,freq,amplitude)
# motor_speed = MotorSpeed()
motor.SetMotorFunction(motor_speed)
mysystem.Add(motor)

body = fea.ChMesh()

num_div_x = 55
num_div_z = 4
num_node_x = num_div_x+1
num_node_z = num_div_z+1

num_elements = num_div_x*num_div_z
num_nodes = num_node_x*num_node_z

dx = beam_length/num_div_x
dy = beam_thickness
dz = beam_height/num_div_z # rad

# Nodes
for k in range(num_node_z):
    for i in range(num_node_x):
        # Position of node 
        x = i*dx-beam_length/2
        y = 0
        z = k*dz
        
        dirX = 0
        dirY = -1
        dirZ = 0
        
        # If nodes added to element in CCW then -y
        node = fea.ChNodeFEAxyzD(
            chrono.ChVectorD(x,y,z),
            chrono.ChVectorD(dirX,dirY,dirZ),
        )
        node.SetMass(0)
    
        if i <= 1: 
            joint = fea.ChLinkPointFrameGeneric(True,True,True)
            joint.Initialize(node,slider)
            mysystem.Add(joint)
        else:
            joint = fea.ChLinkPointFrameGeneric(False,True,False)
            joint.Initialize(node,slider)
            mysystem.Add(joint)            
    
        body.AddNode(node)
    

# Elements
mat = fea.ChMaterialShellANCF(rho, E, nu)
for k in range(num_div_z):
    for i in range(num_div_x):
        nodeA = i+k*num_node_x
        nodeB = i+k*num_node_x+1
        nodeC = i+(k+1)*num_node_x+1
        nodeD = i+(k+1)*num_node_x

        element = fea.ChElementShellANCF()
        element.SetNodes(
            cast_node(body.GetNode(nodeA)),
            cast_node(body.GetNode(nodeB)),
            cast_node(body.GetNode(nodeC)),
            cast_node(body.GetNode(nodeD))
        )
        element.SetDimensions(dx, dz)
        element.AddLayer(dy, 0*chrono.CH_C_DEG_TO_RAD, mat)
        
        element.SetAlphaDamp(mb)
        element.SetGravityOn(False)
        
        body.AddElement(element)

mysystem.Add(body)
        
# Visuals
vbody = fea.ChVisualizationFEAmesh(body)
vbody.SetFEMdataType(fea.ChVisualizationFEAmesh.E_PLOT_SURFACE)
vbody.SetWireframe(True)
body.AddAsset(vbody)

# Solver and stepper
mkl_solver = mkl.ChSolverPardisoMKL()
mysystem.SetSolver(mkl_solver)

hht_stepper = chrono.ChTimestepperHHT(mysystem)
hht_stepper.SetStepControl(False)
mysystem.SetTimestepper(hht_stepper)

application = chronoirr.ChIrrApp(mysystem, "Curve beam", chronoirr.dimension2du(1024, 768), chronoirr.VerticalDir_Z)
application.AddTypicalSky()
application.AddTypicalLights()
application.AddTypicalCamera(chronoirr.vector3df(0, 0.5, 0.5),chronoirr.vector3df(0, 0, 0))
application.AssetBindAll()
application.AssetUpdateAll()
# application.SetShowInfos(True)
# application.SetVideoframeSaveInterval(int(1/step/25)) # N frame per unit time
# application.SetVideoframeSave(True)

application.SetTimestep(step)

t = []
f = []
while application.GetDevice().run():
    print('time: {:.4f}'.format(mysystem.GetChTime()))
    t.append(mysystem.GetChTime())
    motor_react_force = motor.Get_react_force()
    f.append([motor_react_force.x,motor_react_force.y,motor_react_force.z])
    
    application.BeginScene()
    application.DrawAll()
    
    # Draw axis for scale and orientation
    chronoirr.drawSegment(application.GetVideoDriver(),chrono.ChVectorD(0,0,0),chrono.ChVectorD(1,0,0),chronoirr.SColor(1,255,0,0))
    chronoirr.drawSegment(application.GetVideoDriver(),chrono.ChVectorD(0,0,0),chrono.ChVectorD(0,1,0),chronoirr.SColor(1,0,255,0))
    chronoirr.drawSegment(application.GetVideoDriver(),chrono.ChVectorD(0,0,0),chrono.ChVectorD(0,0,1),chronoirr.SColor(1,0,0,255))
    
    application.DoStep()
    application.EndScene()
    
    if mysystem.GetChTime() > tfinal: # in system seconds
          application.GetDevice().closeDevice()
          
t = np.array(t)
f = np.array(f)

plt.figure()
plt.plot(t,f[:,0])