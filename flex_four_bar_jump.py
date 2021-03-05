# -*- coding: utf-8 -*-
import pychrono as chrono
import pychrono.fea as fea
import pychrono.irrlicht as chronoirr
import matplotlib.pyplot as plt
import numpy as np
from scipy.optimize import minimize
import pychrono.pardisomkl as mkl

chrono.SetChronoDataPath('../../../miniconda3/pkgs/pychrono-6.0.0-py37_0/Library/data/')

M_TO_L = 1e2 # convert length unit
KG_TO_W = 1e3 # convert weight unit
S_TO_T = 1e2
PI = np.pi

# Solve for a initial pose of a four bar linkage
floor_length = 0.5*M_TO_L
floor_width = 0.03*M_TO_L
floor_thickness = 0.01*M_TO_L

body_length = 0.05*M_TO_L
body_width = 0.01*M_TO_L
body_thickness = 0.015*M_TO_L

crank_length = 0.02*M_TO_L
coupler_length = 0.03*M_TO_L
output_length = 0.01*M_TO_L
output_front_length = 0.02*M_TO_L
ground_length = 0.03*M_TO_L
link_thickness = 0.002*M_TO_L
link_width = 0.01*M_TO_L
rho = 1000*KG_TO_W/(M_TO_L**3) # (1000kg/m^3)

g = 9.81*M_TO_L/S_TO_T**2 # (9.81m/s^2)
k = 5e-3*KG_TO_W*M_TO_L**2/S_TO_T**2 # kg*m/s^2*m/rad
# k = 0
b = 5e-5*KG_TO_W*M_TO_L**2/S_TO_T # kg*m/s^2*m/(rad/s)

E = 2.0e9*KG_TO_W/M_TO_L/S_TO_T**2 # kg*m/s^2/m^2
Efloor = 2.0e8*KG_TO_W/M_TO_L/S_TO_T**2 # kg*m/s^2/m^2
nu = 0.3
mb = 0.01

t_joint = 0.02*KG_TO_W*M_TO_L**2/S_TO_T**2
# t_joint = 0
step_joint = 0.05*S_TO_T

# Adjust if units change
chrono.ChCollisionModel.SetDefaultSuggestedEnvelope(1e-2)
chrono.ChCollisionModel.SetDefaultSuggestedMargin(1e-3)

col_mat = chrono.ChMaterialSurfaceSMC()
col_mat.SetYoungModulus(Efloor)
col_mat.SetFriction(1.0)
# col_mat.SetPoissonRatio(nu)
# col_mat.SetRestitution(0.0)
# col_mat.SetAdhesion(0.0)

step = 2e-5*S_TO_T

x0 = np.array([0,0,crank_length/2,crank_length/2,output_length,-ground_length,0,-ground_length])

def fourbar_solver(pts):
    pts = pts.reshape((-1,2))
    vcrank = pts[1,:]-pts[0,:]
    vcoupler = pts[2,:]-pts[1,:]
    voutput = pts[2,:]-pts[3,:] # Flip direction so initial constraint value is around 0
    vground = pts[3,:]-pts[0,:]
    
    error = []
    
    # Length
    error.append(vcrank.dot(vcrank)-crank_length**2)
    error.append(vcoupler.dot(vcoupler)-coupler_length**2)
    error.append(voutput.dot(voutput)-output_length**2)
    error.append(vground.dot(vground)-ground_length**2)
    
    # Angle
    error.append(np.arctan2(vground[1],vground[0])+PI/4)
    error.append(np.arctan2(voutput[1],voutput[0]))
    # error.append(np.arctan2(vcrank[1],vcrank[0])) # Constraint crank angle can casuse the minimize to fail
    
    # Starting point
    error.append(pts[0,0])
    error.append(pts[0,1])
    
    error = np.array(error)
    
    return error.dot(error)

# fourbar_solver(x0)

result = minimize(fourbar_solver,x0)
pts = result.x.reshape((-1,2))
print('points')
print(pts)
print('error', result.fun)

# Plot the initial pose
pts = np.vstack([pts,pts[0,:]])
plt.figure()
plt.plot(pts[:,0],pts[:,1],'ko-')
plt.axis('equal')

# Chrono Simulation
mysystem = chrono.ChSystemSMC()
mysystem.Set_G_acc(chrono.ChVectorD(0,-g,0))

def link_center(p1,p2):
    return chrono.ChVectorD(
        (p1[0]+p2[0])/2, 
        (p1[1]+p2[1])/2, 
        0
    )

def link_angle(p1,p2):
    return np.arctan2(p2[1]-p1[1],p2[0]-p1[0])

class JointForce(chrono.TorqueFunctor):
    def __init__(self):
        super(JointForce, self).__init__()
    def __call__(self, time, angle, vel, link):
        return -k*angle-b*vel

class MotorForce(chrono.ChFunction):
    def __init__(self):
          chrono.ChFunction.__init__(self)
    def Get_y(self,x):
        if x < step_joint:
            return t_joint
        else:
            return 0

def CastNode(nb):
    feaNB = fea.CastToChNodeFEAbase(nb)
    nodeFead = fea.CastToChNodeFEAxyzrot(feaNB)
    return nodeFead

# Create body
mfloor = chrono.ChBodyEasyBox(floor_length, floor_thickness, floor_width, rho, True, True, col_mat)
mfloor.SetPos(chrono.ChVectorD(pts[3,0],pts[3,1]-floor_thickness/2-link_thickness/2,0))
mfloor.SetBodyFixed(True)
mysystem.Add(mfloor)

mground = chrono.ChBodyEasyBox(ground_length, link_thickness, link_width, rho)
mground.SetPos(link_center(pts[0,:],pts[3,:]))
mground.SetRot(chrono.Q_from_AngAxis(link_angle(pts[3,:],pts[0,:]), chrono.VECT_Z))
# mground.SetBodyFixed(True)
mysystem.Add(mground)


X_cf = chrono.ChFrameD(chrono.VNULL, chrono.Q_from_AngAxis(link_angle(pts[0,:],pts[1,:]), chrono.VECT_Z))
mbody = chrono.ChBodyEasyBox(body_length, body_thickness, body_width, rho)
mbody.SetPos(X_cf*chrono.ChVectorD(crank_length-body_length/2,body_thickness/2+link_thickness/2,0))
mbody.SetRot(chrono.Q_from_AngAxis(link_angle(pts[0,:],pts[1,:]), chrono.VECT_Z))
mysystem.Add(mbody)

mcrank = chrono.ChBodyEasyBox(crank_length, link_thickness, link_width, rho)
mcrank.SetPos(link_center(pts[0,:],pts[1,:]))
mcrank.SetRot(chrono.Q_from_AngAxis(link_angle(pts[0,:],pts[1,:]), chrono.VECT_Z))
mysystem.Add(mcrank)

moutput = chrono.ChBodyEasyBox(output_length+output_front_length, link_thickness, link_width, rho, True, True, col_mat)
moutput.SetPos(link_center(pts[2,:],pts[3,:])-chrono.ChVectorD(output_front_length/2,0,0))
moutput.SetRot(chrono.Q_from_AngAxis(link_angle(pts[2,:],pts[3,:]), chrono.VECT_Z))
mysystem.Add(moutput)

mjointA = chrono.ChLinkMateGeneric(True,True,True,True,True,True)
mjointA.Initialize(mcrank,mbody,chrono.ChFrameD(chrono.VNULL))
mysystem.Add(mjointA)

mmotor = chrono.ChLinkMotorRotationTorque()
mmotor.Initialize(mcrank,
                  mground,
                  chrono.ChFrameD(chrono.VNULL)) # where to create the motor in abs.space
mmotor_force = MotorForce()
mmotor.SetMotorFunction(mmotor_force)
mysystem.Add(mmotor)

mjointC = chrono.ChLinkMateGeneric(True,True,True,True,True,False)
mjointC.Initialize(moutput,mground,chrono.ChFrameD(chrono.ChVectorD(pts[3,0],pts[3,1],0)))
mysystem.Add(mjointC)

mforceC = chrono.ChLinkRotSpringCB()
mforceC.Initialize(moutput,mground,chrono.ChCoordsysD(chrono.ChVectorD(pts[3,0],pts[3,1],0)))
jointC_force = JointForce()
mforceC.RegisterTorqueFunctor(jointC_force)
mysystem.Add(mforceC)

# mjointD = chrono.ChLinkMateGeneric(False,False,True,True,True,False)
# mjointD.Initialize(mfloor,mground,chrono.ChFrameD(chrono.VNULL))
# mysystem.Add(mjointD)

mmesh = fea.ChMesh()

offset_z = link_width/2

num_div_x = 8
num_div_z = 2

num_node_x = num_div_x+1
num_node_y = 1
num_node_z = num_div_z+1

num_elements = num_div_x*num_div_z
num_nodes = num_node_x*num_node_z

dx = coupler_length/num_div_x
dy = link_thickness
dz = link_width/num_div_z

coupler_frame = chrono.ChFrameD(chrono.ChVectorD(pts[1,0],pts[1,1],0),chrono.Q_from_AngAxis(link_angle(pts[1,:],pts[2,:]), chrono.VECT_Z))

# Nodes
for i in range(num_nodes):
    # Position of node 
    x = i%num_node_x*dx
    y = 0
    z = (i//num_node_x)%num_node_z*dz-offset_z
    
    # If nodes added to element in CCW then -y
    node = fea.ChNodeFEAxyzrot(chrono.ChFrameD(
        coupler_frame*chrono.ChVectorD(x,y,z),
        chrono.Q_from_AngAxis(np.pi, chrono.VECT_Y)
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
    row = i//num_div_x
    col = i%num_div_x
    nodeA = row*num_node_x + col
    nodeB = row*num_node_x + col+1
    nodeC = row*num_node_x + col+1+num_node_x
    nodeD = row*num_node_x + col+num_node_x

    element = fea.ChElementShellReissner4()
    element.SetNodes(CastNode(mmesh.GetNode(nodeA)),
                      CastNode(mmesh.GetNode(nodeB)),
                      CastNode(mmesh.GetNode(nodeC)),
                      CastNode(mmesh.GetNode(nodeD)))
    element.AddLayer(dy, 0*chrono.CH_C_DEG_TO_RAD, mat)
    
    mmesh.AddElement(element)

mysystem.Add(mmesh)

for i in range(num_nodes):
    # Position of node 
    x = i%num_node_x
    
    if x == 0: 
        node_root = CastNode(mmesh.GetNode(i))
        
        mjoint = chrono.ChLinkMateGeneric(True,True,True,True,True,False)
        mjoint.Initialize(mcrank,node_root,node_root.Frame())
        mysystem.Add(mjoint)
        
    if x == num_node_x-1: 
        node_root = CastNode(mmesh.GetNode(i))
        mjoint = chrono.ChLinkMateGeneric(True,True,True,True,True,False)
        mjoint.Initialize(moutput,node_root,node_root.Frame())
        mysystem.Add(mjoint)        
        

# Visuals
mcrank_color = chrono.ChColorAsset()
mcrank_color.SetColor(chrono.ChColor(1,0,0))
mcrank.AddAsset(mcrank_color)

moutput_color = chrono.ChColorAsset()
moutput_color.SetColor(chrono.ChColor(0,0,1))
moutput.AddAsset(moutput_color)

vmeshA = fea.ChVisualizationFEAmesh(mmesh)
vmeshA.SetFEMdataType(fea.ChVisualizationFEAmesh.E_PLOT_SURFACE)
vmeshA.SetWireframe(True)
mmesh.AddAsset(vmeshA)

# Solver and stepper
mysystem.SetSolverForceTolerance(1e-12)

# mysystem.SetTimestepperType(chrono.ChTimestepper.Type_HHT)
# mysystem.SetSolverType(chrono.ChSolver.Type_MINRES)
mkl_solver = mkl.ChSolverPardisoMKL()
mysystem.SetSolver(mkl_solver)

hht_stepper = chrono.ChTimestepperHHT(mysystem)
hht_stepper.SetStepControl(False)
mysystem.SetTimestepper(hht_stepper)

application = chronoirr.ChIrrApp(mysystem, "Flex four bar", chronoirr.dimension2du(1024, 768))
application.AddTypicalSky()
application.AddTypicalLights()
application.AddTypicalCamera(chronoirr.vector3df(0, output_length*3, -link_width*15),chronoirr.vector3df(0, output_length*3, 0))
application.AssetBindAll()
application.AssetUpdateAll()
# application.SetShowInfos(True)
application.SetVideoframeSaveInterval(int(1/step/10)) # 10 frame per unit time
application.SetVideoframeSave(True)

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

t = []
theta = []
h = []
d = []
while application.GetDevice().run():
    # torque.append(mmotor.GetMotorTorque()/g) # Convert to g*cm for easier comparison
    
    t.append(mysystem.GetChTime())    
    theta.append(chrono.Q_to_Euler123(moutput.GetRot()).z)
    h.append(mbody.GetPos().y)
    d.append(mbody.GetPos().x)

    application.BeginScene()
    application.DrawAll()
    drawSysFrame()    
    application.DoStep()
    application.EndScene()
    
    if mysystem.GetChTime() > 0.3*S_TO_T: # in system seconds
          application.GetDevice().closeDevice()    

h = np.array(h)-h[0]
d = -np.array(d)+d[0]
plt.figure()
plt.plot(t,h,label='height')
plt.plot(t,d,label='distance')
plt.ylabel('h,d[cm]')
plt.xlabel('t[0.01s]')
plt.legend()