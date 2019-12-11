#!/usr/bin/env python
import argparse

parser = argparse.ArgumentParser(description='Generate a walking pattern in PyBullet')
parser.add_argument('n_steps', type=int, help='number of steps to be taken')
parser.add_argument('--step_duration', type=float, default=10, help='duration of a step in seconds')
parser.add_argument('--step_length',   type=float, default=0.05, help='distance travel bya one step')
parser.add_argument('--step_height',   type=float, default=0.1, help='max height of the foot in flight')
parser.add_argument('--render', type=bool, default=1, help='create the pybullet server in GUI mode.')
parser.add_argument('--test', type=bool, default=0, help='WIP record a value along the walk and print a graph at the end (default: saggital hip torque)')
parser.add_argument('--save_csv', type=bool, default=0, help='save the gait as a csv containing the leg joint values trajectories')
parser.add_argument('--print_joint_traj', type=bool, default=0, help='save graphs of the joint values as jpegs')
args = parser.parse_args()

import sys
import pybullet as p
import numpy as np
np.set_printoptions(precision=5)
import pybullet_data
import time
from math import radians, degrees
import os, inspect
currentdir = os.path.dirname(os.path.abspath(inspect.getfile(inspect.currentframe())))
parentdir = os.path.dirname(currentdir)
urdf_root = os.path.join(parentdir,"models")

ik_root = os.path.join(parentdir,"source/build")
sys.path.insert(0, ik_root)

import IPython





render = args.render
test = args.test


if render == False:
    physicsClient = p.connect(p.DIRECT)
else:
    physicsClient = p.connect(p.GUI)
   
p.setAdditionalSearchPath(pybullet_data.getDataPath())
p.setGravity(0,0,-9.79983)

planeId = p.loadURDF("plane.urdf")
teoStartOrientation = p.getQuaternionFromEuler([0,0,0])
teoStartPosition = [0,0,0.86]
urdf_path = os.path.join(urdf_root,"TEO_wobble.urdf")
teoId = p.loadURDF(urdf_path,teoStartPosition ,teoStartOrientation)

l0 = 0.1932
l1 = 0.305
l2 = 0.1625
l3 = 0.059742
l4 = 0.037508
l5 = 0.34692
l6 = 0.32992
l7 = 0.215
l8 = 0.090

l9 = 0.092
l10 = 0.330
l11 = 0.300
l12 = 0.123005

l13 = 0.146
l14 = 0.018
l15 = 0.026
l16 = 0.0175

n_joints = p.getNumJoints(teoId)
name2id={}
for i in range(n_joints):
    name2id[p.getJointInfo(teoId,i)[1]]=i

init_len = 80.85 
init_f = 0.

max_len = 110.5 
max_f = 50. 

mid_len = 92.23
mid_f = 10.

min_len = 51.2  
min_f = -50.

leg_length=l9+l10+l11+l12

def torque_spring(angle):
    x=angle*leg_length
    F=-1.47842*x - 0.00861885*x**3
    torque=F*leg_length
    return torque


n_joints = p.getNumJoints(teoId) 

jointIndices = [0] * n_joints
targetPos = [0.0] * n_joints
maxVelocities = [radians(1)] * n_joints
maxForces = [0.0] * n_joints

targetPosOffset = [0.0]*n_joints
targetPosModified = [0.0]*n_joints


for i in range(n_joints):
    jointIndices[i] = i
    maxForces[i] = p.getJointInfo(teoId,i)[10]
    targetPos[i] = 0
    
mode = p.POSITION_CONTROL

p.setJointMotorControlArray(teoId,
                            jointIndices, 
                            controlMode=mode,
                            forces=maxForces,
                            targetVelocities = maxVelocities,
                            targetPositions = targetPos
                           )
                           
                           
timestep = 1/240
p.setTimeStep(timestep)   



pos_l = np.asarray(p.getLinkState(teoId,name2id[b'l_sole_joint'])[0])
pos_r = np.asarray(p.getLinkState(teoId,name2id[b'r_sole_joint'])[0])


#
#       Step Parameters
#

# All units in S.I.
step_duration     =  args.step_duration
step_height       =  args.step_height
step_length       =  args.step_length
hip_step_width    =  (l13-l16)-0.02
com_height_offset = -0.005
com_forward_ofset = -0.01
n_steps = args.n_steps

freq= (240/20)
 
print_joint_traj = True

foot_width=(70+70)/1000
foot_length=(180+70)/1000



p.enableJointForceTorqueSensor(teoId, name2id[b'l_sole_joint'],1)
p.enableJointForceTorqueSensor(teoId, name2id[b'r_sole_joint'],1)


test_sent = []
test_measured = []


'''def test_function():
    measures = p.getJointStates(teoId, [name2id[b'l_sole_joint'], name2id[b'r_sole_joint']])
    left_foot_pos  = p.getLinkState(teoId,name2id[b'l_sole_joint'])[0]
    right_foot_pos = p.getLinkState(teoId,name2id[b'r_sole_joint'])[0]
    
    ftleft  = measures[0][2]
    ftright = measures[1][2]
        
    #   ZMP_i   = Moment_i / Force_z    
    left_zmp_x  = ftleft[4]/ftleft[2]
    left_zmp_y  = ftleft[5]/ftleft[2]
    right_zmp_x = ftright[4]/ftright[2]
    right_zmp_y = ftright[5]/ftright[2]
    
    return [ [left_zmp_x+left_foot_pos[0],left_zmp_y+left_foot_pos[1]],
             [right_zmp_x+right_foot_pos[0],right_zmp_y+right_foot_pos[1]],
             [left_foot_pos[0] ,left_foot_pos[1]],
             [right_foot_pos[0],right_foot_pos[1]]
             ]
             
   test_legend = ["left_zmp", "right_zmp", "left_foot_pos", "right_foot_pos"]

             
             '''

def test_function():
    l_hip_roll_torque = p.getJointState(teoId, name2id[b'l_hip_roll'])[3]
    r_hip_roll_torque = p.getJointState(teoId, name2id[b'r_hip_roll'])[3]
    return [l_hip_roll_torque, r_hip_roll_torque]
  
  
test_legend = ["l_hip_roll_torque", "r_hip_roll_torque"]


# Percentage in length of the suport polygon to plan the ZMP in
stable_percentage = 0.3 # 80%


ori_quat=p.getLinkState(teoId,name2id[b'r_sole_joint'])[1]


def trajectory_pol(xi,xti,xtti,xf,xtf,xttf,T):
    return lambda t:(6*t**5*xf)/T**5 - (15*t**4*xf)/T**4 + (10*t**3*xf)/T**3 + xi - (6*t**5*xi)/T**5 + (15*t**4*xi)/T**4 - (10*t**3*xi)/T**3 - (3*t**5*xtf)/T**4 + (7*t**4*xtf)/T**3 - (4*t**3*xtf)/T**2 + t*xti - (3*t**5*xti)/T**4 + (8*t**4*xti)/T**3 - (6*t**3*xti)/T**2 + (t**5*xttf)/(2.*T**3) - (t**4*xttf)/T**2 + (t**3*xttf)/(2.*T) + (t**2*xtti)/2. - (t**5*xtti)/(2.*T**3) + (3*t**4*xtti)/(2.*T**2) - (3*t**3*xtti)/(2.*T)
    
    
    
import screwIK as ik
    
leftLegIndices =[name2id[b'l_ankle_roll'],
                 name2id[b'l_ankle_pitch'],
                 name2id[b'l_knee_pitch'],
                 name2id[b'l_hip_pitch'],
                 name2id[b'l_hip_roll'],
                 name2id[b'l_hip_yaw']]
                 
                
leftArmIndices =[name2id[b'l_wrist_yaw'],
                 name2id[b'l_wrist_pitch'],
                 name2id[b'l_elbow_pitch'],
                 name2id[b'l_shoulder_pitch'],
                 name2id[b'l_shoulder_roll'],
                 name2id[b'l_shoulder_yaw']]
                 
rightLegIndices=[name2id[b'r_ankle_roll'],
                 name2id[b'r_ankle_pitch'],
                 name2id[b'r_knee_pitch'],
                 name2id[b'r_hip_pitch'],
                 name2id[b'r_hip_roll'],
                 name2id[b'r_hip_yaw']]
                 
rightArmIndices=[name2id[b'r_wrist_yaw'],
                 name2id[b'r_wrist_pitch'],
                 name2id[b'r_elbow_pitch'],
                 name2id[b'r_shoulder_pitch'],
                 name2id[b'r_shoulder_roll'],
                 name2id[b'r_shoulder_yaw']]
                 
                 
                 
                 
def printProgressBar (iteration, total, prefix = '', suffix = '', decimals = 1, length = 100, fill = '█', printEnd = "\r"):
    """
    Call in a loop to create terminal progress bar
    @params:
        iteration   - Required  : current iteration (Int)
        total       - Required  : total iterations (Int)
        prefix      - Optional  : prefix string (Str)
        suffix      - Optional  : suffix string (Str)
        decimals    - Optional  : positive number of decimals in percent complete (Int)
        length      - Optional  : character length of bar (Int)
        fill        - Optional  : bar fill character (Str)
        printEnd    - Optional  : end character (e.g. "\r", "\r\n") (Str)
    """
    if iteration == total:
        suffix = "Complete         "
    percent = ("{0:." + str(decimals) + "f}").format(100 * (iteration / float(total)))
    filledLength = int(length * iteration // total)
    bar = fill * filledLength + '-' * (length - filledLength)
    print('\r%s |%s| %s%% %s' % (prefix, bar, percent, suffix), end = printEnd)
    # Print New Line on Complete
    if iteration == total: 
        print()  
         
def reset_sim(com_0 = np.array([0,0,l9+l10+l11+l12-0.001]), orientation=np.array([0,0,0]), g=-9.79983):
    p.resetSimulation()
    p.setTimeStep(timestep)
    p.setGravity(0,0,g)
    planeId = p.loadURDF("plane.urdf")
    teoStartOrientation = p.getQuaternionFromEuler(-orientation)
    
    teoStartPosition = [0,0,com_0[2]]
    
    urdf_path = os.path.join(urdf_root,"TEO.urdf")
    teoId = p.loadURDF(urdf_path,teoStartPosition ,teoStartOrientation)
    #print("TEO pos =",p.getBasePositionAndOrientation(teoId)[0])
    n_joints = p.getNumJoints(teoId)
    name2id={}

    for i in range(n_joints):
        name2id[p.getJointInfo(teoId,i)[1]]=i
    n_joints = p.getNumJoints(teoId)

    jointIndices = [0] * n_joints
    targetPos = [0.0] * n_joints
    maxVelocities = [radians(1)] * n_joints
    maxForces = [0.0] * n_joints


    for i in range(n_joints):
        jointIndices[i] = i
        maxForces[i] = p.getJointInfo(teoId,i)[10]
        targetPos[i] = 0


    rightLegAngles= ik.rl_com_from_foot(com_0 , orientation)
    leftLegAngles = ik.ll_com_from_foot(com_0 , orientation)

    '''
    print("com    =",com_0)
    print("com_f  =",com_0+np.array([0,l16,-l12]))
    print()
    print("Hst0_r =",[0,-l16,l9+l10+l11+l12])
    print("rightLegP(com_0) =",rightLegP(com_0))
    print("rightLegAngles =",[degrees(q) for q in rightLegAngles])
    print()
    print("Hst0_l =",[0,l16,l9+l10+l11+l12])
    print("leftLegP(com_0) =",leftLegP(com_0))
    print("leftLegAngles =",[degrees(q) for q in leftLegAngles])
    '''
    
    for i, index in enumerate(rightLegIndices):
        if (not np.isnan(rightLegAngles[i])) and (rightLegAngles[i] > p.getJointInfo(teoId, index)[8] and rightLegAngles[i] < p.getJointInfo(teoId, index)[9]):
            targetPos[index]=rightLegAngles[i]


    for i, index in enumerate(leftLegIndices):
        if (not np.isnan(leftLegAngles[i])) and ( leftLegAngles[i] > p.getJointInfo(teoId, index)[8] and leftLegAngles[i] < p.getJointInfo(teoId, index)[9]):
            targetPos[index]=leftLegAngles[i]
            

    for jointIndex in range(n_joints):
        p.resetJointState(teoId,
                          jointIndex,
                          targetPos[jointIndex])

    mode = p.POSITION_CONTROL
    p.setJointMotorControlArray(teoId,
                                jointIndices, 
                                controlMode=mode,
                                forces=maxForces,
                                targetPositions = targetPos
                               )
    return teoId

def repos(teoId, com_0 = np.array([0,0,l9+l10+l11+l12-0.001]), orientation=np.array([0,0,0])):
    
    p.resetBasePositionAndOrientation(teoId, com_0+np.array([8.67, -0.03, 24.06])/1000,  p.getQuaternionFromEuler(orientation))
    
    n_joints = p.getNumJoints(teoId)
    name2id={}
    
    for i in range(n_joints):
        name2id[p.getJointInfo(teoId,i)[1]]=i
    n_joints = p.getNumJoints(teoId)

    jointIndices = [0] * n_joints
    targetPos = [0.0] * n_joints
    maxVelocities = [radians(1)] * n_joints
    maxForces = [0.0] * n_joints

    
    

    for i in range(n_joints):
        jointIndices[i] = i
        maxForces[i] = p.getJointInfo(teoId,i)[10]
        targetPos[i] = 0


    rightLegAngles= ik.rl_com_from_foot(com_0 , orientation)
    leftLegAngles = ik.ll_com_from_foot(com_0 , orientation)

    for i, index in enumerate(rightLegIndices):
        if (not np.isnan(rightLegAngles[i])) and (rightLegAngles[i] > p.getJointInfo(teoId, index)[8] and rightLegAngles[i] < p.getJointInfo(teoId, index)[9]):
            targetPos[index]=rightLegAngles[i]

    for i, index in enumerate(leftLegIndices):
        if (not np.isnan(leftLegAngles[i])) and ( leftLegAngles[i] > p.getJointInfo(teoId, index)[8] and leftLegAngles[i] < p.getJointInfo(teoId, index)[9]):
            targetPos[index]=leftLegAngles[i]

    for jointIndex in range(n_joints):
        p.resetJointState(teoId,
                          jointIndex,
                          targetPos[jointIndex])

    mode = p.POSITION_CONTROL
    p.setJointMotorControlArray(teoId,
                                jointIndices, 
                                controlMode=mode,
                                forces=maxForces,
                                targetPositions = targetPos
                               )

offset_torque = maxForces[name2id[b'l_hip_roll']]
max_offset =  1.2581265067842706
min_offset =  0.7139893352971366


def torqueToPosOffs(torque):
    #print("torque",torque)
    #print("MaxTorque", maxForces[name2id[b'l_hip_roll']])
    #return 0
    offs = -np.tanh(torque/offset_torque)
    if offs > 0:
        offs = offs*max_offset
    else:
        offs = offs*min_offset
    #print("offs",offs)
    return radians(offs)


frictioncoff= 0.9
p.changeDynamics(planeId, -1,
                 lateralFriction=frictioncoff*10000,
                 spinningFriction=0.1,
                 rollingFriction=0.01 )









    
#p.setRealTimeSimulation(True)

com_0 = np.array([com_forward_ofset,0,l12+l11+l10+l9+com_height_offset])
repos(teoId,com_0=com_0)

n_timesteps = int((step_duration/4)/timestep)


r_angles_plot = []
l_angles_plot = []
r_sent_angles_plot = []
l_sent_angles_plot = []

lfoot_plot = []
rfoot_plot = []
com_plot = []




com_x_trajectory = 0
com_y_trajectory = trajectory_pol(0,0,0,(hip_step_width-0),0,0,n_timesteps)
com_z_trajectory = 0

# in euler angles
orientation = [0,0,0]



############          Change weight to left foot           ################

for t in range(n_timesteps):
    p.stepSimulation()
    
    test_measured.append(test_function())
    
    
    printProgressBar(t+1, n_timesteps, prefix = 'Change weight to left foot:', length = 20)

    com_trajectory = np.array([0,
                               com_y_trajectory(t),
                               0])
    
    com_r_t = com_0 - com_trajectory
    com_l_t = com_0 - com_trajectory
    
    
    leftLegAngles =ik.ll_com_from_foot(com_r_t ,orientation)
    rightLegAngles=ik.rl_com_from_foot(com_l_t ,orientation)
    
    l_sent_angles_plot.append([targetPos[index] for index in leftLegIndices])
    r_sent_angles_plot.append([targetPos[index] for index in rightLegIndices])
    com_plot.append(com_l_t)
    
    
    for i, index in enumerate(rightLegIndices):
        if (not np.isnan(rightLegAngles[i])) and (rightLegAngles[i] > p.getJointInfo(teoId, index)[8] and rightLegAngles[i] < p.getJointInfo(teoId, index)[9]):
            targetPos[index]=rightLegAngles[i]
            

    for i, index in enumerate(leftLegIndices):
        if (not np.isnan(leftLegAngles[i])) and ( leftLegAngles[i] > p.getJointInfo(teoId, index)[8] and leftLegAngles[i] < p.getJointInfo(teoId, index)[9]):
            targetPos[index]=leftLegAngles[i]
            
            
        #IPython.embed()

    
    p.setJointMotorControlArray(teoId,
                                jointIndices, 
                                controlMode=mode,
                                forces=maxForces,
                                targetVelocities = maxVelocities,
                                targetPositions = targetPos
                               )
 
    l_angles_plot.append([a[0] for a in p.getJointStates(teoId,  leftLegIndices)])
    r_angles_plot.append([a[0] for a in p.getJointStates(teoId, rightLegIndices)])

    if render == True: 
                time.sleep(timestep)
   
   
com_r_0 = com_r_t
com_l_0 = com_l_t






############          Half-Step with right foot           ################


#                                   f(0) f'(0)  f''(0)             f(T)               f'(T)                      f''(T)       T          
rfoot_x_trajectory = trajectory_pol( 0,   0,      0,         step_length/2,  step_length/int(n_timesteps/2),       0,  int(n_timesteps/2))
rfoot_y_trajectory = trajectory_pol( 0,   0,      0,                0,                  0,                         0,  int(n_timesteps/2))
rfoot_z_trajectory = trajectory_pol( 0,   0,      0,         -step_height,              0,                         0,  int(n_timesteps/2))


for t in range(int(n_timesteps/2)):
    p.stepSimulation()
    
    printProgressBar(t+1, n_timesteps, prefix = ' Half-Step with right foot:', suffix = 'Raising Foot', length = 20)
    
    rfoot_trajectory = np.array([rfoot_x_trajectory(t),
                                 rfoot_y_trajectory(t),
                                 rfoot_z_trajectory(t)])
    
    # Since we want the foot to move, we have to move the com
    # in the opposite direction.
    com_r_t = com_r_0 + rfoot_trajectory
    
    rightLegAngles=ik.rl_com_from_foot(com_r_t ,orientation)
    
    l_sent_angles_plot.append([targetPos[index] for index in leftLegIndices])
    r_sent_angles_plot.append([targetPos[index] for index in rightLegIndices])
    
    rfoot_plot.append(rfoot_trajectory)
    com_plot.append(com_r_t)
    
    test_measured.append(test_function())
    
    for i, index in enumerate(rightLegIndices):
        if (not np.isnan(rightLegAngles[i])) and (rightLegAngles[i] > p.getJointInfo(teoId, index)[8] and rightLegAngles[i] < p.getJointInfo(teoId, index)[9]):
            targetPos[index]=rightLegAngles[i]

    
    p.setJointMotorControlArray(teoId,
                                jointIndices, 
                                controlMode=mode,
                                forces=maxForces,
                                targetVelocities = maxVelocities,
                                targetPositions = targetPos
                               )
                               
    l_angles_plot.append([a[0] for a in p.getJointStates(teoId,  leftLegIndices)])
    r_angles_plot.append([a[0] for a in p.getJointStates(teoId, rightLegIndices)])
    
    if render == True: 
                time.sleep(timestep)

rfoot_trajectory0 = rfoot_trajectory

com_r_1 = com_r_t
com_l_1 = com_l_t

#                                  f(0) f'(0)                                f''(0)        f(T)      f'(T) f''(T)       T          
rfoot_x_trajectory = trajectory_pol(0, step_length/int(n_timesteps/2),         0,      step_length/2, 0,  0,  int(n_timesteps/2))
rfoot_y_trajectory = trajectory_pol(0,    0,                                   0,            0,       0,  0,  int(n_timesteps/2))
rfoot_z_trajectory = trajectory_pol(0,    0,                                   0,      step_height,   0,  0,  int(n_timesteps/2))

for t in range(int(n_timesteps/2)):
    p.stepSimulation()
    
    rfoot_trajectory = np.array([rfoot_x_trajectory(t),
                                 rfoot_y_trajectory(t),
                                 rfoot_z_trajectory(t)])

    
    printProgressBar(n_timesteps/2+t+1, n_timesteps, prefix = ' Half-Step with right foot:', suffix = 'Lowering Foot', length = 20)
    
    com_r_t = com_r_1 + rfoot_trajectory
    
    rightLegAngles=ik.rl_com_from_foot(com_r_t ,orientation)
    
    l_sent_angles_plot.append([targetPos[index] for index in leftLegIndices])
    r_sent_angles_plot.append([targetPos[index] for index in rightLegIndices])
    rfoot_plot.append(rfoot_trajectory + rfoot_trajectory)
    com_plot.append(com_r_t)
    
    test_measured.append(test_function())
    
    for i, index in enumerate(rightLegIndices):
        if (not np.isnan(rightLegAngles[i])) and (rightLegAngles[i] > p.getJointInfo(teoId, index)[8] and rightLegAngles[i] < p.getJointInfo(teoId, index)[9]):
            targetPos[index]=rightLegAngles[i]

    l_hip_roll_torque = p.getJointState(teoId, name2id[b'l_hip_roll'])[3]
    r_hip_roll_torque = p.getJointState(teoId, name2id[b'r_hip_roll'])[3]
   
    
    targetPosOffset[name2id[b'l_hip_roll']] = torqueToPosOffs(l_hip_roll_torque)
    targetPosOffset[name2id[b'r_hip_roll']] = torqueToPosOffs(r_hip_roll_torque)
    
    targetPosModified = np.array(targetPos) + np.array(targetPosOffset)
    

    p.setJointMotorControlArray(teoId,
                                jointIndices, 
                                controlMode=mode,
                                forces=maxForces,
                                targetVelocities = maxVelocities,
                                targetPositions = targetPos
                               )


    l_angles_plot.append([a[0] for a in p.getJointStates(teoId,  leftLegIndices)])
    r_angles_plot.append([a[0] for a in p.getJointStates(teoId, rightLegIndices)])
    
    if render == True: 
                time.sleep(timestep)

printProgressBar(1, 1, prefix = ' Half-Step with right foot:', suffix = 'Complete' , length = 20)

com_r_0 = com_r_t
com_l_0 = com_l_t





for i in range(n_steps):
    ########          Change weight to right foot           ################

    n_timesteps = int((step_duration/4)/timestep)

    #                                   f(0) f'(0)  f''(0)             f(T)               f'(T)                 f''(T)       T          
    rfoot_x_trajectory = trajectory_pol( 0,   0,      0,          -step_length,             0,                    0,  int(n_timesteps))
    rfoot_y_trajectory = trajectory_pol( 0,   0,      0,           2*hip_step_width,              0,                    0,  int(n_timesteps))
    rfoot_z_trajectory = trajectory_pol( 0,   0,      0,                0,                  0,                    0,  int(n_timesteps))

    #                                   f(0) f'(0)  f''(0)             f(T)               f'(T)                 f''(T)       T          
    lfoot_x_trajectory = trajectory_pol( 0,   0,      0,          -step_length,             0,                    0,  int(n_timesteps))
    lfoot_y_trajectory = trajectory_pol( 0,   0,      0,           2*hip_step_width,              0,                    0,  int(n_timesteps))
    lfoot_z_trajectory = trajectory_pol( 0,   0,      0,                 0,                 0,                    0,  int(n_timesteps))


    for t in range(n_timesteps):
        p.stepSimulation()
        
        printProgressBar(t+1, n_timesteps, prefix ='Change weight to right foot:', suffix = 'Complete', length = 20)

        rfoot_trajectory = np.array([rfoot_x_trajectory(t),
                                     rfoot_y_trajectory(t),
                                     rfoot_z_trajectory(t)])
        
        lfoot_trajectory = np.array([lfoot_x_trajectory(t),
                                     lfoot_y_trajectory(t),
                                     lfoot_z_trajectory(t)])
        
        com_r_t = com_r_0 + rfoot_trajectory
        com_l_t = com_l_0 + lfoot_trajectory
        
        rightLegAngles = ik.rl_com_from_foot(com_r_t ,orientation)
        leftLegAngles  = ik.ll_com_from_foot(com_l_t ,orientation)

        l_sent_angles_plot.append([targetPos[index] for index in leftLegIndices])
        r_sent_angles_plot.append([targetPos[index] for index in rightLegIndices])
        com_plot.append(com_l_t)
        
        test_measured.append(test_function())

        for i, index in enumerate(rightLegIndices):
            if (not np.isnan(rightLegAngles[i])) and (rightLegAngles[i] > p.getJointInfo(teoId, index)[8] and rightLegAngles[i] < p.getJointInfo(teoId, index)[9]):
                targetPos[index]=rightLegAngles[i]

        for i, index in enumerate(leftLegIndices):
            if (not np.isnan(leftLegAngles[i])) and ( leftLegAngles[i] > p.getJointInfo(teoId, index)[8] and leftLegAngles[i] < p.getJointInfo(teoId, index)[9]):
                targetPos[index]=leftLegAngles[i] 
        

       
        

        
        targetPosModified = np.array(targetPos) + np.array(targetPosOffset)
        
        p.setJointMotorControlArray(teoId,
                                    jointIndices, 
                                    controlMode=mode,
                                    forces=maxForces,
                                    targetVelocities = maxVelocities,
                                    targetPositions = targetPos
                                   )
        l_angles_plot.append([a[0] for a in p.getJointStates(teoId,  leftLegIndices)])
        r_angles_plot.append([a[0] for a in p.getJointStates(teoId, rightLegIndices)])
    
        if render == True: 
                     time.sleep(timestep)

    com_r_0 = com_r_t
    com_l_0 = com_l_t
    
    printProgressBar(1, 1, prefix = 'Change weight to right foot:', suffix = 'Complete' , length = 20)



    ########        Step with left foot          ################
    n_timesteps = int((step_duration/4)/timestep)

    #                                   f(0) f'(0)  f''(0)             f(T)               f'(T)                 f''(T)       T          
    lfoot_x_trajectory = trajectory_pol( 0,   0,      0,           step_length,  step_length/int(n_timesteps/2),  0,  int(n_timesteps/2))
    lfoot_y_trajectory = trajectory_pol( 0,   0,      0,                0,                  0,                    0,  int(n_timesteps/2))
    lfoot_z_trajectory = trajectory_pol( 0,   0,      0,          -step_height,             0,                    0,  int(n_timesteps/2))


    for t in range(int(n_timesteps/2)):
        p.stepSimulation()
        
        printProgressBar(t+1, n_timesteps, prefix = 'Swing Left Foot:', suffix = 'Raising Foot', length = 20)
        
        lfoot_trajectory = np.array([lfoot_x_trajectory(t),
                                     lfoot_y_trajectory(t),
                                     lfoot_z_trajectory(t)])
        
        # Since we want the foot to move, we have to move the com
        # in the opposite direction.
        com_l_t = com_l_0 + lfoot_trajectory
        
        leftLegAngles=ik.ll_com_from_foot(com_l_t ,orientation)
        
        l_sent_angles_plot.append([targetPos[index] for index in leftLegIndices])
        r_sent_angles_plot.append([targetPos[index] for index in rightLegIndices])
        rfoot_plot.append(rfoot_trajectory)
        com_plot.append(com_r_t)
        
        test_measured.append(test_function())
        
        for i, index in enumerate(leftLegIndices):
            if (not np.isnan(leftLegAngles[i])) and ( leftLegAngles[i] > p.getJointInfo(teoId, index)[8] and leftLegAngles[i] < p.getJointInfo(teoId, index)[9]):
                targetPos[index]=leftLegAngles[i]

        

       
        

        
        targetPosModified = np.array(targetPos) + np.array(targetPosOffset)
        
        p.setJointMotorControlArray(teoId,
                                    jointIndices, 
                                    controlMode=mode,
                                    forces=maxForces,
                                    targetVelocities = maxVelocities,
                                    targetPositions = targetPos
                                   )

        l_angles_plot.append([a[0] for a in p.getJointStates(teoId,  leftLegIndices)])
        r_angles_plot.append([a[0] for a in p.getJointStates(teoId, rightLegIndices)])

        if render == True: 
              time.sleep(timestep)

    lfoot_trajectory0 = lfoot_trajectory

    com_r_1 = com_r_t
    com_l_1 = com_l_t

    #                                  f(0) f'(0)                       f''(0)   f(T)      f'(T)  f''(T)       T          
    lfoot_x_trajectory = trajectory_pol(0, step_length/int(n_timesteps/2),0,   step_length,  0,    0,  int(n_timesteps/2))
    lfoot_y_trajectory = trajectory_pol(0,    0,                          0,        0,       0,    0,  int(n_timesteps/2))
    lfoot_z_trajectory = trajectory_pol(0,    0,                          0,   step_height, 0,    0,  int(n_timesteps/2))
    

    for t in range(int(n_timesteps/2)):
        p.stepSimulation()
        
        lfoot_trajectory = np.array([lfoot_x_trajectory(t),
                                     lfoot_y_trajectory(t),
                                     lfoot_z_trajectory(t)])

        
        printProgressBar(n_timesteps/2+t+1, n_timesteps, prefix = 'Swing Left Foot:', suffix = 'Lowering Foot', length = 20)
        
        com_l_t = com_l_1 + lfoot_trajectory
        
        leftLegAngles=ik.ll_com_from_foot(com_l_t ,orientation)

        l_sent_angles_plot.append([targetPos[index] for index in leftLegIndices])
        r_sent_angles_plot.append([targetPos[index] for index in rightLegIndices])
        rfoot_plot.append(rfoot_trajectory + rfoot_trajectory0)
        com_plot.append(com_r_t)
        
        test_measured.append(test_function())
        
        for i, index in enumerate(leftLegIndices):
            if (not np.isnan(leftLegAngles[i])) and ( leftLegAngles[i] > p.getJointInfo(teoId, index)[8] and leftLegAngles[i] < p.getJointInfo(teoId, index)[9]):
                targetPos[index]=leftLegAngles[i]


       
        

        
        targetPosModified = np.array(targetPos) + np.array(targetPosOffset)
        
        p.setJointMotorControlArray(teoId,
                                    jointIndices, 
                                    controlMode=mode,
                                    forces=maxForces,
                                    targetVelocities = maxVelocities,
                                    targetPositions = targetPos
                                   )

        l_angles_plot.append([a[0] for a in p.getJointStates(teoId,  leftLegIndices)])
        r_angles_plot.append([a[0] for a in p.getJointStates(teoId, rightLegIndices)])
        if render == True:
            time.sleep(timestep)

    com_r_0 = com_r_t
    com_l_0 = com_l_t
    printProgressBar(1, 1, prefix = 'Swing Left Foot:', suffix = 'Complete' , length = 20)





    ########          Change weight to left foot           ################

    n_timesteps = int((step_duration/4)/timestep)

    #                                   f(0) f'(0)  f''(0)             f(T)               f'(T)                 f''(T)       T          
    rfoot_x_trajectory = trajectory_pol( 0,   0,      0,          -step_length,             0,                    0,  int(n_timesteps))
    rfoot_y_trajectory = trajectory_pol( 0,   0,      0,          -2*hip_step_width,              0,                    0,  int(n_timesteps))
    rfoot_z_trajectory = trajectory_pol( 0,   0,      0,                0,                  0,                    0,  int(n_timesteps))

    #                                   f(0) f'(0)  f''(0)             f(T)               f'(T)                 f''(T)       T          
    lfoot_x_trajectory = trajectory_pol( 0,   0,      0,          -step_length,             0,                    0,  int(n_timesteps))
    lfoot_y_trajectory = trajectory_pol( 0,   0,      0,          -2*hip_step_width,              0,                    0,  int(n_timesteps))
    lfoot_z_trajectory = trajectory_pol( 0,   0,      0,                 0,                 0,                    0,  int(n_timesteps))


    for t in range(n_timesteps):
        p.stepSimulation()
        
        printProgressBar(t+1, n_timesteps, prefix ='Change weight to left foot:', length = 20)

        rfoot_trajectory = np.array([rfoot_x_trajectory(t),
                                     rfoot_y_trajectory(t),
                                     rfoot_z_trajectory(t)])
        
        lfoot_trajectory = np.array([lfoot_x_trajectory(t),
                                     lfoot_y_trajectory(t),
                                     lfoot_z_trajectory(t)])
        
        com_r_t = com_r_0 + rfoot_trajectory
        com_l_t = com_l_0 + lfoot_trajectory
        
        rightLegAngles = ik.rl_com_from_foot(com_r_t ,orientation)
        leftLegAngles  = ik.ll_com_from_foot(com_l_t ,orientation)

        l_sent_angles_plot.append([targetPos[index] for index in leftLegIndices])
        r_sent_angles_plot.append([targetPos[index] for index in rightLegIndices])
        com_plot.append(com_l_t)
        
        test_measured.append(test_function())

        for i, index in enumerate(rightLegIndices):
            if (not np.isnan(rightLegAngles[i])) and (rightLegAngles[i] > p.getJointInfo(teoId, index)[8] and rightLegAngles[i] < p.getJointInfo(teoId, index)[9]):
                targetPos[index]=rightLegAngles[i]

        for i, index in enumerate(leftLegIndices):
            if (not np.isnan(leftLegAngles[i])) and ( leftLegAngles[i] > p.getJointInfo(teoId, index)[8] and leftLegAngles[i] < p.getJointInfo(teoId, index)[9]):
                targetPos[index]=leftLegAngles[i] 
        

       
        

        
        targetPosModified = np.array(targetPos) + np.array(targetPosOffset)
        
        p.setJointMotorControlArray(teoId,
                                    jointIndices, 
                                    controlMode=mode,
                                    forces=maxForces,
                                    targetVelocities = maxVelocities,
                                    targetPositions = targetPos
                                   )
                                   
        l_angles_plot.append([a[0] for a in p.getJointStates(teoId,  leftLegIndices)])
        r_angles_plot.append([a[0] for a in p.getJointStates(teoId, rightLegIndices)])
        if render == True: 
                time.sleep(timestep)

    com_r_0 = com_r_t
    com_l_0 = com_l_t
    printProgressBar(1, 1, prefix ='Change weight to left foot:', suffix = 'Complete', length = 20)

    ########        Step with right foot          ################
    n_timesteps = int((step_duration/4)/timestep)

    #                                   f(0) f'(0)  f''(0)             f(T)               f'(T)                 f''(T)       T          
    rfoot_x_trajectory = trajectory_pol( 0,   0,      0,           step_length,  step_length/int(n_timesteps/2),  0,  int(n_timesteps/2))
    rfoot_y_trajectory = trajectory_pol( 0,   0,      0,                0,                  0,                    0,  int(n_timesteps/2))
    rfoot_z_trajectory = trajectory_pol( 0,   0,      0,          -step_height,             0,                    0,  int(n_timesteps/2))


    for t in range(int(n_timesteps/2)):
        p.stepSimulation()
        
        printProgressBar(t+1, n_timesteps, prefix = 'Swing Right Foot:', suffix = 'Raising Foot', length = 20)
        
        rfoot_trajectory = np.array([rfoot_x_trajectory(t),
                                     rfoot_y_trajectory(t),
                                     rfoot_z_trajectory(t)])
        
        # Since we want the foot to move, we have to move the com
        # in the opposite direction.
        com_r_t = com_r_0 + rfoot_trajectory
        
        rightLegAngles=ik.rl_com_from_foot(com_r_t ,orientation)
        
        l_sent_angles_plot.append([targetPos[index] for index in leftLegIndices])
        r_sent_angles_plot.append([targetPos[index] for index in rightLegIndices])
        rfoot_plot.append(rfoot_trajectory)
        com_plot.append(com_r_t)
        
        test_measured.append(test_function())
        
        for i, index in enumerate(rightLegIndices):
            if (not np.isnan(rightLegAngles[i])) and (rightLegAngles[i] > p.getJointInfo(teoId, index)[8] and rightLegAngles[i] < p.getJointInfo(teoId, index)[9]):
                targetPos[index]=rightLegAngles[i]

        targetPosModified = np.array(targetPos) + np.array(targetPosOffset)
        
        p.setJointMotorControlArray(teoId,
                                    jointIndices, 
                                    controlMode=mode,
                                    forces=maxForces,
                                    targetVelocities = maxVelocities,
                                    targetPositions = targetPos
                                   )
                                   
        l_angles_plot.append([a[0] for a in p.getJointStates(teoId,  leftLegIndices)])
        r_angles_plot.append([a[0] for a in p.getJointStates(teoId, rightLegIndices)])
        if render == True: 
                time.sleep(timestep)

    rfoot_trajectory0 = rfoot_trajectory

    com_r_1 = com_r_t
    com_l_1 = com_l_t

    #                                  f(0) f'(0)                       f''(0)   f(T)      f'(T)  f''(T)       T          
    rfoot_x_trajectory = trajectory_pol(0, step_length/int(n_timesteps/2),0,   step_length,  0,    0,  int(n_timesteps/2))
    rfoot_y_trajectory = trajectory_pol(0,    0,                          0,        0,       0,    0,  int(n_timesteps/2))
    rfoot_z_trajectory = trajectory_pol(0,    0,                          0,   step_height, 0,    0,  int(n_timesteps/2))

    for t in range(int(n_timesteps/2)):
        p.stepSimulation()
        
        rfoot_trajectory = np.array([rfoot_x_trajectory(t),
                                     rfoot_y_trajectory(t),
                                     rfoot_z_trajectory(t)])

        
        printProgressBar(n_timesteps/2+t+1, n_timesteps, prefix = 'Swing Right Foot:', suffix = 'Lowering Foot', length = 20)
        
        com_r_t = com_r_1 + rfoot_trajectory
        
        rightLegAngles=ik.rl_com_from_foot(com_r_t ,orientation)

        l_sent_angles_plot.append([targetPos[index] for index in leftLegIndices])
        r_sent_angles_plot.append([targetPos[index] for index in rightLegIndices])
        rfoot_plot.append(rfoot_trajectory + rfoot_trajectory0)
        com_plot.append(com_r_t)
        
        test_measured.append(test_function())
        
        for i, index in enumerate(rightLegIndices):
            if (not np.isnan(rightLegAngles[i])) and (rightLegAngles[i] > p.getJointInfo(teoId, index)[8] and rightLegAngles[i] < p.getJointInfo(teoId, index)[9]):
                targetPos[index]=rightLegAngles[i]


       
        

        
        targetPosModified = np.array(targetPos) + np.array(targetPosOffset)
        
        p.setJointMotorControlArray(teoId,
                                    jointIndices, 
                                    controlMode=mode,
                                    forces=maxForces,
                                    targetVelocities = maxVelocities,
                                    targetPositions = targetPos
                                   )
                                   
        l_angles_plot.append([a[0] for a in p.getJointStates(teoId,  leftLegIndices)])
        r_angles_plot.append([a[0] for a in p.getJointStates(teoId, rightLegIndices)])
        if render == True: 
                time.sleep(timestep)

    com_r_0 = com_r_t
    com_l_0 = com_l_t
    printProgressBar(1, 1, prefix ='Swing Right Foot:', suffix = 'Complete', length = 20)

import matplotlib.pyplot as plt
import matplotlib.patches as patches
from mpl_toolkits.mplot3d import Axes3D

def signal_filter(signal_array,mode='median',filter_size=5):
    if mode == 'mean':
        filter_type = np.mean
    elif mode == 'median':
        filter_type = np.median
    signal_array_filtered = np.zeros(signal_array.shape - np.array([filter_size,0]))
    for j in range(signal_array_filtered.shape[1]): 
        for i in range(signal_array_filtered.shape[0]): 
            signal_array_filtered[i,j] = filter_type(signal_array[i:i+filter_size,j]) 
    return signal_array_filtered


l_angles_plot = np.array(l_angles_plot)
r_angles_plot = np.array(r_angles_plot)

l_sent_angles_plot = signal_filter(np.array(l_sent_angles_plot),mode='mean', filter_size=10)
r_sent_angles_plot = signal_filter(np.array(r_sent_angles_plot),mode='mean', filter_size=10)

l_sent_angles_plot = signal_filter(np.array(l_sent_angles_plot),mode='median', filter_size=10)
r_sent_angles_plot = signal_filter(np.array(r_sent_angles_plot),mode='median', filter_size=10)

l_sent_angles_plot = signal_filter(np.array(l_sent_angles_plot),mode='median', filter_size=10)
r_sent_angles_plot = signal_filter(np.array(r_sent_angles_plot),mode='median', filter_size=10)
  
if args.print_joint_traj:


    l_angles_plot = np.array(l_angles_plot)
    r_angles_plot = np.array(r_angles_plot)
    l_sent_angles_plot = np.array(l_sent_angles_plot)
    r_sent_angles_plot = np.array(r_sent_angles_plot)

    for i in range(r_angles_plot.shape[1]):
        printProgressBar(i, r_angles_plot.shape[1]*2, prefix = 'Printing Graphs:', length = 20)
        plt.plot(r_angles_plot[:,i])
        plt.plot(r_sent_angles_plot[:,i])
        plt.plot([p.getJointInfo(teoId, rightLegIndices[i])[8]]*r_angles_plot.shape[0],'r--')
        plt.plot([p.getJointInfo(teoId, rightLegIndices[i])[9]]*r_angles_plot.shape[0],'r--')
        plt.legend([['r_ankle_roll',
                    'r_ankle_pitch',
                    'r_knee_pitch',
                    'r_hip_pitch',
                    'r_hip_roll',
                    'r_hip_yaw'][i]+' Measured',
                    ['r_ankle_roll',
                    'r_ankle_pitch',
                    'r_knee_pitch',
                    'r_hip_pitch',
                    'r_hip_roll',
                    'r_hip_yaw'][i]+' Sent',
                    "joint limit"])
        plt.grid()
        fig = plt.gcf()
        fig.set_size_inches(13.5, 10.5)
        fig.savefig('r_joint_trajectory_'+str(i)+'.png', dpi=100)
        plt.clf()

    for i in range(l_angles_plot.shape[1]):
        printProgressBar(r_angles_plot.shape[1]+i, r_angles_plot.shape[1]*2, prefix = 'Printing Graphs:', length = 20)
        plt.plot(l_angles_plot[:,i])
        plt.plot(l_sent_angles_plot[:,i])
        plt.plot([p.getJointInfo(teoId, leftLegIndices[i])[8]]*l_angles_plot.shape[0],'r--')
        plt.plot([p.getJointInfo(teoId, leftLegIndices[i])[9]]*l_angles_plot.shape[0],'r--')
        plt.legend([['l_ankle_roll',
                    'l_ankle_pitch',
                    'l_knee_pitch',
                    'l_hip_pitch',
                    'l_hip_roll',
                    'l_hip_yaw'][i]+' Measured',
                    ['l_ankle_roll',
                    'l_ankle_pitch',
                    'l_knee_pitch',
                    'l_hip_pitch',
                    'l_hip_roll',
                    'l_hip_yaw'][i]+' Sent',
                    "joint limit"])
        plt.grid()
        fig = plt.gcf()
        fig.set_size_inches(13.5, 10.5)
        fig.savefig('l_joint_trajectory_'+str(i)+'.png', dpi=100)
        plt.clf()



if test:
    test_measured = np.array(test_measured)
    
    test_measured = signal_filter(test_measured, mode='median', filter_size=10)
    
    
    '''
    left_zmp  = test_measured[:,0]
    right_zmp = test_measured[:,1]
    left_foot_pos  = test_measured[:,2]
    right_foot_pos = test_measured[:,3]
    fig = plt.figure()
    ax = plt.axes(xlim=(-step_length*2, step_length*2), ylim=(-l13, l13))
    plt.plot( left_zmp[:,0]     ,      left_zmp[:,1], '-')
    plt.plot(right_zmp[:,0]     ,     right_zmp[:,1], '-')
    plt.plot( left_foot_pos[:,0], left_foot_pos[:,1], 'or')
    plt.plot(right_foot_pos[:,0],right_foot_pos[:,1], 'ok')
    '''
    for i in range(test_measured.shape[1]):
        plt.plot(test_measured[:,i])
    
    plt.legend(test_legend)
    plt.show()

p.disconnect()

l_sent_angles_plot = l_sent_angles_plot*(360/(2*np.pi))
r_sent_angles_plot = r_sent_angles_plot*(360/(2*np.pi))

if args.save_csv:
    np.savetxt("gait_traj.csv", np.concatenate([l_sent_angles_plot, r_sent_angles_plot], axis=1),  delimiter=',')


    print()
    print("Left Leg Start Pos:")
    print("set poss ({0[0]:.9f} {0[1]:.9f} {0[2]:.9f} {0[3]:.9f} {0[4]:.9f} {0[5]:.9f})".format(l_sent_angles_plot[0][::-1]))  
    print()
    print("Right Leg Start Pos:")
    print("set poss ({0[0]:.9f} {0[1]:.9f} {0[2]:.9f} {0[3]:.9f} {0[4]:.9f} {0[5]:.9f})".format(r_sent_angles_plot[0][::-1]))


