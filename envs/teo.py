import os, inspect
currentdir = os.path.dirname(os.path.abspath(inspect.getfile(inspect.currentframe())))
parentdir = os.path.dirname(os.path.dirname(currentdir))

import pybullet as p
import numpy as np
import copy
import math
import pybullet_data


class TEO:

  def __init__(self, urdfRootPath=pybullet_data.getDataPath(), timeStep=0.01):
    self.urdfRootPath = urdfRootPath
    self.timeStep = timeStep
    self.useInverseKinematics = 1
    self.useSimulation = 1
    self.useOrientation = 1

    #lower limits for null space
    self.ll = [-.967, -2, -2.96, 0.19, -2.96, -2.09, -3.05]
    #upper limits for null space
    self.ul = [.967, 2, 2.96, 2.29, 2.96, 2.09, 3.05]
    #joint ranges for null space
    self.jr = [5.8, 4, 5.8, 4, 5.8, 4, 6]
    #restposes for null space
    self.rp = [0, 0, 0, 0.5 * math.pi, 0, -math.pi * 0.5 * 0.66, 0]
    #joint damping coefficents
    self.jd = [
        0.00001, 0.00001, 0.00001, 0.00001, 0.00001, 0.00001, 0.00001, 0.00001, 0.00001, 0.00001,
        0.00001, 0.00001, 0.00001, 0.00001
    ]
    self.reset()

  def reset(self):

    # measured at root waist
    teoStartOrientation = p.getQuaternionFromEuler([0,0,0])
    teoStartPosition = [0,0,0.820932]
    self.teoUid = p.loadURDF(os.path.join(self.urdfRootPath, "TEO.urdf"), teoStartPosition ,teoStartOrientation)

    self.numJoints = p.getNumJoints(self.teoUid)
    

    # set up the robot parameters
    self.jointIndices = list(range(self.numJoints))
    self.jointPositions = [0.0] * self.numJoints
    self.maxVelocities = [math.radians(1)] * self.numJoints
    self.maxForces = [0.0] * self.numJoints
    self.jointPositions = [0] * self.numJoints  

    for i in range(self.numJoints):
        self.maxForces[i] = p.getJointInfo(self.teoUid,i)[10]

    # Set the starting position/action
    for jointIndex in self.jointIndices:
      p.resetJointState(self.teoUid, jointIndex, self.jointPositions[jointIndex])
      p.setJointMotorControl2(self.teoUid,
                              jointIndex,
                              p.POSITION_CONTROL,
                              targetPosition=self.jointPositions[jointIndex],
                              force=self.maxForces[jointIndex])

    self.motorNames = []
    self.motorIndices = []


    # set motor names and indices
    for i in self.jointIndices:
      jointInfo = p.getJointInfo(self.teoUid, i)
      qIndex = jointInfo[3]
      if qIndex > -1:
        #print("motorname")
        #print(jointInfo[1])
        self.motorNames.append(str(jointInfo[1]))
        self.motorIndices.append(i)
    


  def getActionDimension(self):
    return len(self.motorIndices)


  def getObservationDimension(self):
    return len(self.getObservation())


  def getObservation(self):
    observation = []

    linear, angular = p.getBaseVelocity(self.teoUid)

    p.getBaseVelocity

    observation.extend(list(linear))
    observation.extend(list(angular))

    return observation


  def applyAction(self, motorCommands):
    for action in range(len(motorCommands)):
        motor = self.motorIndices[action]
        p.setJointMotorControl2(self.teoUid,
                                motor,
                                p.POSITION_CONTROL,
                                targetPosition=motorCommands[action],
                                force=self.maxForces[action])






















