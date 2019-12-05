import os, inspect
currentdir = os.path.dirname(os.path.abspath(inspect.getfile(inspect.currentframe())))
parentdir = os.path.dirname(currentdir)

import math
import pickle
import gym
from gym import spaces
from gym.utils import seeding
import numpy as np
import time
import pybullet as p
import teo
import random
import pybullet_data
from pkg_resources import parse_version

largeValObservation = 100

RENDER_HEIGHT = 720
RENDER_WIDTH = 960


class TeoGymEnv(gym.Env):
  metadata = {'render.modes': ['human', 'rgb_array'], 'video.frames_per_second': 30}
  def __init__(self,
               urdfRoot=os.path.join(parentdir,"models"),
               actionRepeat=1,
               isEnableSelfCollision=True,
               renders=True,
               isDiscrete=False,
               maxSteps=1000):

    # Simulation parameters
    self._isDiscrete = isDiscrete
    self._timeStep = 1. / 240.
    self._urdfRoot = urdfRoot
    self._actionRepeat = actionRepeat
    self._isEnableSelfCollision = isEnableSelfCollision
    self._observation = []
    self._envStepCounter = 0
    self._renders = renders
    self._maxSteps = maxSteps
    self.terminated = 0

    # Camera parameters to render view
    self._cam_dist = 1.3
    self._cam_yaw = 180
    self._cam_pitch = -40


    # Select Physics servers (Normally p.DIRECT, meaning running locally with no gui)
    self._p = p
    if self._renders:
      cid = p.connect(p.GUI)
      p.resetDebugVisualizerCamera(1.3, 180, -41, [0.52, -0.2, -0.33])
    else:
      p.connect(p.DIRECT)

    p.setAdditionalSearchPath(pybullet_data.getDataPath())


    # Setup the models in the simulation
    self.seed()
    self.reset()
    observationDim = len(self.getExtendedObservation())

    observation_high = np.array([largeValObservation] * observationDim)

    # Whether or not the action space can only take specific values
    if (self._isDiscrete):
      self.action_space = spaces.Discrete(30)
    else:
      action_dim = 28

      with open ('action_lows', 'rb') as fp:
        action_lows = pickle.load(fp)
      with open ('action_highs', 'rb') as fp:
        action_highs = pickle.load(fp)

      self.action_space = spaces.Box(low=np.array(action_lows), high=np.array(action_highs), dtype=np.float32)

    self.observation_space = spaces.Box(-observation_high, observation_high)
    self.viewer = None

    self.observation_space = spaces.Box(-observation_high, observation_high)
    self.viewer = None



  def reset(self):
    self.terminated = 0
    p.resetSimulation()
    p.setPhysicsEngineParameter(numSolverIterations=150)
    p.setTimeStep(self._timeStep)
    p.setGravity(0, 0, -10)

    p.loadURDF("plane.urdf")

    # defined in teo.py
    self._teo = teo.TEO(urdfRootPath=self._urdfRoot, timeStep=self._timeStep)

    self._envStepCounter = 0
    p.stepSimulation()
    self._observation = self.getExtendedObservation()
    return np.array(self._observation)

  
  def __del__(self):
    p.disconnect()

  # randomness for better generalization
  def seed(self, seed=None):
    self.np_random, seed = seeding.np_random(seed)
    return [seed]

  # no need for extended observation at the moment
  def getExtendedObservation(self):
    self._observation = self._teo.getObservation()
    return self._observation


  def step(self, action):
    if (self._isDiscrete):
      # to be implemented, discrete makes no sense in position control
      pass
    else:
      #print("action[0]=", str(action[0]))
      realAction = action


    return self.step2(realAction)


  def step2(self, action):
    # actions may spann miltiple timesteps
    for i in range(self._actionRepeat):
      self._teo.applyAction(action)
      p.stepSimulation()
      if self._termination():
        break
      self._envStepCounter += 1
    
    # if we are rendering we need to pause the image onscreen
    if self._renders:
      time.sleep(self._timeStep)
    self._observation = self.getExtendedObservation()

    #print("self._envStepCounter")
    #print(self._envStepCounter)

    done = self._termination()
    npaction = np.array(
        [action[i] for i in [4,5,6,7,8,9,10,11,12,13,14,15]]
    )  # penalize arm movement to prevent flailing
 
    actionCost = np.linalg.norm(npaction) * 10. + (0.820932-self._p.getBasePositionAndOrientation(self._teo.teoUid)[0][2])*100
    #print("actionCost")
    #print(actionCost)
    reward = self._reward() - actionCost
    #print("reward")
    #print(reward)

    #print("len=%r" % len(self._observation))

    return np.array(self._observation), reward, done, {}


  def render(self, mode="rgb_array", close=False):
    if mode != "rgb_array":
      return np.array([])

    base_pos, orn = self._p.getBasePositionAndOrientation(self._teo.teoUid)
    view_matrix = self._p.computeViewMatrixFromYawPitchRoll(cameraTargetPosition=base_pos,
                                                            distance=self._cam_dist,
                                                            yaw=self._cam_yaw,
                                                            pitch=self._cam_pitch,
                                                            roll=0,
                                                            upAxisIndex=2)

    proj_matrix = self._p.computeProjectionMatrixFOV(fov=60,
                                                     aspect=float(RENDER_WIDTH) / RENDER_HEIGHT,
                                                     nearVal=0.1,
                                                     farVal=100.0)

    (_, _, px, _, _) = self._p.getCameraImage(width=RENDER_WIDTH,
                                              height=RENDER_HEIGHT,
                                              viewMatrix=view_matrix,
                                              projectionMatrix=proj_matrix,
                                              renderer=self._p.ER_BULLET_HARDWARE_OPENGL)
    #renderer=self._p.ER_TINY_RENDERER)

    rgb_array = np.array(px, dtype=np.uint8)
    rgb_array = np.reshape(rgb_array, (RENDER_HEIGHT, RENDER_WIDTH, 4))

    rgb_array = rgb_array[:, :, :3]
    return rgb_array

  def _termination(self):
    hipPos, _ = p.getBasePositionAndOrientation(self._teo.teoUid)
    return hipPos[2] < 0.6 or self._envStepCounter >= self._maxSteps


  def _reward(self):

    #rewards is distance in x traveled

    reward = -1000
    reward = reward + p.getBasePositionAndOrientation(self._teo.teoUid)[0][0]*1000

    return reward

  if parse_version(gym.__version__) < parse_version('0.9.6'):
    _render = render
    _reset = reset
    _seed = seed
    _step = step

