#add parent dir to find package. Only needed for source code build, pip install doesn't need it.
import os, inspect
currentdir = os.path.dirname(os.path.abspath(inspect.getfile(inspect.currentframe())))
parentdir = os.path.dirname(os.path.dirname(currentdir))

from teoGymEnv import TeoGymEnv
import time

from stable_baselines.common.policies import MlpPolicy
from stable_baselines.common.vec_env import DummyVecEnv
from stable_baselines import PPO2


def main():

  environment = TeoGymEnv(renders=True, isDiscrete=False, maxSteps=10000000)
  env = DummyVecEnv([lambda: environment])  # The algorithms require a vectorized environment to run

  model = PPO2(MlpPolicy, env, verbose=1)
  model.learn(total_timesteps=20000)

if __name__ == "__main__":
  main()
