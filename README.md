# pybullet-tests
Version controll and sync for personal vesion of the teo-pybullet-version

## Requirements

Jupyter Notebooks
```bash
pip install jupyterlab
```

PyBullet
```bash
pip install pybullet
```

CMake
```bash
sudo apt install cmake
```

Eigen
```bash
cd
mkdir -p repos && cd repos
git clone https://gitlab.com/libeigen/eigen.git
cd eigen
cd build_dir
cmake ..
sudo make install
```

Pytest
```bash
pip install pytest
```

pybind11
```bash
cd
mkdir -p repos && cd repos
git clone https://github.com/pybind/pybind11.git
cd pybind
cd build
cmake ..
sudo make install
```


## Install

Clone the repo
```bash
mkdir -p repos && cd repos
git clone https://github.com/imontesino/pybullet-tests.git
```


Compile the inverse kinematics python module
```bash
cd source
mkdir build
cd build
cmake ..
sudo make install
´´´

## Usage

Examples of jupyter notebooks are inside the `/notebooks` folder.

To run a walking simulation use the script `/scripts/walk.py`

```bash
cd scripts
chmod +x walk.py
./walk.py
```

The command `./walk.py -h` lists the parametes and options.
