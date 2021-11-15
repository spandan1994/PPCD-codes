# PPCD-codes

### Description
Codes for verifying absolute and almost sure stability of 2-dimensional Probabilistic Piecewise Continuous Derivative systems.

### System Dependencies
The code was tested on Ubuntu 20.04 LTS.

1. Install latest version of Python (3.9 or higher) and make it default:
```
sudo add-apt-repository ppa:deadsnakes/ppa
sudo apt-get update
apt list | grep python3.9
sudo apt-get install python3.9
sudo update-alternatives --install /usr/bin/python3 python3 /usr/bin/python3.8 1
sudo update-alternatives --install /usr/bin/python3 python3 /usr/bin/python3.9 2
sudo update-alternatives --config python3
```
Now type 2 and hit enter to select Python3.9 as default. Test the current version of Python using `python3 -V`.

2. Install latest Python packages:
```
sudo apt-get install python3.9-dev
```

3. Install pip:
```
sudo apt install python3-pip
```

4. Install NetworkX:
```
pip install networkx
```

5. Install pplpy:
```
sudo apt-get install m4
sudo apt-get install libgmp-dev libmpfr-dev libmpc-dev libppl-dev cython
pip install cysignals --user
pip install gmpy2 --pre --user
pip install pplpy
```

6. Install Matplotlib:
```
sudo apt-get install libjpeg-dev zlib1g-dev
pip install -U Pillow
pip install matplotlib
```


