# minis
'minis' is a software for electrophysiological data analysis.
It allows you to:\
(1) detect and analyse spontaneous postsynaptic potetnials/currents in whole-cell patch clamp recordings (Dervinis and Major, 2023);\
(2) estimate the quantal size in the central synapse (Dervinis and Major, in preparation).

'minis' can be used as a standalone application with a graphical user interface (GUI), as a Matlab (Mathworks) function, or a Python (Python Software Foundation) package. 'minis' has been fully written in Matlab and its code base is open source.

**Key points**
1. You can run 'minis' in parallel fashion. To do so, go to the Optimisation options menu and click on Change options. In the Options window change the Number of parallel cores parameter to another value (e.g., 4). Click OK.
2. 'minis' can also be executed on a computing cluster. The ```source_code``` directory contains a Matlab script file called ```minisCluster.m``` that can be used to launch 'minis' on a cluster.
3. If the Down-going tick box is marked, the recording trace will be inverted.
4. Software documentation is available [here](https://github.com/dervinism/minis/blob/main/minis_documentation.pdf).
5. Read all of the instructions in the README file before you use minis.

## Installation and launch instructions
### Standalone desktop app
Run the installer located inside the minisStandalone folder (Windows: minisInstaller_web.exe, Linux: minisInstaller_web.install, macOS: minisInstaller_web.app). Simply follow the installer instructions. Uninstall the minis the same way you uninstall any other regular app.

On Windows and Linux you can launch the minis app in the same way you launch other apps. On macOS you launch the minis app by navigating to the application folder and issuing the command below:
```
./run_minis.sh <matlab_runtime_directory path>
```
Typically this command should work:
```
./run_minis.sh /Applications/MATLAB/MATLAB_Runtime/v912
```
If you launch the macOS app in the regular way, the log file will not be generated.

### Matlab functions
If you want to use minis within the Matlab environment, download the source_code folder and add it to your Matlab path. You can launch the app by typing ```minis``` in the Matlab console. To run minis without the GUI, use the ```minisHeadless``` function. You can also run 'minis' on a cluster using the script ```minisCluster.m```

### Python package on Windows
To install follow these steps:
1. Install Python 3.9 (in a separate environment if needed).
2. Open the minisPy_Windows/installer folder and run the minisPyInstaller_web.exe. Follow the installation instructions and install Matlab Runtime as part of them.
3. Open the minisPy_Windows/python_files folder and execute the following line in your terminal ```python setup.py install```.
4. Install Axon Binary File format python utility by executing the following line in your terminal ```pip install pyabf```. You are all set.

You can use [p131c_0011_sw6-10.abf](https://github.com/dervinism/minis/blob/main/p131c_0011_sw6-10.abf) file with the testPython.py script to test your installation. Make sure you adapt the script to load your files. It is important that you always import the minisPy package before you import the matlab package. [Here](https://uk.mathworks.com/help/compiler_sdk/python/initialize-the-matlab-runtime.html) you can find further info on how to initialise Matlab Runtime and minisPY.

### Python package on Linux
To install follow these steps:
1. Install Python 3.9 (in a separate environment if needed).
2. Open the minisPy_Linux/installer folder and run the minisPyInstaller_web.install. Follow the installation instructions and install Matlab Runtime as part of them.
3. Open the minisPy_Linux/python_files folder and execute the following line in your terminal ```python setup.py install```.
4. Install Axon Binary File format python utility by executing the following line in your terminal ```pip install pyabf```.
5. Update Matlab Runtime path by adding the lines below to your ~./bashrc file:
    ```
    export LD_LIBRARY_PATH="${LD_LIBRARY_PATH:+${LD_LIBRARY_PATH}:}\
    /usr/local/MATLAB/MATLAB_Runtime/v912/runtime/glnxa64:\
    /usr/local/MATLAB/MATLAB_Runtime/v912/bin/glnxa64:\
    /usr/local/MATLAB/MATLAB_Runtime/v912/sys/os/glnxa64:\
    /usr/local/MATLAB/MATLAB_Runtime/v912/extern/bin/glnxa64"
    export LD_PRELOAD="${LD_PRELOAD:+${LD_PRELOAD}:}\
    /usr/local/MATLAB/MATLAB_Runtime/v912/bin/glnxa64/glibc-2.17_shim.so"
    ```
    Make sure you also execute these commands in your terminal. You are all set.

You can use [p131c_0011_sw6-10.abf](https://github.com/dervinism/minis/blob/main/p131c_0011_sw6-10.abf) file with the testPython.py script to test your installation. Make sure you adapt the script to load your files. It is important that you always import the minisPy package before you import the matlab package. [Here](https://uk.mathworks.com/help/compiler_sdk/python/initialize-the-matlab-runtime.html) you can find further info on how to initialise Matlab Runtime and minisPY.

### Python package on macOS
To install follow these steps:
1. Install Python 3.9 (in a separate environment if needed).
2. Open the minisPy_macOS/installer folder and run the minisPyInstaller_web.app. Follow the installation instructions and install Matlab Runtime as part of them.
3. Open the minisPy_macOS/python_files folder and execute the following line in your terminal ```python3 setup.py install```.
4. Install Axon Binary File format python utility by executing the following line in your terminal ```pip3 install pyabf```.
5. Update Matlab Runtime path by adding the lines below to your ~./bash_profile and ~./zshenv files:
    ```
    export DYLD_LIBRARY_PATH="${DYLD_LIBRARY_PATH:+${DYLD_LIBRARY_PATH}:}\
    /Applications/MATLAB/MATLAB_Runtime/v912/runtime/maci64:\
    /Applications/MATLAB/MATLAB_Runtime/v912/bin/maci64:\
    /Applications/MATLAB/MATLAB_Runtime/v912/sys/os/maci64:\
    /Applications/MATLAB/MATLAB_Runtime/v912/extern/bin/maci64"
    ```
    Make sure you also execute this command in your terminal. \
6. Change python executable associated with mwpython by executing the lines below in your terminal (replace angle brackets with an actual path):
    ```
    echo 'export PYTHONHOME=<python3 installation directory>' >> ~/.zshenv
    source ~/.zshrc
    echo 'export PYTHONHOME=<python3 installation directory>' >> ~/.bash_profile
    source ~/.bash_profile
    ```
    PYTHONHOME needs to be set to an actual directory that has a bin subdirectory that contains the python3.9 executable. Typically the correct directory is /Library/Frameworks/Python.framework/Versions/3.9.
7. You are all set. Just make sure to use mwpython to run any code that involves loading Python packages compiled in Matlab. For example,
    ```
    /Applications/MATLAB/MATLAB_Runtime/v912/bin/mwpython testPython.py
    ```
You can use [p131c_0011_sw6-10.abf](https://github.com/dervinism/minis/blob/main/p131c_0011_sw6-10.abf) file with the testPython.py script to test your installation. Make sure you adapt the script to load your files. It is important that you always import the minisPy package before you import the matlab package. [Here](https://uk.mathworks.com/help/compiler_sdk/python/initialize-the-matlab-runtime.html) you can find further info on how to initialise Matlab Runtime and minisPY.

## Documentation
Software user documentation file [minis_documentation.pdf](https://github.com/dervinism/minis/blob/main/minis_documentation.pdf) is available for a detailed explanation of how to use the software graphical user interface. Examples on how to use programming interfaces in Matlab and Python are given in [testMatlab.m](https://github.com/dervinism/minis/blob/main/testMatlab.m), [testMatlab_preload.m](https://github.com/dervinism/minis/blob/main/testMatlab_preload.m), [testPython.py](https://github.com/dervinism/minis/blob/main/testPython.py), and [testPython_preload.py](https://github.com/dervinism/minis/blob/main/testPython_preload.py) files.

## References
Dervinis, M, Major, G (2023) bioRxiv 2022.03.20.485046; doi: https://doi.org/10.1101/2022.03.20.485046
