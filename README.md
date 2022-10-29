# minis

'minis' is a software for electrophysiological data analysis.
It allows you to:\
(1) detect and analyse spontaneous postsynaptic potetnials/currents in whole-cell patch clamp recordings (Dervinis and Major, 2022);\
(2) estimate the quantal size in the central synapse (Dervinis and Major, in preparation).

Download the trial version of the software here: https://github.com/dervinism/minis/releases/tag/v1.1.0  \
If downloading the executable is blocked by your web browser, make sure that the executable file is kept and not discarded. To install the trial version on your computer, please follow the installer (Windows: minisTrialInstaller_web.exe, Linux: minisTrialInstaller_web.install, macOS: minisTrialInstaller_web.app) instructions. As part of the installation process, you will have to install Matlab R2022a Runtime. The trial version comes with inbuilt data, so you can try out running the software in various modes: Detect, Detect and Compare, Automatic Distribution Fitting, and Simulate. You will not be able to load your own data.

It is a proprietary software that is currently provided free of charge. If you would like to get a software copy for your use, please email Martynas Dervinis (martynas.dervinis@gmail.com). You will need to provide the serial number of your hard-drive (Windows and Linux) or your computing device serial number (macOS). In the email please also indicate which operating system you are using. 'minis' has been fully written in Matlab (Mathworks) and is distributed as an application programming interface in the form of a packaged Matlab app or a Python (Python Software Foundation) package or as a compiled standalone desktop application with a graphical user interface.

**Instructions on how to get the required serial number**\
On Windows open your command prompt and type in ```wmic diskdrive get model,serialnumber```. \
Any internal disk will work.

On macOS open your terminal and type in ```system_profiler SPHardwareDataType | grep Serial```.

On Linux open your terminal and type in ```ls -la /dev/disk/by-uuid | grep sda1```. If your hard drive name is other than sda1, use the serial number of that hard drive instead.

**Installation instructions: Standalone desktop app**\
When you get hold of the fully functional software, you can install the standalone version of the software by running the installer inside the minis folder. Simply follow the installer instructions. Uninstall the minis the same way you unistall any other regular app.

**Installation instructions: Matlab packaged app**\
In order to install minis as a Matlab packaged app, double click minisMatlab.mlappinstall inside the minisMatlab folder and follow instructions inside Matlab. To uninstall, navigate to Matlab Apps section, right-click minisMatlab under MY APPS subsection, and uninstall it.

**Installation instructions: Python package on Windows**\
To install follow these steps:
1. Install Python 3.9 (in a separate environment if needed).
2. Open the minisPy/installer folder and run the minisPyInstaller_web.exe. Follow the installation instructions and install Matlab Runtime as part of them.
3. Open the minisPy/python_files folder and execute the following line in your terminal ```python setup.py install```.
4. Install Axon Binary File format python utility by executing the following line in your terminal ```pip install pyabf```. You are all set.

You can use [p131c_0011_sw6-10.abf](https://github.com/dervinism/minis/blob/main/p131c_0011_sw6-10.abf) file with the testPython.py script to test your installation. Make sure you adapt the script to load your files. It is important that you always import the minisPy package before you import the matlab package. [Here](https://uk.mathworks.com/help/compiler_sdk/python/initialize-the-matlab-runtime.html) you can find further info on how to initialise Matlab Runtime and minisPY.

**Installation instructions: Python package on Linux**\
To install follow these steps:
1. Install Python 3.9 (in a separate environment if needed).
2. Open the minisPy/installer folder and run the minisPyInstaller_web.install. Follow the installation instructions and install Matlab Runtime as part of them.
3. Open the minisPy/python_files folder and execute the following line in your terminal ```python setup.py install```.
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
&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp; Make sure you also execute these commands in your terminal. You are all set. \

You can use [p131c_0011_sw6-10.abf](https://github.com/dervinism/minis/blob/main/p131c_0011_sw6-10.abf) file with the testPython.py script to test your installation. Make sure you adapt the script to load your files. It is important that you always import the minisPy package before you import the matlab package. [Here](https://uk.mathworks.com/help/compiler_sdk/python/initialize-the-matlab-runtime.html) you can find further info on how to initialise Matlab Runtime and minisPY.

**Installation instructions: Python package on macOS**\
To install follow these steps:
1. Install Python 3.9 (in a separate environment if needed).
2. Open the minisPy/installer folder and run the minisPyInstaller_web.app. Follow the installation instructions and install Matlab Runtime as part of them.
3. Open the minisPy/python_files folder and execute the following line in your terminal ```python3 setup.py install```.
4. Install Axon Binary File format python utility by executing the following line in your terminal ```pip3 install pyabf```.
5. Update Matlab Runtime path by adding the lines below to your ~./bash_profile and ~./zshenv files:
```
export DYLD_LIBRARY_PATH="${DYLD_LIBRARY_PATH:+${DYLD_LIBRARY_PATH}:}\
/Applications/MATLAB/MATLAB_Runtime/v912/runtime/maci64:\
/Applications/MATLAB/MATLAB_Runtime/v912/bin/maci64:\
/Applications/MATLAB/MATLAB_Runtime/v912/sys/os/maci64:\
/Applications/MATLAB/MATLAB_Runtime/v912/extern/bin/maci64"
```
&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp; Make sure you also execute this command in your terminal. \
6. Change python executable associated with mwpython by executing the lines below in your terminal (replace angle brackets with an actual path):
```
echo 'export PYTHONHOME=<python3 installation directory>' >> ~/.zshenv
source ~/.zshrc
echo 'export PYTHONHOME=<python3 installation directory>' >> ~/.bash_profile
source ~/.bash_profile
```
&emsp;&emsp;PYTHONHOME needs to be set to an actual directory that has a bin subdirectory that contains the python3.9 \
&emsp;&emsp;executable. Typically the correct directory is /Library/Frameworks/Python.framework/Versions/3.9. \
&nbsp;&nbsp;&nbsp; 7. You are all set. Just make sure to use mwpython to run any code that involves loading Python packages compiled in Matlab. For example,
```
/Applications/MATLAB/MATLAB_Runtime/v912/bin/mwpython testPython.py
```
You can use [p131c_0011_sw6-10.abf](https://github.com/dervinism/minis/blob/main/p131c_0011_sw6-10.abf) file with the testPython.py script to test your installation. Make sure you adapt the script to load your files. It is important that you always import the minisPy package before you import the matlab package. [Here](https://uk.mathworks.com/help/compiler_sdk/python/initialize-the-matlab-runtime.html) you can find further info on how to initialise Matlab Runtime and minisPY.

**Documentation**\
Software user documentation file [minis_documentation.pdf](https://github.com/dervinism/minis/blob/main/minis_documentation.pdf) is available for a detailed explanation of how to use the software graphical user interface. Examples on how to use programming interfaces in Matlab and Python are given in [testMatlab.m](https://github.com/dervinism/minis/blob/main/testMatlab.m), [testMatlab_preload.m](https://github.com/dervinism/minis/blob/main/testMatlab_preload.m), [testPython.py](https://github.com/dervinism/minis/blob/main/testPython.py), and [testPython_preload.py](https://github.com/dervinism/minis/blob/main/testPython_preload.py) files.

**Input files**\
**The software does NOT work with version 2 ABF files**. Please use earlier versions. You can convert them using the pClamp software. [p131c_0011_sw6-10.abf](https://github.com/dervinism/minis/blob/main/p131c_0011_sw6-10.abf) is an example of an ABF file that was used for testing the software.

**References**\
Dervinis, M, Major, G (2022) bioRxiv 2022.03.20.485046; doi: https://doi.org/10.1101/2022.03.20.485046
