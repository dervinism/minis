# minis

'minis' is a software for electrophysiological data analysis.
It allows you to:\
(1) detect and analyse spontaneous postsynaptic potetnials/currents in whole-cell patch clamp recordings (Dervinis and Major, 2022);\
(2) estimate the quantal size in the central synapse (Dervinis and Major, in preparation).

Download the trial version of the software here: https://github.com/dervinism/minis/releases/tag/v1.1.0  \
If downloading the executable is blocked by your web browser, make sure that the executable file is kept and not discarded. To install the trial version on your computer, please follow the installer (Windows: minisTrialInstaller_web.exe, Linux: minisTrialInstaller_web.install, MacOS: minisTrialInstaller_web.app) instructions. As part of the installation process, you will have to install Matlab R2022a Runtime. The trial version comes with inbuilt data, so you can try out running the software in various modes: Detect, Detect and Compare, Automatic Distribution Fitting, and Simulate. You will not be able to load your own data.

It is a proprietary software that is currently provided free of charge. If you would like to get a software copy for your use, please email Martynas Dervinis (martynas.dervinis@gmail.com). You will need to provide the serial number of your hard-drive (Windows) or your computing device serial number (MacOS) or your BIOS serial number (Linux). In the email please also indicate which operating system you are using. 'minis' has been fully written in Matlab (Mathworks) and is distributed as an application programming interface in the form of a packaged Matlab app or a Python (Python Software Foundation) package or as a compiled standalone desktop application with a graphical user interface.

**Instructions on how to get the required serial number**\
On Windows open your command prompt and type in ```wmic diskdrive get model,serialnumber```\
Any internal disk will work.

On MacOS open your terminal and type in ```system_profiler SPHardwareDataType | grep Serial```

On Linux open your terminal and type in ```sudo dmesg | grep "DMI:" | cut -c "6-" | cut -d "," -f "2"``` \
Linux terminal will give you a BIOS number.

**Installation instructions: Standalone desktop app**\
When you get hold of the fully functional software, you can install the standalone version of the software by running the installer inside minisStandalone folder. Simply follow the installer instructions. Uninstall the minis the same way you unistall any other regular app.

**Installation instructions: Matlab packaged app**\
In order to install minis as a Matlab packaged app, double click minisMatlab.mlappinstall inside minisMatlab folder and follow instructions inside Matlab. To uninstall, navigate to Matlab Apps section, right-click minisMatlab under MY APPS subsection, and uninstall it.

**Installation instructions: Python package on Windows**\
Navigate to the minisPy folder and install the setup.py file (follow the instructions outlined inside the GettingStarted.html file; make sure you install Matlab R2022a Runtime before you complete the rest of the installation steps). Prior to installing minisPy package you have to make sure that all required Python package dependencies are also installed. These are:\
python=3.9\
matlab\
pyabf

Execute the commands below in the terminal (if you have correct Matlab Runtime and Python version):
```
pip install matlab
pip install pyabf
cd minisPy
python setup.py install
```

**Installation instructions: Python package on Linux**\
To install follow these steps:
1. Install Python 3.9 (in a separate environment if needed).
2. Install Matlab R2022a Runtime.
3. Update Matlab Runtime path by adding the lines below to your ~./bashrc file:
```
export LD_LIBRARY_PATH="${LD_LIBRARY_PATH:+${LD_LIBRARY_PATH}:}\
/usr/local/MATLAB/MATLAB_Runtime/v912/runtime/glnxa64:\
/usr/local/MATLAB/MATLAB_Runtime/v912/bin/glnxa64:\
/usr/local/MATLAB/MATLAB_Runtime/v912/sys/os/glnxa64:\
/usr/local/MATLAB/MATLAB_Runtime/v912/extern/bin/glnxa64"
export LD_PRELOAD="${LD_PRELOAD:+${LD_PRELOAD}:}\
/usr/local/MATLAB/MATLAB_Runtime/v912/bin/glnxa64/glibc-2.17_shim.so"
```
&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp; Make sure you also execute these commands in your terminal. \
&nbsp;&nbsp;&nbsp; 4. Install minis and its dependencies by executing the lines below in your terminal:
```
pip install matlab
pip install pyabf
cd minisPy
python setup.py install
```

**Installation instructions: Python package on macOS**\
To install follow these steps:
1. Install Python 3.9 (in a separate environment if needed).
2. Install Matlab R2022a Runtime
3. Update Matlab Runtime path by adding the lines below to your ~./bash_profile and ~./zshenv files:
```
export DYLD_LIBRARY_PATH="${DYLD_LIBRARY_PATH:+${DYLD_LIBRARY_PATH}:}\
/Applications/MATLAB/MATLAB_Runtime/v912/runtime/maci64:\
/Applications/MATLAB/MATLAB_Runtime/v912/bin/maci64:\
/Applications/MATLAB/MATLAB_Runtime/v912/sys/os/maci64:\
/Applications/MATLAB/MATLAB_Runtime/v912/extern/bin/maci64"
```
&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp; Make sure you also execute this command in your terminal. \
&nbsp;&nbsp;&nbsp; 4. Change python executable associated with mwpython by executing the lines below in your terminal:
```
echo 'export PYTHONHOME=<python3 installation directory>' >> ~/.zshenv
source ~/.zshrc
echo 'export PYTHONHOME=<python3 installation directory>' >> ~/.bash_profile
source ~/.bash_profile
```
&emsp;&emsp;In R2022a, PYTHONHOME needs to be set to an actual directory that has a bin subdirectory that contains the python3.9 executable.
&emsp;&emsp;Typically the correct directory is /Library/Frameworks/Python.framework/Versions/3.9.

5. Install minis and its dependencies by executing the lines below in your terminal:
```
pip3 install matlab
pip3 install pyabf
cd minisPy
python3 setup.py install
```
6. Use mwpython to run any code that involves loading Python packages compiled in Matlab. For example,
```
/Applications/MATLAB/MATLAB_Runtime/v912/bin/mwpython testPython.py
```

**Documentation**\
Software user documentation file [minis_documentation.pdf](https://github.com/dervinism/minis/blob/main/minis_documentation.pdf) is available for a detailed explanation of how to use the software graphical user interface. Examples on how to use programming interfaces in Matlab and Python are given in [testMatlab.m](https://github.com/dervinism/minis/blob/main/testMatlab.m), [testMatlab_preload.m](https://github.com/dervinism/minis/blob/main/testMatlab_preload.m), [testPython.py](https://github.com/dervinism/minis/blob/main/testPython.py), and [testPython_preload.py](https://github.com/dervinism/minis/blob/main/testPython_preload.py) files.

**Input files**\
The software does not work with version 2 ABF files. Please use earlier versions. You can convert them using the pClamp software.

**References**\
Dervinis, M, Major, G (2022) bioRxiv 2022.03.20.485046; doi: https://doi.org/10.1101/2022.03.20.485046
