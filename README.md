# minis

'minis' is a software for electrophysiological data analysis.
It allows you to:\
(1) detect and analyse spontaneous postsynaptic potetnials/currents in whole-cell patch clamp recordings (Dervinis and Major, 2022);\
(2) estimate the quantal size in the central synapse (Dervinis and Major, in preparation).

Download the trial version of the software here: https://github.com/dervinism/minis/releases  \
If downloading the executable is blocked by your web browser, make sure that the executable file is kept and not discarded. To install the trial version on your PC, please follow the installer instructions. As part of the installation process, you will have to install Matlab Runtime. The trial version comes with inbuilt data, so you can try out running the software in various modes: Detect, Detect and Compare, Automatic Distribution Fitting, and Simulate. You will not be able to load your own data.

It is a proprietary software that is currently provided free of charge. If you would like to get a software copy for your use, please email Martynas Dervinis (martynas.dervinis@gmail.com). You will need to provide the serial number of your hard-drive (Windows) or your computing device serial number (macOS) or your BIOS serial number (Linux). 'minis' has been fully written in Matlab (Mathworks) and is distributed as an application programming interface in the form of a packaged Matlab app or a Python (Python Software Foundation) package or as a compiled standalone desktop application with a graphical user interface.

**Instructions on how to get the required serial number**\
On Windows open your command prompt and type in ```wmic diskdrive get model,serialnumber```\
Any internal disk will work.

On macOS open your terminal and type in ```system_profiler SPHardwareDataType | grep Serial```\

On Linux open your terminal and type in ```dmesg | grep "DMI:" | cut -c "6-" | cut -d "," -f "2"```\

**Installation instructions**\
When you get hold of the fully functional software, you can install the standalone version of the software by running the installer inside minisStandalone folder. Simply follow the installer instructions. Uninstall the minis the same way you unistall any other regular app.

In order to install minis as a Matlab packaged app, double click minisMatlab.mlappinstall inside minisMatlab folder and follow instructions inside Matlab. To uninstall, navigate to Matlab Apps section, right-click minisMatlab under MY APPS subsection, and uninstall it.

If you intend to use minis as a Python package, you can install the setup.py file located inside the minisPy folder. Simply open GettingStarted.html file and follow the instructions outlined inside the file. However, prior to installing minisPy package you have to make sure that all required Python package dependencies are also installed. These are:\
python=3.7\
matlab\
pyabf

If you run Python on Anaconda, you can install minisPy in a separate environment by following these steps (the last line is optional)
```
conda create --name minis-env python=3.7
conda activate minis-env
pip install matlab
pip install pyabf
cd minisPy
python setup.py install
conda install spyder-kernels
```

Software user documentation file minis_documentation.pdf is available for a detailed explanation of how to use the software graphical user interface. Examples on how to use programming interfaces in Matlab and Python are given in testMatlab.m, testMatlab_preload.m, testPython.py, and testPython_preload.py files.

**References**\
Dervinis, M, Major, G (2022) bioRxiv 2022.03.20.485046; doi: https://doi.org/10.1101/2022.03.20.485046
