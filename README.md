# Duplex Modern Physical Simulator

## Getting started

##### Install Environment
- Visual Studio 2022
- CUDA 11.6

##### Install xmake

1. Press `Ctrl+R` and type `powershell`, then enter to open `Powershell` terminal. Note that you don't need to run in admin mode.
2. Type:

        Invoke-Expression (Invoke-Webrequest 'https://xmake.io/psget.text' -UseBasicParsing).Content
3. Type `M` in installing (which means it will install `xmake` for the current user).
4. (Optional, not necessary) Update xmake to dev version

         xmake update -s github:xq114/xmake#dev

##### Compile Project
        $ python make_project.py _tests_reservoir

##### Build Project
Open `\build\_tests_reservoir\vsxmake2022\_tests_reservoir.sln` and compile.

Note: in the output window, the first command should be something like

      $xmake config -y  -P . -p windows -m debug -a x64 -o "..\..\bin\_tests_reservoir"

If it's

      $xmake config -y  -P . -p windows -m debug -a x64 -o "build"
There is an error. But we don't exactly know why. Report when this happens.      
