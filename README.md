# Duplex Modern Physical Simulator

## Getting started

### Install Environment

##### Install xmake

1. Press `Ctrl+R` and type `powershell`, then enter to open `Powershell` terminal. Note that you don't need to run in admin mode.
2. Type:

        $ Invoke-Expression (Invoke-Webrequest 'https://xmake.io/psget.text' -UseBasicParsing).Content
3. Type `M` in installing (which means it will install `xmake` for the current user).


        $ λ xmake project -k vsxmake -m "debug,release" -P proj\compilation_tests\

        $  xmake project -k vsxmake -m "Debug,Release,ReleaseWithDebugInfo" -a "x64" -v -F proj\compilation_tests\xmake.lua ..\..\build\compilation_tests

`-P` stands for project.
