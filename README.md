# MESO Modern Physical Simulator

## Getting started

##### Install Environment
- Visual Studio 2022
- CUDA 11.6
- Paraview 5.10.1-MPI
- Microsoft MPI (should be installed by default)

##### Install xmake

1. Press `Ctrl+R` and type `powershell`, then enter to open `Powershell` terminal. Note that you don't need to run in admin mode.
2. Type:

        Invoke-Expression (Invoke-Webrequest 'https://xmake.io/psget.text' -UseBasicParsing).Content
3. Type `M` in installing (which means it will install `xmake` for the current user).
4. (Optional, not necessary) Update xmake to dev version

         xmake update dev

##### Run Tests

      $ python run_tests.py


You should pass all tests.

##### Generate Solution of a Fluid Simnulator
We use `fluid_euler` as an example here.

      $ python make_project.py fluid_euler

##### Build and Run a Cavity Flow
Open `\build\fluid_euler\vsxmake2022\fluid_euler.sln` and compile in Release mode.
Navigate to `\bin\fluid_euler\windows\x64\release` and copy `docs\fluid_euler\cavity.json`  to the current folder. Run:

      $ fluid_euler.exe cavity.json

An output folder named `cavity` with `.vst` files will be generated.

##### Visualize Results with Paraview
- Open Paraview
- Navigate to the folder containing `.vst` files and open the files.
- Click on the eye symbol to show the results.
- Suggestions on modifying settings:
  -  Coloring: velocity
  - Right click on the file to add filter (for example `glyph`)
  - Orientation array: velocity
  - Masking: all points (continous display of points)
  - Click `apply` to save the settings
- Click &rarr; in the toolbar at top of screen to play the animation.
- `ctrl+s` to save the result.
