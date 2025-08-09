# ug4Projects
This is where I keep my ug4 project codes and explorations. Below are the specifications of the version of `ughub`
I recommending cloning ug4/ughub to your `ROOT` directory and building ug4 there.

```
📂 Repository Information for ../ughub
-----------------------------------------
Commit Hash:        2143e0d7d900a98f4e997d610062758b626d783d
Version Tag:        v2021.amfs
Creation Date:      2015-10-22 11:33:39 +0200
Last Modified:      2024-06-05 08:08:21 +0200
Total Commits:      87
```

I have provided a simple shell script `repo_info.sh` which will print the above information for a repo. Its usage is below

```
./repo_info.sh /path/to/repository
```

---
I have also provided the file structure for my main folder, the main folder contains this repository and `ug4` and `ughub`, the latter two you will need to build separately, by following the directions on the gitHub page for ug4

```
├── ug4
│   ├── .ughub
│   ├── CMakeLists.txt
│   ├── apps
│   ├── bin
│   ├── build
│   ├── externals
│   ├── lib
│   ├── plugins
│   ├── tools
│   └── ugcore
├── ug4Projects
│   ├── LICENSE
│   ├── README.md
│   ├── docs
│   ├── images
│   ├── projects
│   └── repo_info.sh
└── ughub
    ├── .gitignore
    ├── LICENSE
    ├── README.md
    ├── __pycache__
    ├── autocompletions
    ├── docu
    ├── project_templates
    ├── ughub
    ├── ughub.bat
    ├── ughub.py
    ├── ughubHelp.py
    ├── ughubHelpContents.py
    ├── ughubProjectFileGenerator.py
    └── ughubUtil.py
```

You can get the above tree by executing
```
tree -a -L 2 -I '.git|CMake.*'
```

# Install of `ug4` and `ughub`
Please follow the steps for installing `ughub` via the github page for ug4. A few things to note because my installation will most likely be different from yours, but below are some of the notes I have for my installation...
- you will need to execute inside the `ug4` folder
```
../path/to/ughub/ughub init
```
- then execute inside `ug4`
```
../path/to/ughub/ughub install Examples
```
- you will need to make a `build` directory and then `cd` in to it. Once inside the `build` directory execute
```
cmake -DDEBUG=ON -DPARALLEL=OFF -DCMAKE_BUILD_TYPE=Debug -DUSER_LAPACK_LIBRARIES=/usr/lib/x86_64-linux-gnu/lapack/liblapack.so -DUSER_BLAS_LIBRARIES=/usr/lib/x86_64-linux-gnu/libblas.so -DConvectionDiffusion=ON ..
```
- then execute `make -j2`
- the above steps will build `ugshell` inside `ug4/bin/` folder
- still inside `build` execute
```
cmake -DProMesh=ON .. && make -j2
```
- I plan on using ProMesh features in my explorations.
- I am also using codes developed by my gradute school colleagues and doctoral advisor, to include them we execute in `ug4` folder
```
../path/to/ughub/ughub addsource neurobox https://github.com/NeuroBox3D/neurobox-packages.git
```
- then execute
```
../path/to/ughub/ughub install neuro_collection cable_neuron
```
this will install two of the projects to the to ug4 build
- next go back into `build` and execute
```
cmake -Dneuro_collection=ON -Dcable_neuron=ON .. && make -j2
```
