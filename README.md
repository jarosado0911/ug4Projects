# ug4Projects
This is where I keep my ug4 project codes and explorations. Below are the specifications of the version of `ughub`

### 📂 Repository Information for ../ughub
-----------------------------------------
Commit Hash:        2143e0d7d900a98f4e997d610062758b626d783d

Version Tag:        v2021.amfs

Creation Date:      2015-10-22 11:33:39 +0200

Last Modified:      2024-06-05 08:08:21 +0200

Total Commits:      87

---
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
