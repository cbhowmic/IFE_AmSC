# Unified Schema Outline

```
/
├── schema                         # [str] e.g. OpenPMD+X, IMAS
├── schemaVersion                  # [str] version of the unified schema
├── basePath                       # [str] e.g. "/data"
├── meshesPath                     # [str] default "meshes"
├── particlesPath                  # [str] default "particles"
├── iterationEncoding              # [str] groupBased | variableBased | fileBased
├── iterationFormat                # [str] linked to iterationEncoding
│
├── software/
│   ├── softwareName               # [str] RHINO, Laser, IPM
│   ├── softwareVersion            # [str] Version of software (e.g. version of RHINO that generated data)
│   ├── softwareDescription        # [str] Name of the application or code (e.g. RHINO, Laser, IPM). basePath = softwareName?
│   ├── versionControlSoftware     # [str] e.g. GitHub, GitLab
│   ├── softwareCommit             # [str] e.g. Git commit hash
│   └── softwareDocumentation      # [str] Documentation webpage for reproducibility
│
├── provenance/
│   ├── author                     # [str] person who generated the data, e.g. Holly Flynn for RHINO
│   ├── authorAffiliation          # [str]
│   ├── authorEmail                # [str]
│   ├── creationDate               # [str] YYYY-MM-DD HH:mm:ss tz
│   ├── creationTimeUTC            # [str] (Possibly redundant)
│   ├── inputDirectory             # [str] e.g. path to directory where all RHINO .pkl files live
│   ├── inputFiles                 # [list[str]] e.g. each represents relative path (wrt the inputDirectory) of the data files read by schema
│   ├── originalDataDirectory      # [str] e.g. directory where all the original files (non-compliant to this unified schema) live 	
│   ├── originalDataFiles          # [list[str]]
    ├── originalInputDirectory     # [str] e.g. Directory containing input files that are used to run the RHINO exp
    └── originalInputFiles         # [list[str]] e.g. input files that are used to run the RHINO exp
│
├── system/
│   ├── systemIP                   # [str] (for simulation)
│   ├── systemDescription          # [str] description of the machine used to generate input data 
│   └── comment                    # [str] (optional)
│
├── metadata/                       # Application-level metadata (partially filled per code)
│   │
│   ├── species/                    # Relevant for RHINO, IPM, etc.
│   │   ├── names                   # [list[str]] Species names
│   │   ├── description             # [str] Description of species definitions
│   │   │
│   │   ├── species[0]/             # Example: Tritium
│   │   │   ├── name                # [str]
│   │   │   ├── subsystems          # [list[str]] Subsystems species flows through
│   │   │   └── subsystemConnection # Adjacency matrix of subsystem connectivity
│   │   │
│   │   ├── species[1]/             # Example: Deuterium
│   │   └── ...
│   │
│   ├── subsystems/                 # Subsystem definitions (RHINO and others)
│   │   ├── description             # [str] What subsystems represent
│   │   ├── id                      # [list[int]] Encoded IDs (e.g., 0–20 for RHINO)
│   │   ├── labels                  # [list[str]] Names aligned with id
│   │   ├── mapping                 # [dict] {subsystem_name: id}
│   │   ├── connections             # Adjacency matrix of subsystem connectivity
│   │   └── connectionType          # [str] bidirectional | unidirectional | weighted | non-weighted
│   │
│   └── cameras/                    # Relevant for laser data (details TBD)
│
├── data/                           # Iteration-level data begins here
│   ├── time                        # Snapshot time
│   ├── dt                          # Last timestep used to reach snapshot
│   ├── timeUnitLabel               # [str] e.g., "day"
│   ├── timeUnitSI                  # [float] Seconds per time unit
│   │
│   ├── meshes/                     # “Inventories” for RHINO
│   │   ├── meshesNames             # [list[str]] Names of mesh records
│   │   │
│   │   ├── mesh[0]/                # Example: Tritium mesh
│   │   │   ├── description
│   │   │   ├── geometry
│   │   │   ├── geometryParameters
│   │   │   ├── dataOrder
│   │   │   ├── gridUnitDimension
│   │   │   ├── unitSI
│   │   │   ├── timeOffset
│   │   │   ├── axisLabels          # Subsystems if species represented as meshes
│   │   │   ├── gridGlobalOffset
│   │   │   ├── gridSpacing
│   │   │   ├── gridUnitSI
│   │   │   ├── position
│   │   │   ├── particleList
│   │   │   └── comment
│   │   │
│   │   ├── mesh[1]/                # Example: Deuterium mesh
│   │   └── ...
│   │
│   ├── particles/                  # Particle records (DataFrame-like storage)
│   │   ├── group                   # e.g., Tritium
│   │   ├── id
│   │   ├── position/
│   │   │   ├── x
│   │   │   ├── y
│   │   │   └── z
│   │   └── positionOffset/
│   │       ├── x
│   │       ├── y
│   │       └── z
│   │
│   └── arrays/                     # Third record type (non-spatial semantics; proposed)
│
└── internalMetadata/               # Schema-specific, library-managed metadata
                                    # e.g., "__openPMD_internal"
```


## Proposed OpenPMD Extension (0.17.0)

```
- Hierarchical metadata organization
- Third record type
- Additional in-built attributes
```