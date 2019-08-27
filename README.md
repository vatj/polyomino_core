## Forked from ASleonard repository

## Polyomino self-assembly core
[![Codacy Badge](https://api.codacy.com/project/badge/Grade/821948e0aa6c4984b448536609d64e85)](https://www.codacy.com/app/ASLeonard/polyomino_core?utm_source=github.com&amp;utm_medium=referral&amp;utm_content=ASLeonard/polyomino_core&amp;utm_campaign=Badge_Grade)
![License Badge](https://img.shields.io/github/license/ASLeonard/polyomino_core.svg?style=flat)

Base repository for the polyomino model development. This is mainly useful for extending the functionality provided here in a new repository, using this one as a submodule. For an example, see an example [here](https://github.com/ASLeonard/polyomino_interfaces).

### Extending the polyomino core model
Including this repository as a submodule allows for rapid prototyping of the polyomino model.

Generic assembly can be quickly implemented using existing functions. Extending the model **requires**
  - extend base assembly (using the curiously recursive template pattern)
```cpp
class NewAssemblyModel : public PolyominoAssembly<NewAssemblyModel>
```
  - determine the interface type (e.g. integer, string, something more complex, etc.)
```cpp
using subunit_type = ...
using Genotype = std::vector<subunit_type>
```
  - implement the interaction matrix that returns the interaction strength between two interfaces
```cpp
double NewAssemblyModel::InteractionMatrix(subunit A, subunit B)
```
  - optional features to implement may include more specific methods
    - mutations
    - random initilising
    - reversing alignment
    
### Using the core model
The main calls to the core model are to the assembly and classification functions
  - get interaction edges
  ```cpp
  const auto edges = NewAssemblyModel::GetActiveInterfaces(genotype)
  ```
  - assemble polyomino (and track edges used in assembly)
  ```cpp
  auto [assembly_information,interacting_indices] = NewAssemblyModel::AssemblePolyomino(edges)
  ```
  - get phenotype from the raw assembly information
  ```cpp
  Phenotype phen=GetPhenotypeFromGrid(assembly_information)
  ```
  - potentially classify phenotype using a runtime filled table
  
  ```cpp
  Phenotype_ID pid = phenotypeTable->GetPhenotypeID(phen)
  ```
  
More robust examples available in the implementation of _polyomino\_interfaces_ model

### Building
The core library is header only, so need to add an include path to whatever building process used
```make
-Ipath/to/library/polyomino_core
```
 Some of the code makes use of the relatively modern c++17, and so a relatively up to date compiler may be needed. Again a more full example of a makefile used in compiling may be found in the _polyomino\_iterfaces_ repository.

### Visuals
There are some skeleton methods in python3 to plot polyominoes and more elaborate interactive plots to show evolutionary transitions between polyominoes, but these are still under design and development.
