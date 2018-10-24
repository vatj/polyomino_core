# Project Title

Base for polyomino model development

## Getting Started

Some minor diffulties on submodule cloning etc, google is probably more helpful at this point...

### Extending the polyomino core model
Generic assembly can be quickly implemented using existing functions.
Extend the model as
  - extend base assembly 
```
class NewAssemblyModel : public PolyominoAssembly<NewAssemblyModel>
```
  - implement the interaction matrix that returns the interaction strength between two interfaces
```
double InteractionMatrix(typename A, typename B)
```
  - implement any specific methods
    - Mutations
    - interaction matrix related calculations

## Building
The core library is header only, so need to add an include path
  - -Ipath/to/library
