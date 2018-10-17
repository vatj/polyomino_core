# Project Title

Base for polyomino model development

## Getting Started

Some minor diffulties on submodule cloning etc, google is probably more helpful at this point...

### Extending the polyomino core model

Need to implement two basic functions
  
Interaction matrix to determine if two interfaces interact
```
bool InteractionMatrix(const typename face_1,const typename face_2);

```

And if they **do** interact, how strong is their interaction
```
double BindingStrength(const interface_type face_1,const interface_type face_2);
```

## Building
The makefile in your extending directory should include a rule to build the submodule
