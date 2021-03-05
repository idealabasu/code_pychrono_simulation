# PyChrono
Guides for using PyChrono

## Install
Use conda to instal. Please note that this will instal the latest version.
Only the latest version has both `pardisomkl` and `ChMaterialShellReissner`. Hopefully, they can publish a stable release soon.
```
conda install -c projectchrono/label/develop pychrono
```

## Notes
### Unit
- Chrono is unitless but follows SI standard. Adjust the unit system may improve the simulation.

### FEA
- Constraint for fea node seems to only support `ChNodeFEAxyzrot`.
- If combining FEA and constraints, use `ChSolverPardisoMKL` solver and `ChTimestepperHHT` time stepper. Please also `SetStepControl(False)`.

### Collision
- Decreasing time step yields different dynamics, which may be caused by the collision detection algorithm.
- Playing around `SetDefaultSuggestedEnvelope`, `SetDefaultSuggestedMargin`, and `SetYoungModulus` can improve the simulation.
- Adjusting unit of time may also work. For jumping which happens in a short time, one can use 0.1s or 0.01s.

### Visualization
-  The default visualization engine is Irrlicht. It uses LEFT handed coordinate system. Chrono itself uses RIGHT handed coordinate system.
