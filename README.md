# MLFMA

## Design Principles 

- Dimensions are normalized to wavelength.

## Run 

Edit input parameters and run the application using run script.

```bash
export LEAFBOX: Leaf-level box size (in wavelengths)
export MODELFILE: Model file path
export SCALE: Scaling factor

#RUN COMMAND COMES HERE
```

## MODEL FILE BINARY FORMAT

int16 numnode
fp16 x1
fp16 y1
fp16 z1
fp16 x2
fp16 y2
fp16 z2
...
