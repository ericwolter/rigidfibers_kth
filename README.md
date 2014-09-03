rigidfibers
===========

Extensions
----------
- Visualization
- Periodicity
- Brownian motion

Open questions?
---------------

OpenCL vs CUDA?
Library for solver?


Closed question
---------------
Q: External force, only gravity? Why different for each fiber?
A: Yes, for now only gravity. Does not have to be different for each fiber

Q: Equation 20, where is term 3/2(M...)?
A: See page 8 -> All external torques are zero

Q: BVector numeric stability?
A: Non issue, was bug in Gmat calulation
