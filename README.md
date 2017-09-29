#Three Body Problem and Planetary Motion Simulations

Planetary motion simulations use Runge Kutta method and Euler's method to calculate the trajectory of planetary object in 2D and 3D. 

Additional Lagrangian points script allows the calculation and visualisation of the lagrange points in selected 3 body configurations.

## Simulation Details

### Deployment

Units in calculations for three body script are in Au (Astronomical Units), Days and kg. 

### Constraints

Three body calcuations induce error that grows exponentially to the number of steps and will produce wildy innacurate results if velocities, accelerations, positions are too high. Values are optimised for objects in Solar System scale