N-body project: Simulating a cluster of particles

Computational Astrophysics, UvA, 2023/2024

May 6, 2024

1  Introduction

In this project you will simulate a cluster of gravitationally interacting particles using an N-body code that you will write. To properly simulate a cluster we need to have a large number of particles, but as you know by now, the classic n-body problem scales as O(n2), making it extremely computationally expensive for large n. Therefore, in this project, you will implement the [Barnes-Hut algorithm, a ](#_page0_x0.00_y792.00)technique that approximates the N-body problem by grouping distant particles into one aggregated mass. Using this algorithm to simulate a large number of particles, the computational cost now scales as O(n logn), making it highly favourable for large-scale simulations, such as galaxies, clusters, and smoothed particle hydrodynamics (SPH).

2  The N-body problem

For a system of N particles, the force on a particle i from the gravitational attraction of all other particles is given by

¨r = G N m rj − ri (1)

i j − r |3

j=1,j=i |rj i

where r is the position vector, m is the mass, and G is the gravitational constant. In numerical simulations, particles may occasionally get close enough for the distance |rj − ri| ≈ 0, resulting in the force going to infinity. Consequently, a smoothing length is often used to circumvent this issue, which is simply a small number added to the distance calculation to avoid division by zero. Denoting the smoothing length by ε, the equation of motion becomes

N r − r

j i

¨ri = G j=1,j=i mj |rj − ri|3 + ε (2)

3  Simulating using direct summation

Simulating an N-body system of N particles using a direct summation routine simply involves, at each timestep, summing over all particle pairs and calculating the mutual forces, and using this to update the positions and velocities using a numerical scheme such as RK4, Velocity-Verlet, Leapfrog, Euler, etc... As already mentioned, this techniques is effective for small N systems, but quickly slows down the simulation in larger systems.

4  The Barnes-Hut algorithm

The Barnes-Hut (BH) algorithm is a tree-based code, in which the computational grid is recursively divided up into smaller and smaller sections, until each section only contains one or zero particles. When calculating the force on one particle, you can then approximate the force from far-away bodies as<a name="_page0_x0.00_y792.00"></a> the force from a body with the mass equal to the total mass of the distant bodies, and position equal to the center of mass of the group. This means that only particle interactions that are close together are calculated using the direct summation routine, and the rest are approximated as a few, massive bodies. A visual example can be seen in Figure (1).

![](Aspose.Words.d803b2cf-491c-4bed-b2e1-a49aff8ad5c1.001.jpeg)

Figure 1:<a name="_page1_x135.04_y385.66"></a> Example of a computational grid divided into smaller segments using the BH algorithm.

To calculate the force on a particle with position ri and mass mi from a specific node:

1. Calculate the distance between the particle and the center of mass of the node: d = |rcom − ri|
1. If the quantity L/r < θ, where L is the width of the bounding box of the node and θ is a threshold value, then calculate the force using the total mass and center of mass of the node, effectively approximating the force from all the nested particles in the node. If L/r > θ, then the total force becomes a sum over all the children of the node.

The parameter θ determines the resolution of the algorithm. A θ of 0 means that you always sum over all the particles when calculating the force, which reduces the system to a direct-summation routine, while a large θ increases the speed but reduces the accuracy.

BA is a tree-based algorithm, meaning that the grid and the particles are represented as nodes on a tree. Each node can have ”children”, or branches to other nodes. A branch represent a subdivision of a specific subgrid (or boundary box). Each node has a bounding box, denoting the square subdo- main that the node represents. The nodes also have a mass and position property, which defines the mass of all the particles in its nested children and their center of mass position. In this way, the grid can recursively be broken down into smaller and smaller sections by gradually going ”deeper” into the tree. The algorithm goes through all the particles in the simulation, places them in the tree through the process of subdivision until each subgrid contains either one or zero particles. The ”root” of the tree is the node whose children are the four quadrants of the computational grid (or eight octants in 3D), and whose boundary box represents the full grid. The procedure for constructing the tree can be described generally as:

1. Initialize the root of the tree as a node with a bounding box equal to the size of the full grid.
1. Define an insertion procedure for adding a particle to a given node.

[^1]3. Start adding particles to the root node by calling the insertion method. When adding a new

particle, there are three scenarios that can occur:

1) The node contains no particles: In this case, the total mass and center of mass of the node is set to the mass and position of the particle. The node only contains one particle now, so its total mass and position is simply equal to that of the particle.
1) The node already contains a particle, but it has not yet been subdivided: If the node already contains one particle, the node has to be subdivided into smaller subdomains. Initialize the children of the node as a list containing four uninitialized elements. These represents the four quadrants of the bounding box. Then, calculate the quadrant in which the already- occupying particle resides. A new node is then initialized and placed into the list of children, at the index that corresponds to the quadrant, and with a bounding box equal to 1/4th the size of the parent bounding box (its specific boundaries are defined by the quadrant). The previously occupying particle is added to the new node. Repeat the same procedure with the new particle. Then, update the center of mass of the node with the mass and position of the new particle.
1) The node already contains a particle, and it has already been subdivided: The quadrant of the particle is calculated. If there is already a child node at that quadrant, the particle is inserted into this node by calling the insertion method on the child with the particle, otherwise a new node is initialized and the particle is inserted.
   4. Once all the particles have been inserted into the root node, the tree has been constructed, and can then be used to approximate the forces on each particle.
5  Initial conditions
1. Positions

We will assume the density of the star cluster follows a [Plummer model, in ](#_page0_x0.00_y792.00)which the density of stars at a given radius from the centre r is given by Plummer 1911:

3M r2 −5/2

ρ(r) = 1 + (3)

4πa[^2] a2

where M is the total mass of all the stars in the cluster, a is a parameter that determines the size of the core of the cluster. The virial radius is given by

16

r = a. (4)

V 3π

In this project we will be working in [N-body units, whic](#_page0_x0.00_y792.00)h means that the gravitational constant G, the total mass M, and the virial radius rV are all equal to 1.

To generate particles with positions that follow the radial density of the Plummer model, you first sample the inverse probability density function of ρ(r), which is given by

ρ¯(r) = rV~~ , (5)![](Aspose.Words.d803b2cf-491c-4bed-b2e1-a49aff8ad5c1.002.png)

(p− 2 − 1)

3

2. Velocities

For the initial velocities of the particles, there are several options:

1. Random velocities

The particles can be given a random initial velocity equal to

vx = v ∗n1

vy = v ∗n2

√ ![](Aspose.Words.d803b2cf-491c-4bed-b2e1-a49aff8ad5c1.003.png)

where v = 12 2 is the velocity dispersion in N-body units, n1 and n2 are random numbers sampled from the normal distribution with a mean of zero and standard deviation of 1.

2. Rotational velocity

If we want the particles to rotate around the center of the cluster, we can make sure their initial velocity is equal to the tangential rotational velocity. To do this, we simply multiply vx by −sin(θ), and vy by cos(θ), where θ = arctan2(y,x) is the angular position of the particle.

3. Differential rotational velocity

If we want the cluster to move with differential rotation, similar to how galaxies rotate (bodies closer

to the center rotate faster than bodies further out), then the velocity magnitude at a radius r can be calculated using

GM (< r)![](Aspose.Words.d803b2cf-491c-4bed-b2e1-a49aff8ad5c1.004.png)

v(r) = ,

r

where M(< r) is the total mass of all particles with radii smaller than r.

4. Other

You are free to experiment with the initial velocities of the particles. Try combinations of the different methods mentioned above, or try coming up with your own.

6  Tasks
- Implement the Barnes-Hut algorithm in 2D. You may use the given skeleton code as a starting point, but you are also free to write it from scratch. If you want an extra challenge, you may try to implement the 3D version of the algorithm.
- Write a function that initializes a given number of particles using the initial conditions described above.
- Write functions to solve the n-body problem using both direct summation and your Barnes-Hut algorithm. You must write your own integrator, but you are free to choose the method (Euler, Runge-Kutta, Leapfrocy, Velocity-Verlet). Implement both adaptive and fixed timestepping.
- Ensure your code is working correctly by doing the sanity checks.
- Write a method for visualizing your simulations.
7  Sanity checks

With the initial conditions given in Table (1),[ a ](#_page4_x286.08_y123.59)boundary box of size [−0.5,0.5]× [−0.5,0.5], and with θ = 0, you should get the same forces as shown in Table (1) [and](#_page4_x286.08_y123.59) your tree grid should look like Figure ([7).](#_page4_x286.08_y123.59)

|Particle ID x|y||m|<p>F</p><p>x</p>|<p>F</p><p>y</p>|
| - | - | :- | - | - | - |
|1 2 3 4|-0.2 -0.1 0.25 0.25|0\.3 0.25 0.25 -0.25|0\.25 0.25 0.25 0.25|1\.44611 -0.46618 -0.75599 -0.22394|-1.72683 1.56051 -0.34568 0.512|

Table 1:<a name="_page4_x286.08_y123.59"></a> Sanity check table

![](Aspose.Words.d803b2cf-491c-4bed-b2e1-a49aff8ad5c1.005.jpeg)

Figure 2: Sanity check grid

8  Goals

Use your code to simulate a large number of particles. Start with N = 20 and gradually increase the number of particles. Compare the performance of your direct summation routine, and experiment with the value of the θ parameter. Try different initial conditions, and see if you can simulate, for example, a stable cluster, or a galaxy-like object. Plot the total energy of the system over time. Try also changing the initial conditions by adding, for example, a very massive particle in the center, or by shooting a massive particle through the cluster. Be creative! Also, think about the following questions:

- How does the performance of the direct summation technique compare to the Barnes-Hut algo- rithm? Try with different number of particles.
- BH is an approximation algorithm, which means you are losing some accuracy when compared to direct summation. How accurate is the algorithm, how does it relate to θ. How can you compare the accuracy?
- Is the total energy conserved? Why or why not? How does the energy relate to your choice of integrator and/or time step size.
- How many particles can you include in your tree before you run into issues? Can you think of why these issues occur? Can you think of a way to get around them?
- Can you think of ways to optimize this code?
9  Presenting your best simulations

By playing around with initial conditions and simulation parameters, find at least one combination that gives you an interesting result that you would like to present, but if you have multiple exciting simulations then don’t hesitate to show them.
7

[^1]: where p is a number in the range (0,1), generated by sampling a uniform distribution. Once a radius r has been calculated, the cartesian coordinates can be generated by placing the particle in a random position on a circle. This is done by:

    x = r sin(arccos(2p0 − 1))cos(2πp1)

    y = r sin(arccos(2p0 − 1))sin(2πp1) where p0, and p1 are randomly generated numbers in the range (0,1).
[^2]: 
