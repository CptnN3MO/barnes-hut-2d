import numpy as np

class QuadTreeNode:
    """
    A class used to represent the node of a quadtree. 
    A node can represent an individual particle, or a section of the computational grid
    with mass equal to the total mass of all nested particles.

    Methods
    -------
    update_center_of_mass(self, position, mass)
        Update the total mass and COM of the node with the given position and mass.
    
    insert(self, position, mass)
        Insert a body with the given position and mass. If node is already occupied,
        then subdivide into four children, each representing a quadrant within the boundary
        box of the node.

    compute_forces(self, particle, theta, G, eps)
        Compute the forces on the given particle from all the nested particles in the node.
    """
    def __init__(self, x_min, x_max, y_min, y_max):
        """
        TO-DO: initialize the boundary box
        """

        self.bbox = ...
        self.children = None
        self.total_mass = None
        self.center_of_mass = np.zeros(2)


    def update_center_of_mass(self, position, mass):
        """
        TO-DO: write a code for updating the center of mass of the node with the 
        position and mass
        """


    def insert(self, position, mass):
        """
        TO-DO: write a code that adds the given position and mass to the node
        """
        
        if self.total_mass is None:
            # Node is unoccupied: update total mass and COM

        elif self.children is None:
            # Node is occupied, but not yet subdivided: subdivide and distribute old and new particle
            self.children = [None, None, None, None]

        else:
            # Node is occupied and already subdivided: distribute new particle into quadrant

    def subdivide(self, quadrant):
        """
        TO-DO: write a code that adds a new, uninitilized node to the child at the given quadrant
        """

    def compute_forces(self, particle, theta, G, eps):
               """
        TO-DO: write a code that calculates and returns
        the force on a particle from the nested particles in the node.
        The code should check if:
            - the node is uninitialized
            - the particle is the node itself 
            - whether to use the total mass of the node, or a direct summation
              of the nested particles.

        """

        

class QuadTree:
    """
    A class that represents the full tree, which contains all the nodes
    and particles.

    Attributes
    ----------
    particles : array
        An array containing the particles in the simulation. Each particle should have a
        position, velocity, and mass attribute.
    theta : float
        The resolution parameter of the Barnes-Hut algorithm.
    root : QuadTreeNode
        The root node of the tree.

    Methods
    -------
    update_center_of_mass(self, position, mass)
        Update the total mass and COM of the node with the given position and mass.
    
    insert(self)
        Insert and distribute all the particles into the tree.

    """
    def __init__(self, particles, theta):
        """
        TO-DO: initialize the root node and set the particles and theta attributes
        """
        self.root = ...

    def insert(self):
        """
        TO-DO: insert the particles one by one into the root of the tree
        """        
        
        

def get_quadrant(bbox, position):
    """
    TO-DO: write a code that checks which quadrant of the boundary box 
    the given position belongs
    """


def initialize_particles(N):
    """
    TO-DO: write a function that initialized N particles with
    the initial conditions given in the project description.
    """