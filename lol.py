import numpy as np
import matplotlib.pyplot as plt

class Particle:
    def __init__(self, position, velocity, mass):
        self.position = np.array(position)
        self.velocity = np.array(velocity)
        self.mass = mass

class QuadTreeNode:
    def __init__(self, x_min, x_max, y_min, y_max):
        self.bbox = (x_min, x_max, y_min, y_max)  
        self.children = [None, None, None, None]  
        self.total_mass = 0
        self.center_of_mass = np.zeros(2)
        self.particle = None

    def update_center_of_mass(self, position, mass):
        # Update the center of mass and total mass
        if self.total_mass == 0:
            self.center_of_mass = position * mass
        else:
            self.center_of_mass = (self.center_of_mass * self.total_mass + position * mass) / (self.total_mass + mass)
        self.total_mass += mass

    def insert(self, position, mass):
        # Insert a particle into the quadtree
        if self.total_mass == 0:
            self.particle = Particle(position, [0, 0], mass)
            self.update_center_of_mass(position, mass)
        elif self.children[0] is None:
            self.subdivide()
            if self.particle is not None:
                old_particle = self.particle
                self.particle = None
                quadrant = get_quadrant(self.bbox, old_particle.position)
                self.children[quadrant].insert(old_particle.position, old_particle.mass)
            quadrant = get_quadrant(self.bbox, position)
            self.children[quadrant].insert(position, mass)
        else:
            quadrant = get_quadrant(self.bbox, position)
            self.children[quadrant].insert(position, mass)
        self.update_center_of_mass(position, mass)

    def subdivide(self):
        # Subdivide the current node into four quadrants
        x_min, x_max, y_min, y_max = self.bbox
        mid_x = (x_min + x_max) / 2
        mid_y = (y_min + y_max) / 2
        self.children = [
            QuadTreeNode(x_min, mid_x, y_min, mid_y),
            QuadTreeNode(mid_x, x_max, y_min, mid_y),
            QuadTreeNode(x_min, mid_x, mid_y, y_max),
            QuadTreeNode(mid_x, x_max, mid_y, y_max)
        ]

    def compute_forces(self, particle, theta, G, eps):
        # Compute the gravitational forces on a particle from the quadtree
        if self.total_mass == 0:
            return np.zeros(2)
        
        distance = np.linalg.norm(self.center_of_mass - particle.position)
        if (self.children == [None, None, None, None]) or (max(self.bbox[1] - self.bbox[0], self.bbox[3] - self.bbox[2]) / distance < theta):
            if np.array_equal(particle.position, self.center_of_mass):
                return np.zeros(2)
            force_direction = (self.center_of_mass - particle.position) / (distance**3 + eps)
            return G * self.total_mass * particle.mass * force_direction
        else:
            force = np.zeros(2)
            for child in self.children:
                if child is not None:
                    force += child.compute_forces(particle, theta, G, eps)
            return force

class QuadTree:
    def __init__(self, particles, theta):
        x_min, x_max = [-0.5,0.5]
        y_min, y_max = [-0.5,0.5]
        self.root = QuadTreeNode(x_min, x_max, y_min, y_max)
        self.particles = particles
        self.theta = theta

    def insert(self):
        for particle in self.particles:
            self.root.insert(particle.position, particle.mass)
def get_quadrant(bbox, position):
    # Determine which quadrant a position belongs to within a boundary box
    x_min, x_max, y_min, y_max = bbox
    mid_x = (x_min + x_max) / 2
    mid_y = (y_min + y_max) / 2
    x, y = position
    if x <= mid_x and y <= mid_y:
        return 0
    elif x > mid_x and y <= mid_y:
        return 1
    elif x <= mid_x and y > mid_y:
        return 2
    else:
        return 3

# def initialize_particles():
#     # Initialize particles with predefined positions and masses
#     particles = [
#         Particle([-0.2, 0.3], [0, 0], 0.25),
#         Particle([-0.1, 0.15], [0, 0], 0.25),
#         Particle([0.25, 0.25], [0, 0], 0.25),
#         Particle([0.25, -0.25], [0, 0], 0.25)
#     ]
#     return particles
def initialize_particles(N):
    particles = []
    for _ in range(N):
        p0 = np.random.random()
        p1 = np.random.random()
        r = (16 / 3 * np.pi)**(1/3) * np.sqrt(p0**(-2/3) - 1)
        theta = 2 * np.pi * p1
        phi = np.arccos(2 * np.random.random() - 1)
        x = r * np.sin(phi) * np.cos(theta)
        y = r * np.sin(phi) * np.sin(theta)
        mass = 1.0 / N
        particles.append(Particle(np.array([x, y]), np.zeros(2), mass))
    return particles

def velocity_verlet(particles, dt, theta, G, eps, positions):
    tree = QuadTree(particles, theta)
    tree.insert()
    
    for particle in particles:
        force = tree.root.compute_forces(particle, theta, G, eps)
        acceleration = force / particle.mass
        
        particle.position += particle.velocity * dt + 0.5 * acceleration * dt**2
        
        tree = QuadTree(particles, theta)
        tree.insert()
        
        new_force = tree.root.compute_forces(particle, theta, G, eps)
        new_acceleration = new_force / particle.mass
        
        particle.velocity += 0.5 * (acceleration + new_acceleration) * dt
    
    # Append new positions to the history
    for i, particle in enumerate(particles):
        positions[i].append((particle.position[0], particle.position[1]))
    
    return tree, positions


def draw_quadtree(node, ax):
    # Recursively draw the quadtree structure on the plot
    if node is None:
        return
    x_min, x_max, y_min, y_max = node.bbox
    rect = plt.Rectangle((x_min, y_min), x_max - x_min, y_max - y_min, fill=False, edgecolor='blue', linewidth=0.5)
    ax.add_patch(rect)
    for child in node.children:
        if child is not None:
            draw_quadtree(child, ax)
import numpy as np
import matplotlib.pyplot as plt
from matplotlib.animation import FuncAnimation

# Function to draw quadtree
def draw_quadtree(node, ax):
    if node is None:
        return
    x_min, x_max, y_min, y_max = node.bbox
    rect = plt.Rectangle((x_min, y_min), x_max - x_min, y_max - y_min, fill=False, edgecolor='blue', linewidth=0.5)
    ax.add_patch(rect)
    for child in node.children:
        if child is not None:
            draw_quadtree(child, ax)

# Initialize particles
N = 15
particles = initialize_particles(N)
theta = 0.5
dt = 0.05  # Time step
steps = 100  # Number of steps to simulate

positions = [[(p.position[0], p.position[1])] for p in particles]  # Initialize positions list for each particle

# Setup the plot
fig, ax = plt.subplots()
ax.set_xlim(-2, 2)
ax.set_ylim(-2, 2)
ax.set_facecolor('lightgrey')

# Animation function
def animate(step):
    global particles, positions, ax
    
    ax.clear()
    ax.set_xlim(-2, 2)
    ax.set_ylim(-2, 2)
    ax.set_facecolor('lightgrey')
    
    # Update particles and get the updated tree
    tree, positions = velocity_verlet(particles, dt, theta, G=1.0, eps=0.001, positions=positions)
    
    # Plot the traces
    for pos in positions:
        x_vals = [p[0] for p in pos]
        y_vals = [p[1] for p in pos]
        ax.plot(x_vals, y_vals, 'r-', alpha=0.5)
    
    # Plot the particles
    scatter = ax.scatter([p.position[0] for p in particles], [p.position[1] for p in particles], c='red', s=10)
    
    # Draw the quadtree
    draw_quadtree(tree.root, ax)
    
    ax.set_title(f"Step {step}")

# Create animation
anim = FuncAnimation(fig, animate, frames=steps, interval=50)

# Save animation
anim.save('nbody_simulation.mp4', writer='ffmpeg')

plt.show()
