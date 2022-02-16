from astropy.coordinates import SkyCoord
import csv
from os import path
import numpy as np
import matplotlib.pyplot as plt
import astropy.units as u 
import astropy.constants as c
import sys
import numpy as np
import matplotlib.animation as manimation
import time
from tqdm import tqdm
import mpl_toolkits.mplot3d.axes3d as p3
#Sources: https://prappleizer.github.io/Tutorials/RK4/RK4_Tutorial.html

class Body:
    
    def __init__(self, pos_vec, vel_vec, mass, name):
        self.pos = pos_vec.cgs.value
        self.vel = vel_vec.cgs.value
        self.mass = mass.cgs.value
        self.name = name

    def return_vector(self):

        return np.concatenate((self.pos, self.vel))

    def return_masses(self):

        return self.mass

    def __sub__(self, other):
        pos = self.pos - other.pos
        return pos


class System:

    def __init__(self, data, step_size=1, G=9.81):
        self.step_size = step_size
        self.acceleration_constant = G * u.m/(u.s**2)
        self.bodies = []

        with open(data) as csvfile:
            reader = csv.reader(csvfile, delimiter=',')
            for row in reader:
                pos_vec = np.array([float(i) for i in row[0].split()])
                vel_vec = np.array([float(i) for i in row[1].split()])
                mass = float(row[2])
                name = str(row[3])
                body = Body(pos_vec * u.AU, vel_vec *u.km/u.s, mass*c.M_earth, name)
                self.bodies.append(body)

        state_vectors = np.array([body.return_vector() for body in self.bodies], dtype=object)
        self.system_state_vector = np.concatenate(state_vectors, dtype=object)
        self.system_mass_vector = np.array([body.return_masses() for body in self.bodies])

    def gravity_evaluator(self, t, y, m):
        N_bodies = int(len(y) / 6)
        solved_vector = np.zeros(y.size)
        for i in range(N_bodies):
            ioffset = i * 6 
            for j in range(N_bodies):
                joffset = j*6
                solved_vector[ioffset] = y[ioffset+3]
                solved_vector[ioffset+1] = y[ioffset+4]
                solved_vector[ioffset+2] = y[ioffset+5]
                if i != j:
                    dx = y[ioffset] - y[joffset]
                    dy = y[ioffset+1] - y[joffset+1]
                    dz = y[ioffset+2] - y[joffset+2] 
                    r = (dx**2+dy**2+dz**2)**0.5
                    ax = (-c.G.cgs * m[j] / r**3) * dx
                    ay = (-c.G.cgs * m[j] / r**3) * dy
                    az = (-c.G.cgs * m[j] / r**3) * dz
                    ax = ax.value
                    ay = ay.value
                    az = az.value
                    solved_vector[ioffset+3] += ax
                    solved_vector[ioffset+4] += ay
                    solved_vector[ioffset+5] += az    

        return solved_vector 

    def rk4(self, t, dt):
        func = self.gravity_evaluator
        y = self.system_state_vector
        m = self.system_mass_vector
        k1 = dt * func(t, y, m) 
        k2 = dt * func(t + 0.5*dt, y + 0.5*k1, m)
        k3 = dt * func(t + 0.5*dt, y + 0.5*k2, m)
        k4 = dt * func(t + dt, y + k3, m)

        new = y + (1/6.)*(k1+ 2*k2 + 2*k3 + k4)

        return new

    def run_simulation(self, T, dt, t0=0):
        self.total_time = T 
        self.time_step = dt
        clock = (t0 * T.unit).cgs.value
        T = T.cgs.value
        dt = dt.cgs.value
        num_steps = int((T-t0)/dt)
        self.past = []
        print("Running Simulation:")
        for step in tqdm(range(num_steps)):
            self.past.append(self.system_state_vector)
            self.system_state_vector = self.rk4(clock, dt)
            clock += dt

    def transpose_history(self):
        transposed = []
        for y in self.past:
            split = np.split(y, len(y) / 3)
            if len(transposed) == 0:
                transposed = [[[], [], []] for i in range(len(split) // 2)]

            for i in range(0, len(split), 2):
                transposed[i // 2][0].append(split[i][0])
                transposed[i // 2][1].append(split[i][1])
                transposed[i // 2][2].append(split[i][2])

        return transposed


    def plot(self, history, box_length):
        ax = plt.axes(projection='3d')
        ax.set_xlim(-1*box_length, box_length)
        ax.set_ylim(-1*box_length, box_length)
        ax.set_zlim(-1*box_length, box_length)
        ax.set_xlabel('Centimeters')
        ax.set_ylabel('Centimeters')
        ax.set_zlabel('Centimeters')
        ax.view_init(20, 0)
        for val in history:
            ax.plot3D(val[0], val[1], val[2])

        plt.show()
        plt.clf()

    def save_video(self, history, box_length, video_name, rotate=0, fixed='no'):
        fig = plt.figure()
        ax = plt.axes(projection ='3d')

        ax.set_xlabel('Centimeters')
        ax.set_ylabel('Centimeters')
        ax.set_zlabel('Centimeters')
        labels = [i.name for i in self.bodies]
        FFMpegWriter = manimation.writers['ffmpeg']
        metadata = dict(title=f'{video_name}', 
                        artist='Matplotlib',
                        comment='gravity simulation')
        writer = FFMpegWriter(fps=15, metadata=metadata)
        with writer.saving(fig, f"{video_name}.mp4", 100):
            print("Creating Video:")
            for frame_num in tqdm(range(len(self.past))):
                for i, body in enumerate(history):
                    ax.plot3D(body[0][:frame_num], body[1][:frame_num], body[2][:frame_num])
                    ax.scatter(body[0][frame_num], body[1][frame_num], body[2][frame_num])
                    ax.text3D(body[0][:frame_num+1][-1], 
                              body[1][:frame_num+1][-1], 
                              body[2][:frame_num+1][-1], 
                              labels[i], 
                              color='black',
                              size = 15)
                    if fixed == 'yes':
                        ax.set_xlim(-1*box_length, box_length)
                        ax.set_ylim(-1*box_length, box_length)
                        ax.set_zlim(-1*box_length, box_length)
                    ax.set_xlabel('Centimeters')
                    ax.set_ylabel('Centimeters')
                    ax.set_zlabel('Centimeters')
                    if rotate == 'right':
                        coeff = -1
                    if rotate == 'none':
                        coeff = 0
                    if rotate == 'left':
                        coeff = 1
                    ax.view_init(20, coeff*frame_num/3)
                day = round(frame_num * self.time_step.unit.to('yr'), 2)
                ax.set_title(f'{day} years')
                writer.grab_frame()
                plt.cla()

    # def make_slider(self, history):
    #     fig = plt.figure()
    #     ax = plt.axes(projection='3d')
    #     for i in history:
        
    #         graphs = [ax.plot3D([], [], []) for i in range(len(history))]

    #     def animate(i):
    #         for k in range(len(graphs)):
    #             graphs[k].set_data(history[k][0][:i+1], history[k][1][:i+1], history[k][2][:i+1])
        
    #     ani = FuncAnimation(fig, animate, frames=10, interval=200)
    #     plt.show()

        
