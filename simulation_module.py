from astropy.coordinates import SkyCoord
import csv
from os import path
import sys
import numpy as np
import matplotlib.pyplot as plt
import astropy.units as u 
import astropy.constants as c
import matplotlib.animation as manimation
from tqdm import tqdm
import mpl_toolkits.mplot3d.axes3d as p3
from ctypes import c_double, c_int, CDLL

library_path = path.join(sys.path[0], 'gravity_evaluator.so')
try:
    gravity_evaluator = CDLL(library_path)
except:
    exit("Gravity evaluator module not found.")  
grav_eval = gravity_evaluator.gravity_evaluator
grav_eval.restype = None

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

    def evaluator(self, t, y, m):
        n = len(y)
        y = list(y)
        m = list(m)
        c_arr_in = (c_double * n)(*y)
        c_arr_in_2 = (c_double * n)(*m)
        c_arr_out = (c_double * n)()
        grav_eval(n, c_arr_in, c_arr_in_2, c_arr_out)
        return np.asarray(list(c_arr_out[:]))

    def rk4(self, t, dt):
        func = self.evaluator
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

    def save_video(self, history, box_length, video_name, rotate=0, fixed='no', frame_skip=1):
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
            for frame_num in tqdm(range(0, len(self.past), frame_skip)):
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
                ax.view_init(20, rotate*frame_num* 1/frame_skip)
                time = round(frame_num * self.time_step.unit.to('yr'), 2)
                ax.set_title(f'{time} years')
                writer.grab_frame()
                plt.cla()
