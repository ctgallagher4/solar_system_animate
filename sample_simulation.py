from simulation_module import *

simulation_system = System(path.join(sys.path[0], 'sample_binary_star.csv'))

simulation_system.run_simulation(4*u.year, 1*u.day)
history = simulation_system.transpose_history()
simulation_system.plot(history, (1*u.AU).cgs.value)
simulation_system.save_video(history,
                             (1*u.AU).cgs.value, 
                             'Earth-Moon', 
                              rotate=1/3,
                              fixed='yes',
                              frame_skip=10)