# Proposals

Project idea #1: 
I would like to build a 3D n-body orbit simulator in python using runge-kutta methods and matplotlib. It would be able to take input for locations, masses, and velocities of objects from a csv spreadsheet and then plot the orbits with a matplotlib animation. One example use would be modeling a solar system. I think this would be really interesting and give me experience with Runge-Kutta methods in python which would be very helpful for my future career aspirations should I choose to go into GNC software development.

Project idea #2:

I am the CTO of a startup called ZERO for men (www.zeroformen.com). We are currently entered in the College New Venture Challenge (https://polsky.uchicago.edu/programs-events/new-venture-challenge/college-new-venture-challenge/). We would like to run a micro-influencer marketing campaign at some point, and the cheapest way to do this is to find a lot of small youtube channels that fit our brand and send them products for free. These youtube channels then review the product and it is a win-win for everyone. Unfortunately, finding youtube channels that are large enough to drive conversions, but small enough that they do not expect us to pay for a video is challenging and requires a lot of clicking. Obviously, larger channels are prioritized and so these “micro-influencers” can get buried. Solutions certainly exist (https://dovetale.com), but the quality solutions cost money. We do not have money. I would like to make a youtube scraper that will provide me with a list of channels and their contact emails if the channel fits a list of requirements. For example, that they have done a product review before, or that they have a subscriber count around 10,000 to 20,000.

Project idea #3 (Optional):

I keep myself organized with a running to-do list on a piece of paper. This requires wasting a lot of paper. Each day I copy over anything that was not completed from the day before on to the next day’s to do list. Last quarter in Immersion Programming, we made a task manager in python. It is ugly and requires a lot of typing. I would like to make a task manager with a gui in python.

# To-Do List (will be updated as I go)

* Implement RK-4 methods in python or C
* Interface with C if necessary for speed
* Implement celestial body object classes
* Implement a collision model
* Implement matplotlib animations
* Implment read from csv functionality


# Timeline for final project idea #1

* By the end of week 4, I will have read up on implementing Runge-Kutta methods in Python and C.

* By the end of week 5, I will have implemented object classes for celestial bodies which contain attributes such as position, velocity, and mass and the ability to read these from a csv file.

* By the end of week 6, I will have decided what to do about collisions. Obviously this is a complex problem physically so I will probably just have two collision types like "sticky" and "non-sticky". From there I will need to implement some sort of momentum transfer. This will take some research and some remembering basic pool table physics.

* By the end of week 7, I will have finished implementing Runge-Kutta methods.

* By the end of week 8, I will have a matplotlib animation of the entire system.

Finally, I will structure all of these files in a good way and clean things up.








