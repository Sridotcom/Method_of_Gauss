# Summer Science Program - Srivishal Sudharsan
Hello! I wrote this code during my time at the Summer Science Program in Astrophysics at New Mexico State University. For more information about the technical details please refer to the [Orbit Report](https://github.com/Sridotcom/Method_of_Gauss/blob/7fb87df81e02f0775e7f534a1a8cdf380b190f9e/Copy%20of%20SSP%20Orbit%20Report%20(2).pdf). In short, my group and I utilized an orbit determination technique known as the Method of Gauss to determine the orbital mechanics of our chosen asteroid given 3 observations and the position of our asteroid at those observations. 

## What is the Method of Gauss?
The Method of Gauss is a mathematical technique used in astronomy to determine the orbits of celestial objects like asteroids and planets. It relies on  observations of the object's position in the sky at different times and uses various mathematical techniques to calculate the object's orbital parameters, such as its path shape and how it moves around a central body, like the Sun. By fitting observations to these mathematical rules, it helps astronomers predict where the object will be in the future. This method, named after mathematician Carl Friedrich Gauss, is crucial for understanding the motion of celestial bodies.

## Why is this research important?
We observed asteroid 1998 MX5, which is classified as a near-earth asteroid. A Near-Earth Asteroid (NEA) is an asteroid that comes in close proximity to Earth's orbit as it travels around the Sun. Specifically, NEAs are asteroids whose closest approach to Earth is within about 1.3 astronomical units (AU) or less, where 1 AU is the average distance between the Earth and the Sun. The study of NEAs is vital for planetary defense as finding when these asteroids could hit the earth will be very useful for planetary defense. Additionally, finding where the asteroid will be at a certain time will allow us to find optimal rocket launch windows, utilize their gravitational pull for space exploration, and potentially mine them for resources. To determine the possible threat of our asteroid, our group utilized Jupyter Notebook and a Python script to predict the orbit of  "clones" of our asteroid. Random chaos will cause our clones to all have slightly varying orbits and analyzing these orbits will give us insight into possible future collisions. The results of our predictive model are shown below. 

![image2](https://github.com/Sridotcom/Summer-Science-Program-Research/assets/66920443/4c993d11-8e0c-49fa-a66a-e6a5337f9c20)













This code takes in an input of 3 observations of an asteroid and the position of the asteroids at these times. Utilizing the Method of Gauss (outlined in the Copy of Orbit Report file) the program will calculate the orbital parameters of the asteroid (ie how the orbit is shaped). 

![0](https://github.com/Sridotcom/Method_of_Gauss/assets/66920443/0791a971-4538-49bb-ab6d-7e1dbc043c73)
