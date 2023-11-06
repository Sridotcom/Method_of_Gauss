import math as m; import numpy as np

# Constants defined
k = 0.0172020989484 # defined by the IAU
ob = m.radians(23.43829194) # obliquity/tilt of earth (As of July 15 2023)
c = 173.145 # speed of light in au/day
E = m.radians(23.43829194)

data = []
sun_vectors = []
file = open('elements.txt','r')
lines = file.read().split("\n")
for line in lines: 
    raw_data = line.split(" ")
    data.append(raw_data)

# Returns orbital elements given position and velocity vectors
def baby_od (position_x, position_y, position_z, velocity_vector_x, velocity_vector_y, velocity_vector_z):
    # Stores inputs into numpy arrays
    position_vector = np.array([position_x, position_y, position_z])
    velocity_vector = np.array([velocity_vector_x, velocity_vector_y, velocity_vector_z])

    position_magnitude = m.sqrt((position_x**2) + (position_y**2) + position_z**2)

    # a represents semi major axis of the asteroids elliptical orbit
    a = ((2/position_magnitude)- np.dot(velocity_vector, velocity_vector))**-1

    # n represents the mean motion in degrees/day which shows how much it moves along the ellipse every day
    n = m.degrees((k/(m.sqrt(a**3))))

    position_cross_velocity = np.cross(position_vector, velocity_vector)

    # e represents eccentricity which shows how stretched out the orbit is                   
    e = m.sqrt (1- ((position_cross_velocity[0]**2) + (position_cross_velocity[1]**2)+ (position_cross_velocity[2]**2))/a)

    position_cross_velocity = np.cross(position_vector, velocity_vector)

    z = position_cross_velocity [2]
    cross_magnitutde = m.sqrt((position_cross_velocity[0]**2) + (position_cross_velocity[1]**2) + (position_cross_velocity[2]**2)) 

    # I is inclintation which is the angle between the asteroids orbital plane and the earths orbital plane (ecliptic)
    I = m.degrees((m.acos((z)/ (cross_magnitutde))))
    h = position_cross_velocity

    hx = h[0]
    hy = h[1]

    h = m.sqrt(a*(1-e**2))

    sin_omega = ((hx)/(h*m.sin(m.radians(I))))
    cos_omega = ((-hy)/(h*m.sin(m.radians(I))))


    omega = m.degrees(m.atan2(sin_omega, cos_omega))


    X = position_vector[0]
    Y = position_vector[1]
    Z = position_vector[2]
    r = position_magnitude

    sin_f_w = Z/(r*m.sin(m.radians(I)))

    cos_f_w = (1/m.cos(m.radians(omega)))*((X/r) + m.cos(m.radians(I))*sin_f_w*sin_omega)

    f_w = m.degrees(m.atan2(sin_f_w, cos_f_w))

    cos_f = (1/e)*(((a*(1-e**2))/(r))-1)

    sin_f = (np.dot(position_vector, velocity_vector)/(e*r))*(m.sqrt(a*(1-e**2)))

    f =  m.degrees(m.atan2(sin_f,cos_f))


    w = f_w - f
    cos_E = (1/e)*(1- (r/a))


    E = (2*m.pi) - m.acos(cos_E)
    f = 360 - f
    if f<0:
        E = -abs(E)
    else:
        E = abs(E)


    M_m = E - e*m.sin(E)

    print("a: ", a, "au")
    print("e: ", e)
    print("n:", n, "deg/day")
    print("i: ", I, "degrees")
    print("Ω: ", omega, "degrees")
    print('ω: ',  w, "degrees")
    print("M: ", m.degrees(M_m), "degrees")
    return(a, e, n, I, omega, w, m.degrees(M_m))


# Function returns f1, f3, g1, and g3 values and can be changed for 3rd, 4th, or closed function based on the flag
def f_g(tau1, tau3, r2,r2dot,flag):
    mew = 1
    r2_mag = m.sqrt(r2[0]**2 + r2[1]**2 + r2[2]**2)
    u = mew/(r2_mag**3)
    z = (np.dot(r2, r2dot))/(r2_mag**2)
    q = ((np.dot(r2dot,r2dot))/(r2_mag**2)) - u
    if flag == "3rd":
        f1 = 1 - (1/2)*(u)*(tau1**2) + (1/2)*(u*z)*(tau1**3) 
        f3 = 1 - (1/2)*(u)*(tau3**2) + (1/2)*(u*z)*(tau3**3)
        g1 = tau1 - (1/6)*(u)*(tau1**3)
        g3 = tau3 - (1/6)*(u)*(tau3**3)
    if flag == "4th":
        f1 = 1 - (1/2)*(u)*(tau1**2) + (1/2)*(u*z)*(tau1**3) + (1/24)*(3*u*q - 15*u*(z**2) + u**2)*(tau1**4)
        f3 = 1 - (1/2)*(u)*(tau3**2) + (1/2)*(u*z)*(tau3**3) + (1/24)*(3*u*q - 15*u*(z**2) + u**2)*(tau3**4)
        g1 = tau1 - (1/6)*(u)*(tau1**3) + (1/4)*(u*z*tau1**4)
        g3 = tau3 - (1/6)*(u)*(tau3**3) + (1/4)*(u*z*tau3**4)

    if flag == "function":
        r2_mag = m.sqrt(r2[0]**2 + r2[1]**2 + r2[2]**2)

        a = baby_od(r2[0], r2[1], r2[2], r2dot[0], r2dot[1], r2dot[2])[0]
        n = baby_od(r2[0], r2[1], r2[2], r2dot[0], r2dot[1], r2dot[2])[1]
        e = baby_od(r2[0], r2[1], r2[2], r2dot[0], r2dot[1], r2dot[2])[2]


        chunk = (np.dot(r2, r2dot))/((n*(a**2)))
        def f_x (x,tau):
            return (x - (1-(r2_mag/a)) * m.sin(x) + chunk*(1-m.cos(x)) - n*tau)

        def f_prime (x):
            return(1- (1- (r2_mag/a))*m.cos(x) + chunk*m.sin(x))
        def x_eq (tau, sign):
            if sign == 1:
                return ((n*tau + 0.85*e) - chunk)
            if sign == -1:
                return ((n*tau - 0.85*e) - chunk)
        def sign_func (tau):
            return ((((chunk))*m.cos(n*tau - chunk)) + (1- (r2_mag/a)*(m.sin(n*tau - chunk))))

        sign1 = np.sign(sign_func(tau1))
        sign3 = np.sign(sign_func(tau3))


        x1_init = x_eq(tau1, sign1)
        x3_init = x_eq(tau3, sign3)

        # Loop for tau 1
        x_old_1 = x1_init + 1
        x_old_3 = x3_init + 1

        while abs(x_old_1 - x1_init) > 1e-12:
            x_old_1 = x1_init
            x1_init = x1_init - ((f_x(x1_init, tau1))/f_prime(x1_init))


        while abs(x_old_3 - x3_init) >  1e-12:
            x_old_3 = x3_init
            x3_init = x3_init - ((f_x(x3_init, tau3))/f_prime(x3_init))
        
        delta_e1 = x1_init
        delta_e3 = x3_init

        f1 = 1- (a/r2_mag)*(1-m.cos(delta_e1))
        f3 = 1- (a/r2_mag)*(1-m.cos(delta_e3))
        g1 = tau1 + (1/n)*(m.sin(delta_e1) -delta_e1)
        g3 = tau3 + (1/n)*(m.sin(delta_e3) -delta_e3)


    return [f1, f3, g1, g3]
# Function used to convert RA and Dec into radians
def convertAngle(degrees, minutes, seconds, radians, normalize ):
    minutes = m.copysign(minutes, degrees) # Takes sign of 2nd argument and moves to 1st argument
    seconds = m.copysign(seconds,degrees)
    minutes = minutes/60
    seconds = seconds/3600
    decimal = (degrees + minutes + seconds)
    if normalize and radians:
        return ((decimal%360) * (m.pi/180))
    if normalize == False and radians == True:
        return (decimal * (m.pi/180))
    if normalize == True and radians == False:
        return ((decimal%360))
    if normalize == False and radians == False:
        return(decimal)

# Converts all dec and ra into radians, and date/time to Julian date
for index in range(len(data)):
    # Converting dec into integers then into radians
    Dec = data[index][5]
    Dec = Dec.split(":")
    degrees = int(Dec[0])
    arcminutes = int(Dec[1])
    arcseconds = float(Dec[2])
    data[index][5] = convertAngle(degrees, arcminutes, arcseconds, True, False)

    # Converting RA values into int/floats then into radian decimals
    RA = data[index][4]
    RA = RA.split(":")
    hours = int(RA[0])
    minutes = int(RA[1])
    seconds = float(RA[2])
    data[index][4] = convertAngle(hours, minutes, seconds, False, False)
    data[index][4] = m.radians(data[index][4]*15)
    data[index][6] = float(data[index][6])
    data[index][7] = float(data[index][7])
    data[index][8] = float(data[index][8])
    data[index][0] = int(data[index][0])
    data[index][1] = int(data[index][1])
    data[index][2] = int(data[index][2])

    # Changes time to decimal time then converts everything into Julian dates
    Y = data[index][0]
    M = data[index][1]
    D = data[index][2]
    time = data[index][3]
    time = time.split(":")
    t_hours = int(time[0])
    t_minutes = int(time[1])
    t_seconds = int(time[2])
    data[index][3] = convertAngle(t_hours, t_minutes, t_seconds, False, False)
    J_0 = 367*(Y) - int((7/4)*(Y + int((M+9)/12))) + int((275*M)/9) + D + 1721013.5
    JD = J_0 + ((data[index][3])/24)
    data[index].append(JD)
    

t1 = data[0][9]   
t2 = data[1][9]
t3 = data[2][9] 

# Converting proper time to gaussian time
tau3 = k*(t3 - t2)
tau1 = k*(t1 - t2)
tau = tau3 - tau1

# Solving for fixed row vectors which represents direction from earth to the asteroid
row_vectors = []
for index in range(len(data)):
    dec = data[index][5]
    ra = data[index][4]
    row_hat = np.array([m.cos(ra)*m.cos(dec), m.sin(ra)*m.cos(dec), m.sin(dec)])
    row_vectors.append(row_hat)


# Converts row hat vectors into equatorial coordinates because asteroid is orbiting the sun
for index in range(len(row_vectors)):
    row_vectors[index] = np.matmul(np.array([[1, 0, 0], [0, m.cos(ob), m.sin(ob)], [0, -m.sin(ob), m.cos(ob)]]), row_vectors[index])

# Sets row hat vectors
phat_1 = row_vectors[0]
phat_2 = row_vectors[1]
phat_3 = row_vectors[2]


# Converting Sun Vector into equatorial
sun_arrays  = []
for index in range(len(data)):
    sun_array = np.array([data[index][6] , data[index][7] , data[index][8]])
    sun_arrays.append(sun_array)
for index in range(len(sun_arrays)):
    sun_arrays[index] = np.matmul(np.array([[1, 0, 0], [0, m.cos(ob), m.sin(ob)], [0, -m.sin(ob), m.cos(ob)]]), sun_arrays[index])


R1 = sun_arrays[0]
R2 = sun_arrays[1]
R3 = sun_arrays[2]


# Calculating D values (stay the same throughout) (checked)
D0 = np.dot(phat_1, np.cross(phat_2, phat_3))
D11 = np.dot(np.cross(R1, phat_2), phat_3)
D12 = np.dot(np.cross(R2, phat_2), phat_3)
D13 = np.dot(np.cross(R3, phat_2), phat_3)
D21 = np.dot(np.cross(phat_1, R1), phat_3)
D22 = np.dot(np.cross(phat_1, R2), phat_3)
D23 = np.dot(np.cross(phat_1, R3), phat_3)
D31 = np.dot(phat_1, np.cross(phat_2, R1))
D32 = np.dot(phat_1, np.cross(phat_2, R2))
D33 = np.dot(phat_1, np.cross(phat_2, R3))

# Calculate initial r value via Keplers 2nd law
# Using Keplers Laws to solve for initial r2 magnitude
# Calculating c values
c1 = tau3/tau
c2 = -1
c3 = -tau1/tau

# Using Equation 103 solve for p scalar values (checked)
p1 = (c1*D11 + c2*D12 + c3*D13)/(c1*D0)
p2 = (c1*D21 + c2*D22 + c3*D23)/(c2*D0)
p3 = (c1*D31 + c2*D32 + c3*D33)/(c3*D0)

# Solving for r vector values utilizing the fundamental triangle
r1_vec = p1*phat_1 - R1
r2_vec = p2*phat_2 - R2
r3_vec = p3*phat_3 - R3

# Finding r2dot via equation 131
r2dot = (r3_vec-r1_vec)/tau

# Getting r2 magnitude value
r2mag = np.linalg.norm(r2_vec)


p2old = p2 + 1

n = 0


# Iterating through to converge on a r vector value
while abs(p2old-p2) > 1e-12 and n<1000:
    p2old = p2
    # Correcting for light travel time
    t1_new = t1 - (p1/c)
    t2_new = t2 - (p2/c)
    t3_new = t3 - (p3/c)
    # Recalcuating tau values using new time values
    tau3 = k*(t3_new-t2_new)
    tau1 = k*(t1_new-t2_new)
    tau = tau3 - tau1

    # Running through f and g functions to get f1, f3, g1, and g3 values using 4th order taylor series
    f_g_elements = ((f_g(tau1, tau3, r2_vec, r2dot, "4th")))
    f1 = f_g_elements[0]
    f3 = f_g_elements[1]
    g1 = f_g_elements[2]
    g3 = f_g_elements[3]
    


    # Calculating new c and d values via equations 98, 99, and 100
    c1 = g3/(f1*g3 - g1*f3)
    c3 = (-g1/(f1*g3-g1*f3))
    d1 = -f3/(f1*g3 - f3*g1)    
    d3 = f1/(f1*g3 - f3*g1)

    # Check D values
    p1 = (c1*D11 + c2*D12 + c3*D13)/(c1*D0)
    p2 = (c1*D21 + c2*D22 + c3*D23)/(c2*D0)
    p3 = (c1*D31 + c2*D32 + c3*D33)/(c3*D0)

    r1_vec = p1*phat_1 - R1
    r2_vec = p2*phat_2 - R2
    r3_vec = p3*phat_3 - R3

    r2mag = np.linalg.norm(r2_vec) 

    r2dot = d1*r1_vec + d3*r3_vec
    r2dotmag = np.linalg.norm(r2dot) 

    n+= 1
   


print(r2_vec, "=", r2mag, "au")
print(k*r2dot, "=", r2dotmag, "au/day")
print(p2, "au")

print(baby_od(r2_vec[0], r2_vec[1], r2_vec[2], r2dot[0], r2dot[1], r2dot[2]))





