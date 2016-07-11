import numpy as np

class rot_3D:
    """Rotation"""

    a = np.array([1, 1, 1])
    b = np.array([1, 1, 1])

    def __init__(self,a,b):
        self.a = a
        self.b = b

    def norm(k):
        np.linalg.norm(k)


    a_x = a_y = a_z = a
    b_x = b_y = b_z = b
    a_x[0] = 0 ; a_y[1] = 0 ; a_z[2] = 0
    b_x[0] = 0 ; b_y[1] = 0 ; b_z[2] = 0

    cos_x = a_x.dot(b_x) / (norm(a_x) * norm(b_x))
    sin_x = (1 - cos_x**2)**(1/2)
    if(abs(cos_x - 1) < 10e-3): sin_x = 0

    cos_y = a_y.dot(b_y) / (norm(a_y) * norm(b_y))
    sin_y = (1 - cos_y**2)**(1/2)
    if(abs(cos_y - 1) < 10e-3): sin_y = 0

    cos_z = a_z.dot(b_z) / (norm(a_z) * norm(b_z))
    sin_z = (1 - cos_z**2)**(1/2)
    if(abs(cos_z - 1) < 10e-3): sin_z = 0

    rot_x = np.array([
        [1,0,0],
        [0,cos_x,-sin_x],
        [0,sin_x,cos_x]
        ])

    rot_y = np.array([
        [cos_y,0,sin_y],
        [0,1,0],
        [-sin_y,0,cos_y]
        ])

    rot_z = np.array([
        [cos_z,-sin_z,0],
        [sin_z,cos_z,0],
        [0,0,1]])


    def rotation(self,P): (self.rot_x * self.rot_y * self.rot_z * P.T).T