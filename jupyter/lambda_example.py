
# coding: utf-8

# In[1]:

import sympy as sp
import scipy.optimize as opt


# In[15]:

x, y, z, a, b, c, d, e, g, h, i, j, k = sp.symbols('x y z a b c d e g h i j k')


# In[16]:

f = x + y + z + a +b+c+d+e+g+h+i+j+k


# In[7]:

f2 = sp.lambdify((x,y,z,a, b, c, d, e, g, h,i,j,k),f,'numpy')


# In[8]:

f2


# In[13]:




# In[10]:

def f3(x,y,z,a,b):
    return lambda c,d,e,g,h,i,j,k: f2(x,y,z,a, b, c, d, e, g, h,i,j,k)


# In[11]:

f3(11,1,1,1,1)(1,1,1,1,1,2,2,2)


# In[12]:

f2(1,1,1,1,1,1,1,1,1,1,1,1,1)


# In[ ]:



