#!/usr/bin/env python
# coding: utf-8

# # Introduction: The Planck Function - Jacob Parnell 44597452

# In this lab you will use python to make plots, and compute integrals and derivatives. Each task involves writing some python code, and perhaps producing a plot. Fill in your code below each task (including any plots), creating new 'cells' if you have to. Feel free to include text that explains what you are doing (for example, using headings to indicate your responses, or using <font color=blue>blue text</font> so I can see your contributions clearly), and to 'comment' your code as much as possible to explain what you are doing. I will give you feedback in <font color=red>red text</font> after marking it.

# In[1]:


# This cell is used to load an image of the Planck function from the web:
from IPython.display import Image
i = Image(url='https://dl.dropboxusercontent.com/u/33496045/Planck.png')
i


# In week 2 we shall discuss how a radiation field is described by its intensity $I_\lambda$ or $I_\nu$,  the energy  crossing unit area per unit time per unit solid angle per unit wavelength or frequency interval.
# 
# The Planck function describes the intensity of a black body radiation field.  It can be defined in terms of wavelength:
# 
# $$B_\lambda (T) = \frac{2hc^2}{\lambda^5} \frac{1}{\exp(hc/\lambda kT)-1}$$
# 
# (W m$^{-3}$ sr$^{-1}$) or frequency, $B_\nu$:
# 
# $$B_\nu (T) = \frac{2h\nu^3}{c^2} \frac{1}{\exp(h\nu/kT)-1}$$
# 
# (W m$^{-3}$ sr$^{-1}$).

# In[2]:


