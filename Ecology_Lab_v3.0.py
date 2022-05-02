########################################################################################################################
#           ECOLOGY LAB V3.0
#           Gerard Solanes
########################################################################################################################
#   Version Notes
#   - Lotka-Volterra equations applied, i.e. no further interaction takes place
#   - GUI is only for visaluzation purposes, but plays no further role.
#   - 2 visualization tools are used: a) Field visualization b) Population time series following L-V equations
########################################################################################################################
#   Import

import numpy as np
import matplotlib.pyplot as plt
import pygame
import random

########################################################################################################################
#   Importing missing modules
# This is a test for commit




########################################################################################################################
#   Ecological Variables Set Up

alpha = 1
beta = 1.
delta = 1.5
gamma = 1.
x0 = 4.
y0 = 2.

########################################################################################################################
#   GUI Set Up
plt.ion()

x_max = 15
y_max = 15
pygame.init()

win = pygame.display.set_mode((500, 500))
pygame.display.set_caption("Ecology Lab")

window_width = 500
window_heigth = 500
width = 15
height = 15
margin_x = (window_width - 20 * x_max) / 2
margin_y = 100

title_text = 'Ecology Lab'
predator_text = 'Predators: '
herbivore_text = 'Herbivores: '
time_text = 'Time: '

font_title = pygame.font.SysFont("menlo", 20, False)
font_body = pygame.font.SysFont("menlo", 15, False)

pred_color = (204, 0, 0)
herb_color = (102, 178, 255)
grass_color = (229, 255, 204)
background_color = (255, 255, 255)


########################################################################################################################
#   Functions

def derivative(X, t, alpha, beta, delta, gamma):
    x, y = X
    dotx = x * (alpha - beta * y)
    doty = y * (-delta + gamma * x)
    return np.array([dotx, doty])

def RK4(func, X0, t, alpha,  beta, delta, gamma):
    """
    Runge Kutta 4 solver.
    """
    dt = t[1] - t[0]
    nt = len(t)
    X  = np.zeros([nt, len(X0)])
    X[0] = X0
    for i in range(nt-1):
        k1 = func(X[i], t[i], alpha,  beta, delta, gamma)
        k2 = func(X[i] + dt/2. * k1, t[i] + dt/2., alpha,  beta, delta, gamma)
        k3 = func(X[i] + dt/2. * k2, t[i] + dt/2., alpha,  beta, delta, gamma)
        k4 = func(X[i] + dt    * k3, t[i] + dt, alpha,  beta, delta, gamma)
        X[i+1] = X[i] + dt / 6. * (k1 + 2. * k2 + 2. * k3 + k4)
    return X

def plotting(predators, herbivores, time):
    plt.clf()
    plt.title("Time Series")
    plt.plot(time, herbivores, color='b', label='Herbivores')
    plt.plot(time, predators, color='r', label="Predators")
    plt.grid()
    plt.xlabel("Time")
    plt.ylabel('Population')
    plt.legend(loc="best")

    plt.draw()
    plt.pause(0.01)


def draw_GUI(animals_mat, num_predators, num_herbivores, time_step):
    for x in range(len(animals_mat)):
        for y in range(len(animals_mat[0])):
            if animals_mat[x][y] == 0:
                pygame.draw.rect(win, grass_color, (margin_x + x * 20, margin_y + y * 20, width, height))
            if animals_mat[x][y] == 1:
                pygame.draw.rect(win, herb_color, (margin_x + x * 20, margin_y + y * 20, width, height))
            elif animals_mat[x][y] == 2:
                pygame.draw.rect(win, pred_color, (margin_x + x * 20, margin_y + y * 20, width, height))

    text = font_title.render(title_text, 1, (0, 0, 0))
    win.blit(text, (1 / 3 * window_width, 10))
    text = font_body.render(time_text + str(time_step), 1, (0, 0, 0))
    win.blit(text, (20, 50))
    text = font_body.render(predator_text + str(num_predators), 1, (0, 0, 0))
    win.blit(text, (20, margin_y + y * 20 + 50))
    text = font_body.render(herbivore_text + str(num_herbivores), 1, (0, 0, 0))
    win.blit(text, (1 / 3 * window_width, margin_y + y * 20 + 50))

    pygame.display.update()


########################################################################################################################
#   Solving equations

dt = 500
tmax = 30.
t = np.linspace(0.,tmax, dt)
X0 = [x0, y0]

Xrk4 = RK4(derivative, X0, t, alpha,  beta, delta, gamma)

########################################################################################################################
#   GUI

###########################################################
# Representation matrix

#   Set Up Field
field = [[0 for i in range(x_max)] for j in range(y_max)]
predators = Xrk4[:, 1]
predators_plot = []
predators_theo = []

herbivores = Xrk4[:, 0]
herbivores_plot = []
herbivores_theo = []

#Adapting theoretical results to plot (*10)
for i in range(len(predators)):
    predators_theo.append(predators[i]*10)
    herbivores_theo.append(herbivores[i]*10)

time = t
time_plot = []

for clock in range(len(time)):
    field = [[0 for i in range(x_max)] for j in range(y_max)]
    num_predators_0 = int(predators[clock]*10)
    num_herbivores_0 = int(herbivores[clock]*10)

    num_predators = 0
    num_herbivores = 0
    while num_predators_0 != 0 and num_herbivores_0 != 0:
        for x in range(len(field)):
            for y in range(len(field[0])):
                p_h = random.randint(1, 100)
                p_p = random.randint(1, 100)
                if p_h > 90 and num_herbivores_0 > 0 and field[x][y] == 0:
                    field[x][y] = 1
                    num_herbivores_0 -= 1
                    num_herbivores += 1
                if p_p > 90 and num_predators_0 > 0 and field[x][y] == 0:
                    field[x][y] = 2
                    num_predators_0 -= 1
                    num_predators += 1

    ######################################################
    # Animated plot
    predators_plot.append(num_predators)
    herbivores_plot.append(num_herbivores)
    time_plot.append(clock)
    plotting(predators_plot, herbivores_plot, time_plot)

    ######################################################
    # Draw GUI
    draw_GUI(field, num_predators, num_herbivores, clock)
    pygame.time.delay(1)
    win.fill(background_color)

    for event in pygame.event.get():
        if event.type == pygame.QUIT:
            run = False

###########################################################
# Final plot
plt.close()
plt.ioff()
plt.title("Time Series")
plt.plot(time, herbivores_plot, color='b', label='Herbivores')
plt.plot(time, predators_plot, color='r', label="Predators")
plt.plot(time, herbivores_theo, color='c', label='Theoretical Herbivores')
plt.plot(time, predators_theo, color='y', label="Theoretical Predators")
plt.grid()
plt.xlabel("Time")
plt.ylabel('Population')
plt.legend(loc="best")

plt.show()
#plotting(Xrk4[:, 0], Xrk4[:, 1], t)
