

x_points = np.linspace(-2,2,100)
y_points = np.linspace(-2,2,100)
X, Y = np.meshgrid(x_points, y_points)   
Z= twoD_Gaussian(X, Y, amplitude=-1.0, xo=0.0, yo=0.0, sigma_x=1.0, sigma_y=1.0, theta=0, offset=0)
CS = plt.contour(X, Y, Z, levels=10)
plt.clabel(CS, inline=1, fontsize=10)
