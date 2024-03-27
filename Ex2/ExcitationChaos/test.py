from matplotlib import pyplot as plt
from matplotlib.patches import Rectangle
from matplotlib.patches import ConnectionPatch


# Generate some example data
x = range(100)
y = [i ** 2 for i in x]

# Create the main plot
fig, ax = plt.subplots()
ax.plot(x, y, label='Main Plot')

# Define the zoom region
xlim = [20, 40]
ylim = [300, 1000]

# Create inset plot
axin = ax.inset_axes([0.6, 0.6, 0.35, 0.35])  # [x, y, width, height]
axin.plot(x, y)
axin.set_xlim(xlim[0], xlim[1])
axin.set_ylim(ylim[0], ylim[1])
axin.set_title('Zoomed Inset')

axin.tick_params(axis='both', which='both', bottom=False, top=False, left=False, right=False, labelbottom=False, labelleft=False)



rect = Rectangle((xlim[0], ylim[0]), xlim[1] - xlim[0], ylim[1] - ylim[0], 
                 linewidth=1, edgecolor='gray', facecolor='none')
ax.add_patch(rect)

# Calculate the center of the rectangle
center_x = (xlim[0] + xlim[1]) / 2
center_y = (ylim[0] + ylim[1]) / 2

# Add an arrow annotation from the center of the rectangle to the zoomed plot
arrow_props = dict(arrowstyle="-", color='black',zorder=10)
arrow = ConnectionPatch((xlim[1], ylim[0]), (99, 6000), "data", "data", **arrow_props)
ax.add_artist(arrow)

plt.savefig('inset_plot.png')