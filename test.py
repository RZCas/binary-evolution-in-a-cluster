from matplotlib import pyplot as plt

fig, (ax1, ax2) = plt.subplots(2)

#standard is that both x- and y-axis are autoscaled
ax1.plot([1, 3, 5], [2, 5, 1], label="autoscale on")
#rendering the current output
fig.draw_without_rendering() 
#turning off autoscale for the x-axis of the upper panel
#the y-axis will still be autoscaled for all following artists
ax1.autoscale(False, axis="x")
ax1.plot([-1, 7], [-2, 4], label="autoscale off")
ax1.legend()

#other axis objects are not influenced
ax2.plot([-2, 4], [3, 1])

plt.show()
