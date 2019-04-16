# coding=utf-8
import numpy as np
import matplotlib.pyplot as plt

import connection_python_server

n_frames = 25002  # SET BY HAND
fig, ax = plt.subplots()
ax.set_xlim([0, n_frames])
ax.set_ylim([0, 1])
point, = ax.plot([0], [0], 'r*', markersize=12)

@connection_python_server.split_byte_stream(int)
def update(frame):
    """
    Update function to move the current point in the graph.

    Parameters
    ----------
    t : int
        frame
    """
    x = int(frame)
    point.set_data([x], [np.sin(x / n_frames * np.pi)])
    point.get_figure().canvas.draw()

try:
    fig.show()
except AttributeError as e:
    print(e)


connection_python_server.run_server(update)
