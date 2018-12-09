import cgi
import cgitb
cgitb.enable()


import numpy as np
import matplotlib.pyplot as plt
import mpld3
from mpld3 import plugins


print("Content-type: text/html")

# chmod 755 script.py

template = "<html><body><h1>Hello {}!</h1></body></html>"


print(template.format("Reader"))


form = cgi.FieldStorage()

# plot line + confidence interval
fig, ax = plt.subplots()
ax.grid(True, alpha=0.3)

ax.plot([1,2,3,4,5,6], [2,3,4,5,3,4], 'o')

# define interactive legend

# handles, labels = ax.get_legend_handles_labels() # return lines and labels
# interactive_legend = plugins.InteractiveLegendPlugin(zip(handles,
#                                                          ax.collections),
#                                                      labels,
#                                                      alpha_unsel=0.5,
#                                                      alpha_over=1.5,
#                                                      start_visible=True)
# plugins.connect(fig, interactive_legend)

ax.set_xlabel('x')
ax.set_ylabel('y')
ax.set_title('Interactive legend', size=20)

mpld3.show()



# try:   # NEW
#     print("Content-type: text/html\n\n")   # say generating html
#     main()
# except:
#     cgi.print_exception()                 # catch and print errors




