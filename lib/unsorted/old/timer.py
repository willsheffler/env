import sys,os
from Tkinter import *
from time import sleep

t = 300
if len(sys.argv) > 1:
    t = int(sys.argv[1]) 
message = "WAKE UP!!!"
if len(sys.argv) > 2:
    message = sys.argv[2]

class App:

    def __init__(self, master):

        frame = Frame(master)
        frame.pack()

        self.button = Button(frame, text="QUIT", fg="red", command=frame.quit)
        self.button.pack(side=LEFT)

        self.hi_there = Button(frame, text=message,
                               font=("Arial", 64 , "bold"))
        self.hi_there.pack(side=LEFT)


root = Tk()

app = App(root)


sleep(t)

root.mainloop()

