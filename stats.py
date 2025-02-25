#First install the following libraries
#pip install pyautogui
#pip install pynput

#to run type: python ./stats.py

import pyautogui
import threading
import datetime
from pynput import mouse

screenSize = pyautogui.size()
manual_movement_detected = False  # Flag to track manual movement

def on_move(x, y):
    """Callback for detecting manual mouse movement."""
    global manual_movement_detected
    manual_movement_detected = True
    # print("Manual movement detected! Exiting...")
    quit()

def moveMouse():
    """Moves the mouse to a corner of the screen."""
    if manual_movement_detected:
        quit()
    # pyautogui.moveTo(5, screenSize[1], duration=1)

def clickMouse():
    """Clicks the mouse."""
    if manual_movement_detected:
        quit()
    # print("click")

    pyautogui.click()
    main()

def main():
    """Main function to control mouse movement and clicking."""
    hour = datetime.datetime.now().hour
    if hour == 17 or hour == 12:
        print("End of day reached")
        quit()
    else:
        threading.Timer(5.0, moveMouse).start()
        threading.Timer(10.0, clickMouse).start()

# Start mouse listener in a separate thread
listener = mouse.Listener(on_move=on_move)
listener.start()

main()
