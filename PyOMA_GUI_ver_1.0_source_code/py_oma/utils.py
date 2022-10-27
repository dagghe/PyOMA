import sys
import os


def resolve_path(path):
    if getattr(sys, "frozen", False):
        # If the 'frozen' flag is set, we are in bundled-app mode!
        resolved_path = os.path.abspath(os.path.join(sys._MEIPASS, path))
    else:
        # Typical mode
        ui_path = os.path.dirname(os.path.abspath(__file__))
        resolved_path = os.path.join(ui_path, path)

    return resolved_path
