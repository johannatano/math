import os
import sys

def get_project_root() -> str:
    """
    Returns the project root by going up one level from the script's location,
    or from the current working directory if running in a notebook.
    """
    if hasattr(sys, "_getframe"):
        try:
            script_dir = os.path.dirname(os.path.abspath(__file__))
        except NameError:
            script_dir = os.getcwd()  # fallback for notebooks
    else:
        script_dir = os.getcwd()

    # Go up TWO levels to get project root
    return os.path.dirname(os.path.dirname(script_dir))
    return os.path.dirname(script_dir)  # go up one level to project root

class Path:
    @staticmethod
    def exports(name: str) -> str:
        """
        Returns the full path to a file inside the 'exports' folder,
        relative to the project root.
        """
        return os.path.join(get_project_root(), "exports", name)
