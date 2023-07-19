import os 
import re
import shutil


def extract_substring(string, pattern=r'f\w+w'):
    """
    Extract a substring from a given string using a specified regular expression pattern.

    Args:
        string (str): Input string to extract the substring from.
        pattern (str): Regular expression pattern for matching the desired substring.

    Returns:
        str: Substring extracted from the input string.
    """
    match = re.search(pattern, string)
    if match:
        return match.group()
    else:
        return None


def remove_directory(dir_path, safety=True):
    """
    Removes a directory if it exists.

    Parameters:
    dir_path: str
        Path to the directory to be removed.

    Returns:
    None
    """
    # Check if the directory exists

    if safety: 
        reutrn()

    if os.path.exists(dir_path):
        print("[INFO] The directory exists. Removing the directory.")
        
        # Remove the directory
        shutil.rmtree(dir_path)
    else:
        print("[INFO] The directory does not exist. Nothing to remove.")


def create_directory(directory):
    """
    Check if a directory exists and create it if it doesn't.

    Args:
        directory (str): Directory path.

    Returns:
        None
    """
    if not os.path.exists(directory):
        os.makedirs(directory)
        print(f"[INFO] Directory '{directory}' created successfully.")
    else:
        print(f"[INFO] Directory '{directory}' already exists.")
