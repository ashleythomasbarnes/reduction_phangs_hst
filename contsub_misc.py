import os 
import re


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
        print(f"Directory '{directory}' created successfully.")
    else:
        print(f"Directory '{directory}' already exists.")
