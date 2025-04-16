import os

def replace_string_in_file(file_path, old_string, new_string):
    """Replace old_string with new_string in the given file."""
    try:
        with open(file_path, 'r') as file:
            content = file.read()
        
        # Replace the target string
        content = content.replace(old_string, new_string)
        
        with open(file_path, 'w') as file:
            file.write(content)
        print(f"Replaced in file: {file_path}")
    except Exception as e:
        print(f"Error processing file {file_path}: {e}")

def replace_string_in_directory(directory, old_string, new_string):
    """Recursively replace old_string with new_string in all files in the given directory and subdirectories."""
    for dirpath, _, filenames in os.walk(directory):
        for filename in filenames:
            file_path = os.path.join(dirpath, filename)
            replace_string_in_file(file_path, old_string, new_string)

if __name__ == "__main__":
    directory_to_search = './'  # Change this to your directory
    old_string = 'CMSSW_14_0_0_pre0/src'
    new_string = 'CMSSW_14_0_0_pre0/src'
    
    replace_string_in_directory(directory_to_search, old_string, new_string)
